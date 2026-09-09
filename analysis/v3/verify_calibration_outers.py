"""Restore independent outer metrics from pinned synthetic arrays, without refits."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .module_ecology import summarize_draws
from .protected_artifacts import require
from .simulate_full_family_tail_pilot import _slot_rows
from .workflow import ROOT, canonical_digest, digest, text_digest

IMPLEMENTATIONS = (
    "analysis/v3/run_final_calibration_outer.py",
    "analysis/v3/run_final_calibration_shard.py",
    "analysis/v3/multicoordinate_calibration.py",
    "analysis/v3/module_ecology.py",
    "analysis/v3/joint_partial_pooling.py",
    "analysis/v3/dependence_resampling.py",
    "analysis/v3/nuisance_design.py",
)


def verify_outer(path, contract, *, root=ROOT):
    """Every failed/non-estimable slot fails closed, never success-conditioned."""
    path = Path(path)
    report = json.loads(path.read_text(encoding="utf-8"))
    execution = report["execution_contract"]
    require(execution == json.loads((path.parent / "execution_contract.json").read_text(encoding="utf-8")),
            "Standalone and embedded execution contracts differ")
    require(execution["specification_canonical_sha256"] == canonical_digest(contract), "Execution contract differs")
    scenario, outer = execution["scenario"], execution["outer_replicate"]
    require(scenario in contract["scenarios"] and type(outer) is int
            and 0 <= outer < contract["sequential_outer_rule"]["maximum_outer_replicates_per_scenario"],
            "Invalid independent outer identity")
    seed = contract["seed"] + 100000 * contract["scenarios"].index(scenario) + outer
    require(execution["outer_seed"] == seed, "Outer seed changed")
    require(execution["implementation_sha256_text_lf"] == {p: text_digest(root / p) for p in IMPLEMENTATIONS},
            "Restored numerical implementation differs; do not combine estimator versions")
    for obj in (report, execution):
        require(obj.get("empirical_trait_environment_values_read") == 0 and obj.get("ecological_models_executed") == 0
                and obj.get("ecological_fitting_authorized") is False, "Synthetic boundary differs")
    blocks = json.loads((root / "analysis/v3/crossed_bootstrap_simulation_contract.json").read_text(encoding="utf-8"))["process_indices"]
    dimensions, replicates = contract["module_dimensions"], contract["bootstrap_replicates"]
    result = {}
    for case in report["cases"]:
        grid = case["grid_degrees"]
        require(grid in contract["grid_degrees"] and grid not in result, "Duplicate or unexpected grid")
        require((case["scenario"], case["outer_replicate"], case["outer_seed"]) == (scenario, outer, seed),
                "Case identity differs from execution")
        require(case["status"] == "FINAL_MULTICOORDINATE_OUTER_COMPLETE_NOT_AGGREGATED"
                and case["planned_bootstrap_replicates"] == replicates
                and case["estimable_shared_bootstrap_replicates"] == replicates
                and case.get("ecological_fitting_authorized") is False,
                "Incomplete outer: preserve failures without replacement or inference")
        records_path = path.parent / f"shared_draw_records_{grid}deg.jsonl"
        arrays_path = path.parent / f"shared_draws_{grid}deg.npz"
        require(digest(records_path) == case["shared_draw_records_sha256"]
                and digest(arrays_path) == case["shared_draws_sha256"], "Returned draw bytes differ")
        records = [json.loads(line) for line in records_path.read_text(encoding="utf-8").splitlines()]
        require(len(records) == replicates and all(
            r["replicate"] == i and r["seed"] == seed and r["grid_degrees"] == grid and r["status"] == "estimated"
            for i, r in enumerate(records)), "Missing, repeated, reordered or failed shared draw")
        summaries, truths, coverage = {}, {}, {}
        with np.load(arrays_path, allow_pickle=False) as saved:
            require(set(saved.files) == {f"{m}_{k}" for m in dimensions for k in ("point", "draws", "truth")},
                    "Returned coordinate inventory differs")
            for module, dim in dimensions.items():
                point, draws, truth = (saved[f"{module}_{k}"] for k in ("point", "draws", "truth"))
                shape = (3, dim, contract["predictors"])
                require(point.shape == truth.shape == shape and draws.shape == (replicates, *shape)
                        and all(np.isfinite(v).all() for v in (point, draws, truth)), "Invalid saved coefficient geometry")
                summaries[module] = summary = summarize_draws(point, draws, blocks)
                truths[module] = truth
                covered = (summary["basic_interval_low"] <= truth) & (truth <= summary["basic_interval_high"])
                coverage[module] = {"covered": int(covered.sum()), "total": int(covered.size), "rate": float(covered.mean())}
        require(case["interval_coverage"] == coverage, "Reported interval coverage differs from saved arrays")
        slots = _slot_rows(summaries, truths, blocks)
        reported = case["family_slots"]
        require(len(slots) == len(reported) == contract["family_slots"], "Incomplete family inventory")
        for computed, old in zip(slots, reported, strict=True):
            require(computed["candidate_status"] == "candidate_tail_estimated_not_calibrated"
                    and computed["candidate_probability"] is not None, "Non-estimable family slot; no calibration admission")
            for key in computed:
                if key in ("candidate_probability", "holm_probability"):
                    require(type(old.get(key)) in (float, int) and np.isfinite(old[key])
                            and np.isclose(computed[key], old[key], rtol=0, atol=1e-12), "Family probability differs")
                else:
                    require(computed[key] == old.get(key), "Family identity, null status, rank or rejection differs")
        false_count = sum(s["generating_null"] and s["holm_reject_point05"] for s in slots)
        true_count = sum(not s["generating_null"] and s["holm_reject_point05"] for s in slots)
        require(case["any_false_holm_rejection"] is (false_count > 0)
                and case["false_holm_rejections"] == false_count and case["true_holm_rejections"] == true_count,
                "Reported family error counts differ")
        total = sum(r["total"] for r in coverage.values())
        count = sum(r["covered"] for r in coverage.values())
        require(total == 351, "Independent-outer coverage denominator changed")
        result[grid] = {"false_family": int(false_count > 0), "covered": count, "total": total,
                        "coverage_fraction": count / total, "false_holm_rejections": false_count,
                        "true_holm_rejections": true_count}
    require(set(result) == set(contract["grid_degrees"]), "Missing frozen grid")
    return (scenario, outer), result


def verify_inventory(paths, expected_manifest, contract, *, root=ROOT):
    expected = {(r["scenario"], r["outer_replicate"]): r["sha256"] for r in expected_manifest}
    require(len(expected) == len(expected_manifest), "Duplicate expected outer identity")
    observed, manifests = {}, []
    for raw in paths:
        path = Path(raw)
        report = json.loads(path.read_text(encoding="utf-8"))
        execution = report["execution_contract"]
        key = (execution["scenario"], execution["outer_replicate"])
        require(key in expected and key not in observed, "Unexpected or duplicate outer report")
        require(digest(path) == expected[key], "Returned report SHA differs from batch manifest")
        _, observed[key] = verify_outer(path, contract, root=root)
        manifests.append({"scenario": key[0], "outer_replicate": key[1], "path": path.as_posix(), "sha256": expected[key]})
    require(set(observed) == set(expected), "Restore every planned outer; no success-conditioned subset")
    return observed, sorted(manifests, key=lambda r: (r["scenario"], r["outer_replicate"]))

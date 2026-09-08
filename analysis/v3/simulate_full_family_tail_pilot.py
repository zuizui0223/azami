"""Tail-resolved 36-slot family mechanics pilot on generated data only.

This module never reads empirical trait/environment values and cannot authorize ecology.
All three synthetic modules share each crossed bootstrap draw before Holm assembly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from .dependence_resampling import cohort_partition, crossed_draw, source_partition
from .environment_model import holm_complete_family, test_family
from .module_ecology import SCALES, fit_matched_module, summarize_draws
from .protected_artifacts import new_json, require
from .simulate_crossed_bootstrap import generate
from .workflow import ROOT, canonical_digest, digest, text_digest

SPEC = ROOT / "analysis/v3/full_family_tail_pilot_contract.json"
PRELIM_SPEC = ROOT / "analysis/v3/crossed_bootstrap_simulation_contract.json"


def _specs():
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    prelim = json.loads(PRELIM_SPEC.read_text(encoding="utf-8"))
    require(spec["modules"] == ["orientation", "visible_colour", "gross_shape"], "Unexpected module order")
    require(spec["processes"] == list(prelim["process_indices"]), "Process order differs from estimator contract")
    require(spec["scales"] == list(SCALES), "Scale order differs from estimator contract")
    require(spec["family_slots"] == len(test_family()) == 36, "Full family is not 36 slots")
    require(abs(spec["minimum_plus_one_probability"] - 1 / (spec["bootstrap_replicates"] + 1)) < 1e-15,
            "Tail resolution does not match planned draws")
    require(spec["minimum_plus_one_probability"] <= spec["first_holm_threshold"],
            "Pilot cannot resolve the first Holm threshold")
    return spec, prelim


def generated_modules(scenario: str, seed: int):
    """Return three response modules and exact generating coefficient arrays."""
    spec, _ = _specs()
    y, x, g, nuisance, components, lat, lon, truth = generate(scenario, seed)
    cfg = spec["module_response_coordinates"]["gross_shape"]
    weights = np.asarray(cfg["linear_combination_of_preliminary_coordinates"], dtype=float)
    require(weights.shape == (2,), "Gross-shape synthetic weights must have length two")
    rng = np.random.default_rng(seed + int(cfg["noise_seed_offset"]))
    gross = y @ weights + rng.normal(0.0, float(cfg["independent_noise_sd"]), len(y))
    responses = {
        "orientation": y[:, [0]],
        "visible_colour": y[:, [1]],
        "gross_shape": gross[:, None],
    }
    truths = {
        "orientation": truth[:, [0], :],
        "visible_colour": truth[:, [1], :],
        "gross_shape": (truth[:, 0, :] * weights[0] + truth[:, 1, :] * weights[1])[:, None, :],
    }
    return responses, truths, x, g, nuisance, components, lat, lon


def _slot_rows(module_summaries, truths, blocks):
    probabilities = {}
    raw_rows = {}
    for module, summary in module_summaries.items():
        for test in summary["candidate_tests"]:
            key = (module, test["process"], test["scale"])
            probabilities[key] = test["candidate_probability"]
            indices = blocks[test["process"]]
            target = truths[module][SCALES.index(test["scale"])][:, indices]
            raw_rows[key] = {
                "module": module,
                "process": test["process"],
                "question": test["scale"],
                "generating_null": bool(np.all(target == 0)),
                "candidate_probability": test["candidate_probability"],
                "candidate_status": test["status"],
                "candidate_rank": test.get("rank"),
            }
    require(set(probabilities) == set(test_family()), "Did not assemble all 36 planned slots")
    adjusted = holm_complete_family(probabilities)
    rows = []
    for key in test_family():
        row = raw_rows[key]
        row["holm_probability"] = adjusted[key]
        row["raw_reject_point05"] = (row["candidate_probability"] is not None
                                      and row["candidate_probability"] < 0.05)
        row["holm_reject_point05"] = (row["holm_probability"] is not None
                                       and row["holm_probability"] < 0.05)
        rows.append(row)
    return rows


def run(out: Path, scenario: str, outer_replicate: int):
    spec, prelim = _specs()
    require(scenario in spec["scenarios"], "Scenario outside pilot contract")
    require(0 <= outer_replicate < spec["replicates_per_scenario"], "Outer replicate outside pilot contract")
    require(not out.exists(), "Preserve earlier pilot run")
    out.mkdir(parents=True)

    scenario_index = spec["scenarios"].index(scenario)
    seed = int(spec["seed"]) + 10000 * scenario_index + int(outer_replicate)
    responses, truths, x, g, nuisance, components, lat, lon = generated_modules(scenario, seed)
    blocks = prelim["process_indices"]

    code_paths = [
        "analysis/v3/dependence_resampling.py",
        "analysis/v3/environment_model.py",
        "analysis/v3/joint_partial_pooling.py",
        "analysis/v3/module_ecology.py",
        "analysis/v3/simulate_crossed_bootstrap.py",
        "analysis/v3/simulate_full_family_tail_pilot.py",
    ]
    execution = {
        "specification": spec,
        "specification_canonical_sha256": canonical_digest(spec),
        "preliminary_generator_contract_canonical_sha256": canonical_digest(prelim),
        "implementation_sha256_text_lf": {path: text_digest(ROOT / path) for path in code_paths},
        "scenario": scenario,
        "outer_replicate": outer_replicate,
        "seed": seed,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    new_json(out / "execution_contract.json", execution)

    case_summaries = []
    for degrees in spec["grid_degrees"]:
        started = time.perf_counter()
        source = source_partition(np.arange(len(x)), g, components, lat, lon, grid_degrees=degrees)
        positions, _, _ = cohort_partition(source, np.arange(len(x)))
        points = {module: fit_matched_module(response, x, g, nuisance)
                  for module, response in responses.items()}
        draws = {
            module: np.full((spec["bootstrap_replicates"], *points[module]["coefficients"].shape), np.nan)
            for module in spec["modules"]
        }
        replicate_records = []
        for replicate in range(spec["bootstrap_replicates"]):
            row = {"replicate": replicate, "seed": seed, "grid_degrees": degrees}
            try:
                indices, copy_taxa = crossed_draw(source, positions, seed=seed, replicate=replicate)
                fitted = {}
                for module in spec["modules"]:
                    fitted[module] = fit_matched_module(
                        responses[module][indices], x[indices], copy_taxa, nuisance[indices]
                    )
                    require(fitted[module]["coefficients"].shape == points[module]["coefficients"].shape,
                            "Shared draw changed coefficient geometry")
                    require(np.isfinite(fitted[module]["coefficients"]).all(),
                            "Shared draw produced non-finite coefficients")
                for module in spec["modules"]:
                    draws[module][replicate] = fitted[module]["coefficients"]
                row.update(status="estimated", sampled_observations=len(indices),
                           sampled_taxon_copies=len(fitted[spec["modules"][0]]["taxa"]))
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                row.update(status="not_estimable", error_type=type(error).__name__, reason=str(error),
                           optimizer_attempts=getattr(error, "optimizer_attempts", []))
            replicate_records.append(row)
            if replicate % 100 == 0 or replicate + 1 == spec["bootstrap_replicates"]:
                print(json.dumps({"scenario": scenario, "outer_replicate": outer_replicate,
                                  "grid_degrees": degrees, "bootstrap_replicate": replicate,
                                  "status": row["status"]}), flush=True)

        estimable = sum(row["status"] == "estimated" for row in replicate_records)
        complete = estimable == spec["bootstrap_replicates"]
        draw_path = out / f"shared_draws_{degrees}deg.npz"
        np.savez_compressed(draw_path, **{f"{module}_draws": draws[module] for module in spec["modules"]},
                            **{f"{module}_point": points[module]["coefficients"] for module in spec["modules"]})
        records_path = out / f"shared_draw_records_{degrees}deg.jsonl"
        with records_path.open("x", encoding="utf-8", newline="\n") as handle:
            for row in replicate_records:
                handle.write(json.dumps(row, allow_nan=False) + "\n")

        summary = {
            "scenario": scenario,
            "outer_replicate": outer_replicate,
            "seed": seed,
            "grid_degrees": degrees,
            "observations": len(x),
            "planned_bootstrap_replicates": spec["bootstrap_replicates"],
            "estimable_shared_bootstrap_replicates": estimable,
            "shared_draws_sha256": digest(draw_path),
            "shared_draw_records_sha256": digest(records_path),
            "elapsed_seconds": time.perf_counter() - started,
            "status": "TAIL_RESOLVED_FULL36_FAMILY_ASSEMBLED_NOT_CALIBRATED" if complete
                      else "SHARED_BOOTSTRAP_INCOMPLETE_NO_FAMILY_INFERENCE",
            "ecological_fitting_authorized": False,
        }
        if complete:
            module_summaries = {
                module: summarize_draws(points[module]["coefficients"], draws[module], blocks)
                for module in spec["modules"]
            }
            slots = _slot_rows(module_summaries, truths, blocks)
            coverage = {}
            for module in spec["modules"]:
                s = module_summaries[module]
                covered = ((s["basic_interval_low"] <= truths[module])
                           & (truths[module] <= s["basic_interval_high"]))
                coverage[module] = {
                    "covered": int(covered.sum()),
                    "total": int(covered.size),
                    "rate": float(covered.mean()),
                }
            null_slots = [row for row in slots if row["generating_null"]]
            nonnull_slots = [row for row in slots if not row["generating_null"]]
            summary.update(
                family_slots=slots,
                interval_coverage=coverage,
                null_slots=len(null_slots),
                nonnull_slots=len(nonnull_slots),
                false_holm_rejections=sum(row["holm_reject_point05"] for row in null_slots),
                any_false_holm_rejection=any(row["holm_reject_point05"] for row in null_slots),
                true_holm_rejections=sum(row["holm_reject_point05"] for row in nonnull_slots),
                tail_resolution_hits=sum(
                    row["candidate_probability"] is not None
                    and row["candidate_probability"] <= spec["first_holm_threshold"] for row in slots
                ),
            )
        new_json(out / f"summary_{degrees}deg.json", summary)
        case_summaries.append(summary)

    public = {
        "status": "FULL36_TAIL_MECHANICS_PILOT_EXECUTED_NO_ECOLOGY",
        "execution_contract": execution,
        "cases": case_summaries,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
        "qualification_boundary": spec["qualification_boundary"],
    }
    new_json(out / "public_report.json", public)
    print(json.dumps({"status": public["status"], "scenario": scenario,
                      "outer_replicate": outer_replicate, "cases": len(case_summaries)}), flush=True)
    return public


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--outer-replicate", type=int, required=True)
    args = parser.parse_args()
    run(args.out, args.scenario, args.outer_replicate)

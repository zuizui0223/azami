"""Execute one frozen final multicoordinate outer calibration dataset.

Synthetic-only. All 1/8/4 module coordinates share every crossed draw before the
complete 36-slot Holm family is assembled on both fixed spatial grids.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from .dependence_resampling import cohort_partition, crossed_draw, source_partition
from .environment_model import test_family
from .module_ecology import fit_matched_module, summarize_draws
from .multicoordinate_calibration import MODULES, generate
from .protected_artifacts import new_json, require
from .run_final_calibration_shard import definition as shard_definition, outer_seed
from .simulate_full_family_tail_pilot import _slot_rows
from .workflow import ROOT, canonical_digest, digest, text_digest

SPEC = ROOT / "analysis/v3/final_module_calibration_contract.json"


def run(out: Path, scenario: str, outer: int):
    spec, prelim = shard_definition()
    require(not out.exists(), "Preserve earlier final calibration outer")
    out.mkdir(parents=True)
    seed = outer_seed(spec, scenario, outer)
    data = generate(scenario, seed)
    x, g, nuisance = data["predictors"], data["taxa"], data["nuisance"]
    blocks = prelim["process_indices"]

    execution = {
        "specification_canonical_sha256": canonical_digest(spec),
        "scenario": scenario,
        "outer_replicate": outer,
        "outer_seed": seed,
        "implementation_sha256_text_lf": {path: text_digest(ROOT / path) for path in [
            "analysis/v3/run_final_calibration_outer.py",
            "analysis/v3/run_final_calibration_shard.py",
            "analysis/v3/multicoordinate_calibration.py",
            "analysis/v3/module_ecology.py",
            "analysis/v3/joint_partial_pooling.py",
            "analysis/v3/dependence_resampling.py",
            "analysis/v3/nuisance_design.py",
        ]},
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    new_json(out / "execution_contract.json", execution)

    cases = []
    for grid in spec["grid_degrees"]:
        started = time.perf_counter()
        source = source_partition(np.arange(len(x)), g, data["components"], data["latitude"],
                                  data["longitude"], grid_degrees=grid)
        positions, _, _ = cohort_partition(source, np.arange(len(x)))
        points = {module: fit_matched_module(data["responses"][module], x, g, nuisance)
                  for module in MODULES}
        for module in MODULES:
            require(points[module]["coefficients"].shape == data["truths"][module].shape,
                    "Point coefficient geometry differs from generating truth")
        draws = {module: np.full((spec["bootstrap_replicates"], *points[module]["coefficients"].shape), np.nan)
                 for module in MODULES}
        records = []
        for replicate in range(spec["bootstrap_replicates"]):
            row = {"replicate": replicate, "seed": seed, "grid_degrees": grid}
            try:
                indices, copy_taxa = crossed_draw(source, positions, seed=seed, replicate=replicate)
                fitted = {}
                for module in MODULES:
                    fitted[module] = fit_matched_module(
                        data["responses"][module][indices], x[indices], copy_taxa, nuisance[indices]
                    )
                    require(fitted[module]["coefficients"].shape == points[module]["coefficients"].shape,
                            "Shared draw changed coefficient geometry")
                    require(np.isfinite(fitted[module]["coefficients"]).all(),
                            "Shared draw produced non-finite coefficients")
                for module in MODULES:
                    draws[module][replicate] = fitted[module]["coefficients"]
                row.update(status="estimated", sampled_observations=len(indices),
                           sampled_taxon_copies=len(fitted[MODULES[0]]["taxa"]))
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                row.update(status="not_estimable", error_type=type(error).__name__, reason=str(error),
                           optimizer_attempts=getattr(error, "optimizer_attempts", []))
            records.append(row)
            if replicate % 100 == 0 or replicate + 1 == spec["bootstrap_replicates"]:
                print(json.dumps({"scenario": scenario, "outer": outer, "grid": grid,
                                  "replicate": replicate, "status": row["status"]}), flush=True)

        estimable = sum(row["status"] == "estimated" for row in records)
        record_path = out / f"shared_draw_records_{grid}deg.jsonl"
        with record_path.open("x", encoding="utf-8", newline="\n") as handle:
            for row in records:
                handle.write(json.dumps(row, allow_nan=False) + "\n")
        draw_path = out / f"shared_draws_{grid}deg.npz"
        np.savez_compressed(draw_path,
                            **{f"{m}_draws": draws[m] for m in MODULES},
                            **{f"{m}_point": points[m]["coefficients"] for m in MODULES},
                            **{f"{m}_truth": data["truths"][m] for m in MODULES})

        complete = estimable == spec["bootstrap_replicates"]
        case = {
            "scenario": scenario,
            "outer_replicate": outer,
            "outer_seed": seed,
            "grid_degrees": grid,
            "observations": len(x),
            "planned_bootstrap_replicates": spec["bootstrap_replicates"],
            "estimable_shared_bootstrap_replicates": estimable,
            "shared_draw_records_sha256": digest(record_path),
            "shared_draws_sha256": digest(draw_path),
            "elapsed_seconds": time.perf_counter() - started,
            "status": "FINAL_MULTICOORDINATE_OUTER_COMPLETE_NOT_AGGREGATED" if complete
                      else "FINAL_MULTICOORDINATE_OUTER_INCOMPLETE_NO_INFERENCE",
            "ecological_fitting_authorized": False,
        }
        if complete:
            summaries = {module: summarize_draws(points[module]["coefficients"], draws[module], blocks)
                         for module in MODULES}
            slots = _slot_rows(summaries, data["truths"], blocks)
            require(len(slots) == len(test_family()) == spec["family_slots"],
                    "Final outer did not assemble the complete family")
            coverage = {}
            for module in MODULES:
                s = summaries[module]
                covered = ((s["basic_interval_low"] <= data["truths"][module])
                           & (data["truths"][module] <= s["basic_interval_high"]))
                coverage[module] = {"covered": int(covered.sum()), "total": int(covered.size),
                                    "rate": float(covered.mean())}
            null_slots = [row for row in slots if row["generating_null"]]
            nonnull_slots = [row for row in slots if not row["generating_null"]]
            case.update(
                family_slots=slots,
                interval_coverage=coverage,
                false_holm_rejections=sum(row["holm_reject_point05"] for row in null_slots),
                any_false_holm_rejection=any(row["holm_reject_point05"] for row in null_slots),
                true_holm_rejections=sum(row["holm_reject_point05"] for row in nonnull_slots),
            )
        new_json(out / f"summary_{grid}deg.json", case)
        cases.append(case)

    report = {
        "status": "FINAL_MULTICOORDINATE_CALIBRATION_OUTER_EXECUTED_NO_ECOLOGY",
        "execution_contract": execution,
        "cases": cases,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
        "qualification_boundary": "One synthetic outer dataset only. Sequential admission requires the frozen 25-per-scenario batches, minimum 100 per scenario and precision/admission rules in the final contract plus later realized-design calibration.",
    }
    new_json(out / "public_report.json", report)
    print(json.dumps({"status": report["status"], "scenario": scenario, "outer": outer,
                      "complete_cases": sum(c["status"].startswith("FINAL_MULTICOORDINATE_OUTER_COMPLETE") for c in cases)}),
          flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--outer", type=int, required=True)
    args = parser.parse_args()
    run(args.out, args.scenario, args.outer)

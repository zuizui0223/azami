"""Execute one frozen final-calibration module x grid x outer shard.

Synthetic generated data only. Sharding is computational: deterministic crossed
resampling is fingerprinted so later aggregation can prove that all three modules
used the same 999 source-factor draws for a scenario/outer/grid family.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .dependence_resampling import cohort_partition, crossed_draw, source_partition
from .environment_model import test_family
from .module_ecology import fit_matched_module, summarize_draws
from .multicoordinate_calibration import MODULES, definition as calibration_definition, generate
from .protected_artifacts import new_json, require
from .simulate_crossed_bootstrap import SPEC as PRELIM_SPEC
from .workflow import ROOT, canonical_digest, digest, text_digest

SPEC = ROOT / "analysis/v3/final_module_calibration_contract.json"


def definition():
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    prelim = json.loads(PRELIM_SPEC.read_text(encoding="utf-8"))
    calibration, env = calibration_definition()
    require(spec["seed"] == 2026090893, "Final calibration seed changed")
    require(spec["bootstrap_replicates"] == 999, "Final calibration draw count changed")
    require(spec["module_dimensions"] == {"orientation": 1, "visible_colour": 8, "gross_shape": 4},
            "Final module geometry changed")
    require(spec["predictors"] == 9 == len(env["variables"]) and spec["nuisance_columns"] == 13,
            "Final predictor/nuisance geometry changed")
    require(spec["family_slots"] == len(test_family()) == 36, "Primary family changed")
    require(spec["scenarios"] == calibration["scenarios"], "Calibration scenario order changed")
    require(spec["grid_degrees"] == calibration["grid_degrees"], "Calibration grids changed")
    require(list(prelim["process_indices"]) == ["wetting_moisture", "radiation", "heat_drying", "mechanical"],
            "Process blocks changed")
    require(not spec["ecological_fitting_authorized"], "Calibration cannot authorize ecology directly")
    return spec, prelim


def outer_seed(spec, scenario, outer):
    require(scenario in spec["scenarios"], "Scenario outside final calibration")
    require(isinstance(outer, int) and not isinstance(outer, bool) and 0 <= outer < spec["sequential_outer_rule"]["maximum_outer_replicates_per_scenario"],
            "Outer replicate outside frozen sequential range")
    return int(spec["seed"]) + 100000 * spec["scenarios"].index(scenario) + outer


def run(out: Path, scenario: str, outer: int, grid: int, module: str):
    spec, prelim = definition()
    require(grid in spec["grid_degrees"], "Grid outside final calibration")
    require(module in MODULES, "Module outside final calibration")
    require(not out.exists(), "Preserve earlier calibration shard")
    out.mkdir(parents=True)
    seed = outer_seed(spec, scenario, outer)
    data = generate(scenario, seed)
    response = data["responses"][module]
    truth = data["truths"][module]
    require(response.shape[1] == spec["module_dimensions"][module], "Generated module dimension changed")

    source = source_partition(np.arange(len(data["taxa"])), data["taxa"], data["components"],
                              data["latitude"], data["longitude"], grid_degrees=grid)
    positions, _, _ = cohort_partition(source, np.arange(len(data["taxa"])))
    point = fit_matched_module(response, data["predictors"], data["taxa"], data["nuisance"])
    require(point["coefficients"].shape == truth.shape, "Point coefficient geometry differs from truth")

    draws = np.full((spec["bootstrap_replicates"], *point["coefficients"].shape), np.nan)
    records = []
    fingerprint = hashlib.sha256()
    started = time.perf_counter()
    for replicate in range(spec["bootstrap_replicates"]):
        row = {"replicate": replicate, "seed": seed, "grid_degrees": grid}
        try:
            indices, copy_taxa = crossed_draw(source, positions, seed=seed, replicate=replicate)
            fingerprint.update(np.asarray(indices, dtype=np.int64).tobytes())
            fingerprint.update(b"\0")
            fingerprint.update("\n".join(copy_taxa.tolist()).encode("utf-8"))
            fingerprint.update(b"\xff")
            fitted = fit_matched_module(response[indices], data["predictors"][indices], copy_taxa,
                                         data["nuisance"][indices])
            require(fitted["coefficients"].shape == point["coefficients"].shape,
                    "Crossed draw changed coefficient geometry")
            require(np.isfinite(fitted["coefficients"]).all(), "Crossed draw produced non-finite coefficients")
            draws[replicate] = fitted["coefficients"]
            row.update(status="estimated", sampled_observations=len(indices), sampled_taxon_copies=len(fitted["taxa"]))
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            row.update(status="not_estimable", error_type=type(error).__name__, reason=str(error),
                       optimizer_attempts=getattr(error, "optimizer_attempts", []))
        records.append(row)
        if replicate % 100 == 0 or replicate + 1 == spec["bootstrap_replicates"]:
            print(json.dumps({"scenario": scenario, "outer": outer, "grid": grid, "module": module,
                              "replicate": replicate, "status": row["status"]}), flush=True)

    estimable = sum(row["status"] == "estimated" for row in records)
    record_path = out / "draw_records.jsonl"
    with record_path.open("x", encoding="utf-8", newline="\n") as handle:
        for row in records:
            handle.write(json.dumps(row, allow_nan=False) + "\n")
    np.savez_compressed(out / "draws.npz", draws=draws, point=point["coefficients"], truth=truth)

    complete = estimable == spec["bootstrap_replicates"]
    report = {
        "status": "FINAL_CALIBRATION_SHARD_COMPLETE_NOT_AGGREGATED" if complete else "FINAL_CALIBRATION_SHARD_INCOMPLETE_NO_INFERENCE",
        "scenario": scenario,
        "outer_replicate": outer,
        "outer_seed": seed,
        "grid_degrees": grid,
        "module": module,
        "module_dimension": spec["module_dimensions"][module],
        "observations": len(response),
        "planned_bootstrap_replicates": spec["bootstrap_replicates"],
        "estimable_bootstrap_replicates": estimable,
        "resampling_sequence_sha256": fingerprint.hexdigest(),
        "draw_records_sha256": digest(record_path),
        "draws_npz_sha256": digest(out / "draws.npz"),
        "elapsed_seconds": time.perf_counter() - started,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    if complete:
        summary = summarize_draws(point["coefficients"], draws, prelim["process_indices"])
        covered = ((summary["basic_interval_low"] <= truth) & (truth <= summary["basic_interval_high"]))
        report.update(
            candidate_tests=summary["candidate_tests"],
            point_coefficients=point["coefficients"].tolist(),
            generating_coefficients=truth.tolist(),
            basic_interval_low=summary["basic_interval_low"].tolist(),
            basic_interval_high=summary["basic_interval_high"].tolist(),
            coefficient_coverage={"covered": int(covered.sum()), "total": int(covered.size), "rate": float(covered.mean())},
        )
    execution = {
        "specification_canonical_sha256": canonical_digest(spec),
        "implementation_sha256_text_lf": {path: text_digest(ROOT / path) for path in [
            "analysis/v3/run_final_calibration_shard.py", "analysis/v3/multicoordinate_calibration.py",
            "analysis/v3/module_ecology.py", "analysis/v3/joint_partial_pooling.py",
            "analysis/v3/dependence_resampling.py", "analysis/v3/nuisance_design.py"]},
        "scenario": scenario, "outer_replicate": outer, "grid_degrees": grid, "module": module,
        "outer_seed": seed, "ecological_fitting_authorized": False,
    }
    new_json(out / "execution_contract.json", execution)
    new_json(out / "public_report.json", report)
    print(json.dumps({"status": report["status"], "scenario": scenario, "outer": outer,
                      "grid": grid, "module": module, "estimable": estimable}), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--outer", type=int, required=True)
    parser.add_argument("--grid", type=int, required=True)
    parser.add_argument("--module", required=True)
    args = parser.parse_args()
    run(args.out, args.scenario, args.outer, args.grid, args.module)

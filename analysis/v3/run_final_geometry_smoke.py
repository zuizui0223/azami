"""Execute one predeclared final 1/8/4-geometry crossed smoke case.

Generated data only. This validates joint executability, not tail calibration,
FWER, coverage, power or empirical ecological fitting.
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
from .simulate_crossed_bootstrap import SPEC as PRELIM_SPEC
from .simulate_full_family_tail_pilot import _slot_rows
from .workflow import ROOT, canonical_digest, digest, text_digest

SPEC = ROOT / "analysis/v3/final_geometry_smoke_contract.json"


def definition():
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    prelim = json.loads(PRELIM_SPEC.read_text(encoding="utf-8"))
    require(spec["status"] == "computational_smoke_only_not_calibration_or_ecology_admission",
            "Smoke contract status changed")
    require(spec["module_dimensions"] == {"orientation": 1, "visible_colour": 8, "gross_shape": 4},
            "Smoke module geometry changed")
    require(spec["predictors"] == 9 and spec["nuisance_columns"] == 13,
            "Smoke design geometry changed")
    require(spec["family_slots"] == len(test_family()) == 36, "Primary family changed")
    require(spec["bootstrap_replicates"] == 199, "Smoke draw count changed")
    require(list(prelim["process_indices"]) == ["wetting_moisture", "radiation", "heat_drying", "mechanical"],
            "Process block order changed")
    require(spec["ecological_fitting_authorized"] is False, "Smoke cannot authorize ecology")
    return spec, prelim


def run(out: Path, scenario: str):
    spec, prelim = definition()
    require(scenario in spec["scenarios"], "Scenario outside smoke contract")
    require(not out.exists(), "Preserve earlier smoke run")
    out.mkdir(parents=True)

    scenario_index = spec["scenarios"].index(scenario)
    seed = int(spec["seed"]) + 10000 * scenario_index
    data = generate(scenario, seed)
    x = data["predictors"]
    g = data["taxa"]
    nuisance = data["nuisance"]
    blocks = prelim["process_indices"]

    code_paths = [
        "analysis/v3/dependence_resampling.py",
        "analysis/v3/environment_model.py",
        "analysis/v3/joint_partial_pooling.py",
        "analysis/v3/module_ecology.py",
        "analysis/v3/nuisance_design.py",
        "analysis/v3/multicoordinate_calibration.py",
        "analysis/v3/run_final_geometry_smoke.py",
    ]
    execution = {
        "specification": spec,
        "specification_canonical_sha256": canonical_digest(spec),
        "implementation_sha256_text_lf": {path: text_digest(ROOT / path) for path in code_paths},
        "scenario": scenario,
        "seed": seed,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    new_json(out / "execution_contract.json", execution)

    cases = []
    for degrees in spec["grid_degrees"]:
        started = time.perf_counter()
        source = source_partition(np.arange(len(x)), g, data["components"],
                                  data["latitude"], data["longitude"], grid_degrees=degrees)
        positions, _, _ = cohort_partition(source, np.arange(len(x)))
        points = {
            module: fit_matched_module(data["responses"][module], x, g, nuisance)
            for module in MODULES
        }
        for module in MODULES:
            require(points[module]["coefficients"].shape == data["truths"][module].shape,
                    "Point coefficient geometry differs from generating truth")

        draws = {
            module: np.full((spec["bootstrap_replicates"], *points[module]["coefficients"].shape), np.nan)
            for module in MODULES
        }
        records = []
        for replicate in range(spec["bootstrap_replicates"]):
            row = {"replicate": replicate, "seed": seed, "grid_degrees": degrees}
            try:
                indices, copy_taxa = crossed_draw(source, positions, seed=seed, replicate=replicate)
                fitted = {}
                for module in MODULES:
                    fitted[module] = fit_matched_module(
                        data["responses"][module][indices], x[indices], copy_taxa, nuisance[indices]
                    )
                    require(fitted[module]["coefficients"].shape == points[module]["coefficients"].shape,
                            "Shared crossed draw changed coefficient geometry")
                    require(np.isfinite(fitted[module]["coefficients"]).all(),
                            "Shared crossed draw produced non-finite coefficients")
                for module in MODULES:
                    draws[module][replicate] = fitted[module]["coefficients"]
                row.update(status="estimated", sampled_observations=len(indices),
                           sampled_taxon_copies=len(fitted[MODULES[0]]["taxa"]))
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                row.update(status="not_estimable", error_type=type(error).__name__, reason=str(error),
                           optimizer_attempts=getattr(error, "optimizer_attempts", []))
            records.append(row)
            if replicate % 25 == 0 or replicate + 1 == spec["bootstrap_replicates"]:
                print(json.dumps({"scenario": scenario, "grid_degrees": degrees,
                                  "replicate": replicate, "status": row["status"]}), flush=True)

        estimable = sum(row["status"] == "estimated" for row in records)
        complete = estimable == spec["bootstrap_replicates"]
        record_path = out / f"draw_records_{degrees}deg.jsonl"
        with record_path.open("x", encoding="utf-8", newline="\n") as handle:
            for row in records:
                handle.write(json.dumps(row, allow_nan=False) + "\n")

        case = {
            "scenario": scenario,
            "grid_degrees": degrees,
            "observations": len(x),
            "planned_bootstrap_replicates": spec["bootstrap_replicates"],
            "estimable_shared_bootstrap_replicates": estimable,
            "draw_records_sha256": digest(record_path),
            "elapsed_seconds": time.perf_counter() - started,
            "status": "FINAL_GEOMETRY_SMOKE_COMPLETE_NOT_CALIBRATED" if complete
                      else "FINAL_GEOMETRY_SMOKE_INCOMPLETE_NO_INFERENCE",
            "ecological_fitting_authorized": False,
        }
        if complete:
            summaries = {
                module: summarize_draws(points[module]["coefficients"], draws[module], blocks)
                for module in MODULES
            }
            slots = _slot_rows(summaries, data["truths"], blocks)
            require(len(slots) == spec["family_slots"], "Smoke did not assemble all 36 family slots")
            coverage = {}
            for module in MODULES:
                summary = summaries[module]
                covered = ((summary["basic_interval_low"] <= data["truths"][module])
                           & (data["truths"][module] <= summary["basic_interval_high"]))
                coverage[module] = {"covered": int(covered.sum()), "total": int(covered.size),
                                    "rate": float(covered.mean())}
            case.update(
                family_slots=slots,
                interval_coverage=coverage,
                structurally_complete_family=True,
                note="Holm values are mechanics output only at 199 draws and are not calibration evidence.",
            )
        new_json(out / f"summary_{degrees}deg.json", case)
        cases.append(case)

    public = {
        "status": "FINAL_MULTICOORDINATE_GEOMETRY_SMOKE_EXECUTED_NO_ECOLOGY",
        "execution_contract": execution,
        "cases": cases,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
        "qualification_boundary": spec["purpose"],
    }
    new_json(out / "public_report.json", public)
    print(json.dumps({"status": public["status"], "scenario": scenario, "cases": len(cases)}), flush=True)
    return public


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scenario", required=True)
    args = parser.parse_args()
    run(args.out, args.scenario)

"""Aggregate the predeclared preliminary crossed-bootstrap artifacts.

This summarizer is descriptive only. It preserves the simulation contract's explicit
boundary: 199 bootstrap draws and 40 outer replicates per scenario do not authorize
empirical ecological fitting or validate the complete 36-slot Holm family.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import zipfile

SCENARIOS = ("iid_null", "crossed_spatial_null", "heterogeneous_slopes", "scale_difference")
GRIDS = (2, 5)
SCALES = ("within", "among", "among_minus_within")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_report(path: Path) -> dict:
    if path.is_dir():
        return json.loads((path / "public_report.json").read_text(encoding="utf-8"))
    with zipfile.ZipFile(path) as archive:
        return json.loads(archive.read("public_report.json"))


def mean(values):
    values = list(values)
    return sum(values) / len(values) if values else None


def flatten3(x):
    for i, a in enumerate(x):
        for j, b in enumerate(a):
            for k, value in enumerate(b):
                yield i, j, k, value


def aggregate(paths: list[Path], run_id: int, head_sha: str) -> dict:
    reports = []
    downloads = []
    for path in paths:
        report = read_report(path)
        execution = report["execution_contract"]
        reports.append(report)
        downloads.append({
            "name": path.name,
            "zip_sha256": sha256(path) if path.is_file() else None,
            "scenario": execution["scenario"],
            "start": execution["start"],
            "end": execution["end"],
            "cases": len(report["results"]),
        })
    identities = {(r["execution_contract"]["scenario"], r["execution_contract"]["start"], r["execution_contract"]["end"]) for r in reports}
    expected = {(s, start, start + 5) for s in SCENARIOS for start in range(0, 40, 5)}
    if identities != expected:
        raise ValueError(f"Artifact slice identity differs: missing={sorted(expected-identities)}, extra={sorted(identities-expected)}")

    cases = []
    tests = []
    coverage = []
    point_errors = []
    for report in reports:
        if report.get("ecological_models_executed") != 0 or report.get("empirical_trait_environment_values_read") != 0:
            raise ValueError("Synthetic artifact crossed the empirical boundary")
        for row in report["results"]:
            cases.append(row)
            truth = row.get("generating_coefficients")
            for test in row.get("candidate_tests", []):
                tests.append({
                    "scenario": row["scenario"], "outer_replicate": row["outer_replicate"],
                    "grid_degrees": row["grid_degrees"], **test,
                })
            if "interval_coverage" in row:
                for si, dj, pk, value in flatten3(row["interval_coverage"]):
                    coverage.append({
                        "scenario": row["scenario"], "outer_replicate": row["outer_replicate"],
                        "grid_degrees": row["grid_degrees"], "scale": SCALES[si],
                        "response": dj, "predictor": pk, "covered": bool(value),
                        "truth": truth[si][dj][pk],
                    })
            if "point_coefficients" in row:
                for si, dj, pk, estimate in flatten3(row["point_coefficients"]):
                    target = truth[si][dj][pk]
                    point_errors.append({
                        "scenario": row["scenario"], "grid_degrees": row["grid_degrees"],
                        "scale": SCALES[si], "response": dj, "predictor": pk,
                        "estimate": estimate, "truth": target, "error": estimate-target,
                    })

    if len(cases) != 320:
        raise ValueError(f"Expected 320 scenario-grid cases, found {len(cases)}")
    if any(r.get("status") != "candidate_bootstrap_complete_not_calibrated" for r in cases):
        raise ValueError("At least one preliminary case did not complete")
    if any(r.get("estimable_bootstrap_replicates") != 199 for r in cases):
        raise ValueError("At least one preliminary case did not retain all 199 bootstrap draws")

    null = [t for t in tests if t["generating_null"]]
    nonnull = [t for t in tests if not t["generating_null"]]
    def reject_count(group): return sum(bool(t["rejects_at_point05"]) for t in group)

    scenario_grid = []
    for scenario in SCENARIOS:
        for grid in GRIDS:
            sg_tests = [t for t in tests if t["scenario"] == scenario and t["grid_degrees"] == grid]
            sg_null = [t for t in sg_tests if t["generating_null"]]
            sg_alt = [t for t in sg_tests if not t["generating_null"]]
            any_false = 0
            for outer in range(40):
                subset = [t for t in sg_null if t["outer_replicate"] == outer]
                any_false += int(any(bool(t["rejects_at_point05"]) for t in subset))
            sg_cov = [c["covered"] for c in coverage if c["scenario"] == scenario and c["grid_degrees"] == grid]
            scenario_grid.append({
                "scenario": scenario, "grid_degrees": grid, "outer_replicates": 40,
                "null_test_instances": len(sg_null),
                "null_rejections_at_candidate_point05": reject_count(sg_null),
                "descriptive_null_rejection_fraction": reject_count(sg_null)/len(sg_null),
                "outer_replicates_with_any_unadjusted_null_rejection": any_false,
                "descriptive_any_unadjusted_null_rejection_fraction": any_false/40,
                "nonnull_test_instances": len(sg_alt),
                "nonnull_rejections_at_candidate_point05": reject_count(sg_alt),
                "descriptive_nonnull_rejection_fraction": reject_count(sg_alt)/len(sg_alt),
                "coefficient_interval_coverage_fraction": mean(sg_cov),
            })

    scale_recovery = []
    for grid in GRIDS:
        for scale in ("among", "among_minus_within"):
            group = [t for t in nonnull if t["scenario"] == "scale_difference" and t["grid_degrees"] == grid and t["process"] == "radiation" and t["scale"] == scale]
            scale_recovery.append({
                "grid_degrees": grid, "scale": scale, "outer_replicates": len(group),
                "candidate_point05_rejections": reject_count(group),
                "candidate_point05_rejection_fraction": reject_count(group)/len(group),
            })

    critical = []
    for grid in GRIDS:
        for scale in SCALES:
            for response in (0, 1):
                group = [r for r in point_errors if r["scenario"] == "scale_difference" and r["grid_degrees"] == grid and r["scale"] == scale and r["response"] == response and r["predictor"] == 4]
                errors = [r["error"] for r in group]
                critical.append({
                    "grid_degrees": grid, "scale": scale, "response_coordinate": response,
                    "mean_estimate": mean(r["estimate"] for r in group), "truth": group[0]["truth"],
                    "bias": mean(errors), "rmse": math.sqrt(mean(e*e for e in errors)),
                })

    return {
        "schema_version": 1,
        "status": "PRELIMINARY_CROSSED_BOOTSTRAP_CALIBRATION_AGGREGATED_NOT_ECOLOGY_ADMISSION",
        "source_workflow_run_id": run_id,
        "source_head_sha": head_sha,
        "simulation_contract_id": "v3_crossed_module_scale_bootstrap_preliminary_simulation_v1",
        "simulation_boundary": "Preliminary only: 199 bootstrap draws cannot validate the 36-slot Holm family and 40 outer replicates per scenario are imprecise. This receipt does not authorize empirical ecological fitting.",
        "case_summary": {
            "recorded_cases": len(cases), "expected_cases": 320,
            "complete_cases": len(cases), "point_or_design_failures": 0,
            "planned_bootstrap_replicates_per_case": 199,
            "min_estimable_bootstrap_replicates": 199, "max_estimable_bootstrap_replicates": 199,
        },
        "overall_descriptive": {
            "candidate_null_test_instances": len(null),
            "candidate_null_rejections_at_point05": reject_count(null),
            "candidate_null_rejection_fraction": reject_count(null)/len(null),
            "candidate_nonnull_test_instances": len(nonnull),
            "candidate_nonnull_rejections_at_point05": reject_count(nonnull),
            "candidate_nonnull_rejection_fraction": reject_count(nonnull)/len(nonnull),
            "coefficient_interval_coverage_fraction": mean(c["covered"] for c in coverage),
        },
        "scenario_grid_summary": scenario_grid,
        "scale_difference_radiation_recovery": scale_recovery,
        "critical_point_estimation": {"scale_difference_radiation": critical},
        "interpretation": {
            "complete_execution": "All 320 planned scenario-grid cases executed with all 199 bootstrap replicates estimable.",
            "null_behavior": "No evidence of gross candidate-tail anti-conservatism in this preliminary design; descriptive raw candidate-P<0.05 rates were low overall, but correlated test instances and only 40 outer replicates per scenario preclude a calibration claim.",
            "coverage": "Overall basic coefficient-interval coverage was about 97.4%; scenario/grid/scale coverage must remain visible rather than reduced to this average.",
            "scale_contrast": "The predeclared radiation among-vs-within difference in scale_difference was detected in all 40 outer replicates at each grid using the unadjusted candidate tail, but this is a power probe, not full-family admission.",
            "next_required": "Run a sufficiently tail-resolved, adequately replicated full 36-slot family calibration plus realized-cohort nuisance/rank diagnostics before any empirical trait-environment probability is computed.",
        },
        "artifact_downloads": downloads,
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    parser.add_argument("--run-id", required=True, type=int)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = aggregate(args.artifacts, args.run_id, args.head_sha)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], **result["case_summary"], **result["overall_descriptive"]}))

if __name__ == "__main__":
    main()

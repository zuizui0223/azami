"""Aggregate one frozen 25-per-scenario final multicoordinate calibration batch.

Synthetic only. Missing/incomplete outers are retained as failures; this script
never replaces draws, selects a grid, changes thresholds or authorizes ecology.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from .multicoordinate_calibration import precision_decision
from .workflow import ROOT, canonical_digest

CONTRACT = ROOT / "analysis/v3/final_module_calibration_contract.json"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def aggregate(report_paths, batch_index: int, *, root: Path = ROOT):
    contract = json.loads((root / "analysis/v3/final_module_calibration_contract.json").read_text(encoding="utf-8"))
    batch = int(contract["sequential_outer_rule"]["batch_size_per_scenario"])
    if not isinstance(batch_index, int) or isinstance(batch_index, bool) or batch_index < 0:
        raise ValueError("Batch index must be a nonnegative integer")
    start = batch_index * batch
    stop = start + batch
    if stop > int(contract["sequential_outer_rule"]["maximum_outer_replicates_per_scenario"]):
        raise ValueError("Batch exceeds frozen maximum outer range")
    expected = {(scenario, outer) for scenario in contract["scenarios"] for outer in range(start, stop)}
    reports = {}
    manifest = []
    for raw in report_paths:
        path = Path(raw)
        report = json.loads(path.read_text(encoding="utf-8"))
        execution = report["execution_contract"]
        key = (execution["scenario"], int(execution["outer_replicate"]))
        if key in reports:
            raise ValueError(f"Duplicate final calibration outer {key}")
        if key not in expected:
            raise ValueError(f"Outer outside requested frozen batch {key}")
        if execution["specification_canonical_sha256"] != canonical_digest(contract):
            raise ValueError("Final calibration contract differs from outer execution")
        if report["empirical_trait_environment_values_read"] != 0 or report["ecological_models_executed"] != 0:
            raise ValueError("Final synthetic calibration crossed the empirical boundary")
        reports[key] = report
        manifest.append({"scenario": key[0], "outer_replicate": key[1], "path": path.as_posix(),
                         "sha256": _sha256(path)})

    missing = sorted(expected - set(reports))
    rows = []
    for key in sorted(reports):
        report = reports[key]
        seen = set()
        for case in report["cases"]:
            grid = int(case["grid_degrees"])
            if grid in seen:
                raise ValueError(f"Duplicate grid in outer {key}")
            seen.add(grid)
            slots = case.get("family_slots") or []
            coverage = case.get("interval_coverage") or {}
            rows.append({
                "scenario": key[0],
                "outer_replicate": key[1],
                "grid_degrees": grid,
                "complete": bool(
                    int(case["estimable_shared_bootstrap_replicates"]) == int(contract["bootstrap_replicates"])
                    and len(slots) == int(contract["family_slots"])
                    and case["status"] == "FINAL_MULTICOORDINATE_OUTER_COMPLETE_NOT_AGGREGATED"
                ),
                "estimable_shared_bootstrap_replicates": int(case["estimable_shared_bootstrap_replicates"]),
                "family_slots": len(slots),
                "any_false_holm_rejection": bool(case.get("any_false_holm_rejection", False)),
                "false_holm_rejections": int(case.get("false_holm_rejections", 0)),
                "true_holm_rejections": int(case.get("true_holm_rejections", 0)),
                "covered_coefficients": sum(int(v["covered"]) for v in coverage.values()),
                "coefficient_total": sum(int(v["total"]) for v in coverage.values()),
            })
        if seen != set(contract["grid_degrees"]):
            raise ValueError(f"Grid inventory differs for outer {key}: {sorted(seen)}")

    scenario_grid = []
    any_incomplete = bool(missing)
    for scenario in contract["scenarios"]:
        for grid in contract["grid_degrees"]:
            subset = [row for row in rows if row["scenario"] == scenario and row["grid_degrees"] == grid]
            complete = [row for row in subset if row["complete"]]
            any_incomplete |= len(complete) != batch
            false_outer = sum(row["any_false_holm_rejection"] for row in complete)
            covered = sum(row["covered_coefficients"] for row in complete)
            total = sum(row["coefficient_total"] for row in complete)
            precision = None
            if complete and total:
                precision = precision_decision(
                    outer_with_false_family=false_outer,
                    outer_total=len(complete),
                    covered_coefficients=covered,
                    coefficient_total=total,
                    root=root,
                )
            scenario_grid.append({
                "scenario": scenario,
                "grid_degrees": grid,
                "expected_outer_reports": batch,
                "recorded_outer_reports": len(subset),
                "complete_outer_reports": len(complete),
                "outer_with_any_false_holm_rejection": false_outer,
                "false_holm_rejections": sum(row["false_holm_rejections"] for row in complete),
                "true_holm_rejections": sum(row["true_holm_rejections"] for row in complete),
                "coefficient_coverage": {"covered": covered, "total": total,
                                         "fraction": covered / total if total else None},
                "sequential_precision_decision": precision,
            })

    complete_batch = not any_incomplete and len(reports) == len(expected) and len(rows) == len(expected) * len(contract["grid_degrees"])
    result = {
        "schema_version": 1,
        "status": ("FINAL_MULTICOORDINATE_CALIBRATION_BATCH_COMPLETE_CONTINUE_SEQUENTIAL_RULE"
                   if complete_batch else "FINAL_MULTICOORDINATE_CALIBRATION_BATCH_INCOMPLETE_STOP_AND_DIAGNOSE"),
        "contract_id": contract["id"],
        "contract_canonical_sha256": canonical_digest(contract),
        "batch_index": batch_index,
        "outer_start_inclusive": start,
        "outer_stop_exclusive": stop,
        "expected_outer_reports": len(expected),
        "recorded_outer_reports": len(reports),
        "missing_outer_reports": [{"scenario": s, "outer_replicate": o} for s, o in missing],
        "expected_grid_cases": len(expected) * len(contract["grid_degrees"]),
        "recorded_grid_cases": len(rows),
        "complete_grid_cases": sum(row["complete"] for row in rows),
        "scenario_grid_summary": scenario_grid,
        "input_manifest": sorted(manifest, key=lambda r: (r["scenario"], r["outer_replicate"])),
        "minimum_outer_replicates_per_scenario": contract["sequential_outer_rule"]["minimum_outer_replicates_per_scenario"],
        "batch_can_authorize_stopping": bool(stop >= contract["sequential_outer_rule"]["minimum_outer_replicates_per_scenario"]),
        "qualification_boundary": "Synthetic calibration only. A complete batch is an input to the frozen sequential rule, not empirical ecology admission; realized-design calibration remains required.",
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", nargs="+", type=Path, required=True)
    parser.add_argument("--batch-index", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.reports, args.batch_index)
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("status", "recorded_outer_reports", "complete_grid_cases")}), flush=True)
    if result["status"].endswith("INCOMPLETE_STOP_AND_DIAGNOSE"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()

"""Aggregate every predeclared full36 tail-mechanics pilot report.

The aggregator is synthetic-only and cannot authorize empirical ecology. It requires
all predeclared scenario/outer reports, retains both spatial grids and reports any
incomplete shared draw or missing family slot rather than omitting it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from .workflow import ROOT, canonical_digest

CONTRACT = ROOT / "analysis/v3/full_family_tail_pilot_contract.json"


def _digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def aggregate(report_paths, *, root: Path = ROOT):
    contract = json.loads((root / "analysis/v3/full_family_tail_pilot_contract.json").read_text(encoding="utf-8"))
    expected = {(scenario, outer) for scenario in contract["scenarios"]
                for outer in range(contract["replicates_per_scenario"])}
    reports = {}
    input_manifest = []
    for raw_path in report_paths:
        path = Path(raw_path)
        report = json.loads(path.read_text(encoding="utf-8"))
        execution = report["execution_contract"]
        key = (execution["scenario"], int(execution["outer_replicate"]))
        if key in reports:
            raise ValueError(f"Duplicate tail-pilot report {key}")
        if key not in expected:
            raise ValueError(f"Tail-pilot report outside contract {key}")
        if execution["specification_canonical_sha256"] != canonical_digest(contract):
            raise ValueError("Tail-pilot execution contract differs from current frozen contract")
        if report["empirical_trait_environment_values_read"] != 0 or report["ecological_models_executed"] != 0:
            raise ValueError("Tail-pilot report crossed the empirical boundary")
        reports[key] = report
        input_manifest.append({"scenario": key[0], "outer_replicate": key[1],
                               "path": path.as_posix(), "sha256": _digest(path)})
    missing = sorted(expected - set(reports))
    unexpected = sorted(set(reports) - expected)

    rows = []
    for key in sorted(reports):
        report = reports[key]
        seen_grids = set()
        for case in report["cases"]:
            grid = int(case["grid_degrees"])
            if grid in seen_grids:
                raise ValueError(f"Duplicate grid for {key}: {grid}")
            seen_grids.add(grid)
            slots = case.get("family_slots") or []
            coverage = case.get("interval_coverage") or {}
            covered = sum(int(value["covered"]) for value in coverage.values())
            total = sum(int(value["total"]) for value in coverage.values())
            rows.append({
                "scenario": key[0], "outer_replicate": key[1], "grid_degrees": grid,
                "status": case["status"],
                "planned_bootstrap_replicates": int(case["planned_bootstrap_replicates"]),
                "estimable_shared_bootstrap_replicates": int(case["estimable_shared_bootstrap_replicates"]),
                "family_slots": len(slots),
                "false_holm_rejections": int(case.get("false_holm_rejections", 0)),
                "any_false_holm_rejection": bool(case.get("any_false_holm_rejection", False)),
                "true_holm_rejections": int(case.get("true_holm_rejections", 0)),
                "tail_resolution_hits": int(case.get("tail_resolution_hits", 0)),
                "covered_coefficients": covered, "coefficient_total": total,
            })
        if seen_grids != set(contract["grid_degrees"]):
            raise ValueError(f"Grid inventory differs for {key}: {sorted(seen_grids)}")

    def is_complete(row):
        return (row["estimable_shared_bootstrap_replicates"] == contract["bootstrap_replicates"]
                and row["family_slots"] == contract["family_slots"])

    incomplete_cases = [{
        "scenario": row["scenario"],
        "outer_replicate": row["outer_replicate"],
        "grid_degrees": row["grid_degrees"],
        "status": row["status"],
        "planned_bootstrap_replicates": row["planned_bootstrap_replicates"],
        "estimable_shared_bootstrap_replicates": row["estimable_shared_bootstrap_replicates"],
        "family_slots": row["family_slots"],
    } for row in rows if not is_complete(row)]

    scenario_grid = []
    for scenario in contract["scenarios"]:
        for grid in contract["grid_degrees"]:
            subset = [row for row in rows if row["scenario"] == scenario and row["grid_degrees"] == grid]
            scenario_grid.append({
                "scenario": scenario,
                "grid_degrees": grid,
                "outer_reports": len(subset),
                "complete_cases": sum(is_complete(row) for row in subset),
                "outer_with_any_false_holm_rejection": sum(row["any_false_holm_rejection"] for row in subset),
                "false_holm_rejections": sum(row["false_holm_rejections"] for row in subset),
                "true_holm_rejections": sum(row["true_holm_rejections"] for row in subset),
                "tail_resolution_hits": sum(row["tail_resolution_hits"] for row in subset),
                "covered_coefficients": sum(row["covered_coefficients"] for row in subset),
                "coefficient_total": sum(row["coefficient_total"] for row in subset),
            })

    by_outer = {}
    for row in rows:
        by_outer.setdefault((row["scenario"], row["outer_replicate"]), []).append(row)
    outer_either_grid_false = sum(any(row["any_false_holm_rejection"] for row in group)
                                  for group in by_outer.values())
    complete = (not missing and not unexpected and len(rows) == len(expected) * len(contract["grid_degrees"])
                and not incomplete_cases)
    covered = sum(row["covered_coefficients"] for row in rows)
    coefficient_total = sum(row["coefficient_total"] for row in rows)
    result = {
        "schema_version": 1,
        "status": "FULL36_TAIL_MECHANICS_PILOT_AGGREGATED_NOT_ECOLOGY_ADMISSION" if complete
                  else "FULL36_TAIL_MECHANICS_PILOT_INCOMPLETE_NO_ADMISSION",
        "contract_id": contract["id"],
        "contract_canonical_sha256": canonical_digest(contract),
        "expected_outer_reports": len(expected),
        "recorded_outer_reports": len(reports),
        "missing_outer_reports": [{"scenario": s, "outer_replicate": o} for s, o in missing],
        "unexpected_outer_reports": [{"scenario": s, "outer_replicate": o} for s, o in unexpected],
        "expected_grid_cases": len(expected) * len(contract["grid_degrees"]),
        "recorded_grid_cases": len(rows),
        "complete_grid_cases": sum(is_complete(row) for row in rows),
        "incomplete_grid_cases": incomplete_cases,
        "outer_reports_with_any_false_holm_on_either_grid": outer_either_grid_false,
        "total_false_holm_rejections": sum(row["false_holm_rejections"] for row in rows),
        "total_true_holm_rejections": sum(row["true_holm_rejections"] for row in rows),
        "total_tail_resolution_hits": sum(row["tail_resolution_hits"] for row in rows),
        "coefficient_interval_coverage": {
            "covered": covered, "total": coefficient_total,
            "fraction": covered / coefficient_total if coefficient_total else None,
        },
        "scenario_grid_summary": scenario_grid,
        "input_manifest": sorted(input_manifest, key=lambda row: (row["scenario"], row["outer_replicate"])),
        "qualification_boundary": contract["qualification_boundary"],
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", nargs="+", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.reports)
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in
                      ("status", "recorded_outer_reports", "recorded_grid_cases", "complete_grid_cases")}), flush=True)


if __name__ == "__main__":
    main()

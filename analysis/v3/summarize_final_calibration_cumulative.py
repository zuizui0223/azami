"""Cumulatively evaluate frozen final multicoordinate calibration batches.

Requires contiguous completed batches from outer 0 and their returned numerical
arrays. Monte Carlo precision uses independent outer datasets under the explicit
20260909 amendment. Historical count-only Wilson decisions remain descriptive.

Synthetic only. No threshold, seed, grid, family, truth, or empirical ecology
choice is made here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from .multicoordinate_calibration import precision_decision as historical_precision_decision
from .calibration_precision import amendment, precision_decision
from .verify_calibration_outers import verify_inventory
from .workflow import ROOT, canonical_digest

COMPLETE_BATCH_STATUS = "FINAL_MULTICOORDINATE_CALIBRATION_BATCH_COMPLETE_CONTINUE_SEQUENTIAL_RULE"
CONTINUE_STATUS = "FINAL_MULTICOORDINATE_CALIBRATION_CUMULATIVE_CONTINUE_NEXT_FROZEN_BATCH"
QUALIFIED_STATUS = "FINAL_MULTICOORDINATE_CALIBRATION_CUMULATIVE_ADMISSION_MET_STOP"
MAX_FAIL_STATUS = "FINAL_MULTICOORDINATE_CALIBRATION_CUMULATIVE_MAX_REACHED_WITHOUT_FULL_QUALIFICATION_STOP_AND_DIAGNOSE"
OUTERS_REQUIRED_STATUS = "FINAL_CALIBRATION_OUTER_REPORTS_REQUIRED_NO_QUALIFICATION"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _expected_manifest(contract: dict, start: int, stop: int):
    return {(scenario, outer) for scenario in contract["scenarios"] for outer in range(start, stop)}


def _validate_batch(report: dict, *, contract: dict, contract_digest: str, path: Path):
    rule = contract["sequential_outer_rule"]
    batch_size = int(rule["batch_size_per_scenario"])
    max_outer = int(rule["maximum_outer_replicates_per_scenario"])
    scenarios = tuple(contract["scenarios"])
    grids = tuple(int(v) for v in contract["grid_degrees"])

    if report.get("status") != COMPLETE_BATCH_STATUS:
        raise ValueError(f"Batch aggregate is not complete: {path}")
    if report.get("contract_id") != contract["id"]:
        raise ValueError(f"Batch contract id differs: {path}")
    if report.get("contract_canonical_sha256") != contract_digest:
        raise ValueError(f"Batch contract digest differs: {path}")
    if report.get("empirical_trait_environment_values_read") != 0:
        raise ValueError(f"Batch crossed empirical trait-environment boundary: {path}")
    if report.get("ecological_models_executed") != 0 or report.get("ecological_fitting_authorized") is not False:
        raise ValueError(f"Batch crossed ecological fitting boundary: {path}")

    batch_index = report.get("batch_index")
    if not isinstance(batch_index, int) or isinstance(batch_index, bool) or batch_index < 0:
        raise ValueError(f"Invalid batch index: {path}")
    start = batch_index * batch_size
    stop = start + batch_size
    if stop > max_outer:
        raise ValueError(f"Batch exceeds frozen maximum outer range: {path}")
    if report.get("outer_start_inclusive") != start or report.get("outer_stop_exclusive") != stop:
        raise ValueError(f"Batch outer range differs from frozen sequential rule: {path}")

    expected_outer = len(scenarios) * batch_size
    expected_grid = expected_outer * len(grids)
    if report.get("expected_outer_reports") != expected_outer or report.get("recorded_outer_reports") != expected_outer:
        raise ValueError(f"Batch outer inventory is incomplete: {path}")
    if report.get("missing_outer_reports") != []:
        raise ValueError(f"Batch has missing outer reports: {path}")
    if report.get("expected_grid_cases") != expected_grid or report.get("recorded_grid_cases") != expected_grid:
        raise ValueError(f"Batch grid inventory differs: {path}")
    if report.get("complete_grid_cases") != expected_grid:
        raise ValueError(f"Batch contains incomplete grid cases: {path}")

    manifest = report.get("input_manifest")
    if not isinstance(manifest, list):
        raise ValueError(f"Batch input manifest is absent: {path}")
    expected_manifest = _expected_manifest(contract, start, stop)
    observed_manifest = []
    for row in manifest:
        key = (row.get("scenario"), row.get("outer_replicate"))
        observed_manifest.append(key)
        if not isinstance(row.get("sha256"), str) or len(row["sha256"]) != 64:
            raise ValueError(f"Batch manifest lacks an outer SHA-256: {path}")
    if len(observed_manifest) != len(set(observed_manifest)):
        raise ValueError(f"Batch manifest contains duplicate outer reports: {path}")
    if set(observed_manifest) != expected_manifest:
        raise ValueError(f"Batch manifest is not the exact frozen outer range: {path}")

    summaries = report.get("scenario_grid_summary")
    if not isinstance(summaries, list):
        raise ValueError(f"Batch scenario-grid summary is absent: {path}")
    expected_cells = {(scenario, grid) for scenario in scenarios for grid in grids}
    cell_rows = {}
    for row in summaries:
        key = (row.get("scenario"), row.get("grid_degrees"))
        if key in cell_rows:
            raise ValueError(f"Duplicate scenario-grid summary: {path}")
        if key not in expected_cells:
            raise ValueError(f"Unexpected scenario-grid summary: {path}")
        if row.get("expected_outer_reports") != batch_size:
            raise ValueError(f"Scenario-grid expected count differs: {path}")
        if row.get("recorded_outer_reports") != batch_size or row.get("complete_outer_reports") != batch_size:
            raise ValueError(f"Scenario-grid batch is incomplete: {path}")
        false_outer = row.get("outer_with_any_false_holm_rejection")
        covered = (row.get("coefficient_coverage") or {}).get("covered")
        total = (row.get("coefficient_coverage") or {}).get("total")
        if not isinstance(false_outer, int) or isinstance(false_outer, bool) or not 0 <= false_outer <= batch_size:
            raise ValueError(f"Invalid false-family count: {path}")
        if not isinstance(covered, int) or isinstance(covered, bool) or not isinstance(total, int) or isinstance(total, bool):
            raise ValueError(f"Invalid coefficient coverage counts: {path}")
        if total <= 0 or not 0 <= covered <= total:
            raise ValueError(f"Coefficient coverage counts are out of range: {path}")
        cell_rows[key] = row
    if set(cell_rows) != expected_cells:
        raise ValueError(f"Batch lacks a frozen scenario-grid cell: {path}")

    return {
        "batch_index": batch_index,
        "start": start,
        "stop": stop,
        "cell_rows": cell_rows,
        "outer_manifest": manifest,
    }


def aggregate(batch_report_paths, *, outer_report_paths=None, root: Path = ROOT):
    contract = json.loads((root / "analysis/v3/final_module_calibration_contract.json").read_text(encoding="utf-8"))
    contract_digest = canonical_digest(contract)
    rule = contract["sequential_outer_rule"]
    batch_size = int(rule["batch_size_per_scenario"])
    min_outer = int(rule["minimum_outer_replicates_per_scenario"])
    max_outer = int(rule["maximum_outer_replicates_per_scenario"])
    scenarios = tuple(contract["scenarios"])
    grids = tuple(int(v) for v in contract["grid_degrees"])

    if not batch_report_paths:
        raise ValueError("At least one completed batch aggregate is required")

    batches = {}
    manifest = []
    for raw in batch_report_paths:
        path = Path(raw)
        report = json.loads(path.read_text(encoding="utf-8"))
        checked = _validate_batch(report, contract=contract, contract_digest=contract_digest, path=path)
        index = checked["batch_index"]
        if index in batches:
            raise ValueError(f"Duplicate batch index {index}")
        batches[index] = checked
        manifest.append({
            "batch_index": index,
            "path": path.as_posix(),
            "sha256": _sha256(path),
            "outer_start_inclusive": checked["start"],
            "outer_stop_exclusive": checked["stop"],
        })

    latest = max(batches)
    required_indices = list(range(latest + 1))
    if sorted(batches) != required_indices:
        raise ValueError("Cumulative calibration batches must be contiguous from batch 0")

    cumulative_outer = (latest + 1) * batch_size
    if cumulative_outer > max_outer:
        raise ValueError("Cumulative calibration exceeds frozen maximum outer count")

    precision_spec = amendment(root)
    outer_values, outer_manifest = None, []
    if outer_report_paths is not None:
        expected = [r for i in required_indices for r in batches[i]["outer_manifest"]]
        outer_values, outer_manifest = verify_inventory(outer_report_paths, expected, contract, root=root)

    cumulative_rows = []
    decisions = []
    for scenario in scenarios:
        for grid in grids:
            false_outer = 0
            false_holm = 0
            true_holm = 0
            covered = 0
            coefficient_total = 0
            for batch_index in required_indices:
                row = batches[batch_index]["cell_rows"][(scenario, grid)]
                false_outer += int(row["outer_with_any_false_holm_rejection"])
                false_holm += int(row.get("false_holm_rejections", 0))
                true_holm += int(row.get("true_holm_rejections", 0))
                coverage = row["coefficient_coverage"]
                covered += int(coverage["covered"])
                coefficient_total += int(coverage["total"])
            historical = historical_precision_decision(
                outer_with_false_family=false_outer,
                outer_total=cumulative_outer,
                covered_coefficients=covered,
                coefficient_total=coefficient_total,
                root=root,
            )
            decision = {"precision_satisfied": False, "admission_satisfied": False,
                        "status": "independent_outer_vectors_required"}
            if outer_values is not None:
                vectors = [outer_values[(scenario, i)][grid] for i in range(cumulative_outer)]
                if (sum(r["false_family"] for r in vectors) != false_outer
                        or sum(r["false_holm_rejections"] for r in vectors) != false_holm
                        or sum(r["true_holm_rejections"] for r in vectors) != true_holm
                        or sum(r["covered"] for r in vectors) != covered
                        or sum(r["total"] for r in vectors) != coefficient_total):
                    raise ValueError("Cumulative batch counts differ from verified independent outer data")
                decision = precision_decision([r["false_family"] for r in vectors],
                                              [r["coverage_fraction"] for r in vectors], root=root)
            decisions.append(decision)
            cumulative_rows.append({
                "scenario": scenario,
                "grid_degrees": grid,
                "cumulative_outer_reports": cumulative_outer,
                "outer_with_any_false_holm_rejection": false_outer,
                "false_holm_rejections": false_holm,
                "true_holm_rejections": true_holm,
                "coefficient_coverage": {
                    "covered": covered,
                    "total": coefficient_total,
                    "fraction": covered / coefficient_total,
                },
                "cumulative_precision_decision": decision,
                "historical_count_only_decision_not_valid_for_qualification": historical,
            })

    minimum_reached = cumulative_outer >= min_outer
    all_precision = minimum_reached and all(row["precision_satisfied"] for row in decisions)
    all_admission = minimum_reached and all(row["admission_satisfied"] for row in decisions)
    synthetic_qualified = all_precision and all_admission

    if outer_values is None:
        status = OUTERS_REQUIRED_STATUS
        next_batch_index = None
    elif synthetic_qualified:
        status = QUALIFIED_STATUS
        next_batch_index = None
    elif cumulative_outer >= max_outer:
        status = MAX_FAIL_STATUS
        next_batch_index = None
    else:
        status = CONTINUE_STATUS
        next_batch_index = latest + 1

    next_outer_range = None
    if next_batch_index is not None:
        next_start = next_batch_index * batch_size
        next_outer_range = {
            "batch_index": next_batch_index,
            "outer_start_inclusive": next_start,
            "outer_stop_exclusive": next_start + batch_size,
        }

    return {
        "schema_version": 1,
        "status": status,
        "contract_id": contract["id"],
        "contract_canonical_sha256": contract_digest,
        "included_batch_indices": required_indices,
        "completed_batches": len(required_indices),
        "cumulative_outer_replicates_per_scenario": cumulative_outer,
        "minimum_outer_replicates_per_scenario": min_outer,
        "maximum_outer_replicates_per_scenario": max_outer,
        "minimum_reached": minimum_reached,
        "all_scenario_grid_precision_satisfied": all_precision,
        "all_scenario_grid_admission_satisfied": all_admission,
        "synthetic_calibration_qualified_for_realized_design_gate": synthetic_qualified,
        "next_frozen_batch": next_outer_range,
        "scenario_grid_cumulative_summary": cumulative_rows,
        "input_batch_manifest": sorted(manifest, key=lambda row: row["batch_index"]),
        "input_outer_manifest": outer_manifest,
        "independent_outer_arrays_verified": outer_values is not None,
        "precision_amendment_id": precision_spec["id"],
        "precision_amendment_canonical_sha256": canonical_digest(precision_spec),
        "decision_rule": "Use independent outer-dataset coverage fractions and false-family indicators with 16 simultaneous, time-uniform confidence sequences. All scenarios/grids must meet unchanged precision/admission thresholds at n>=100; otherwise continue by 25 until n=400. Counts without verified outer arrays never qualify or authorize another batch.",
        "uncertainty_method_changed_by_explicit_amendment": True,
        "batch_local_stopping_decisions_used": False,
        "thresholds_changed": False,
        "seed_changed": False,
        "grid_changed": False,
        "family_changed": False,
        "synthetic_truth_changed": False,
        "qualification_boundary": "A successful cumulative synthetic calibration advances only to the frozen realized-design synthetic-response gate. It never authorizes empirical ecological fitting by itself.",
        "empirical_trait_environment_values_read": 0,
        "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", nargs="+", type=Path, required=True)
    parser.add_argument("--outer-reports", nargs="+", type=Path,
                        help="Returned original outer reports with sibling NPZ arrays and draw records")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.batches, outer_report_paths=args.outer_reports)
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "cumulative_outer_replicates_per_scenario": result["cumulative_outer_replicates_per_scenario"],
        "synthetic_calibration_qualified_for_realized_design_gate": result["synthetic_calibration_qualified_for_realized_design_gate"],
        "next_frozen_batch": result["next_frozen_batch"],
    }), flush=True)
    if result["status"] in (MAX_FAIL_STATUS, OUTERS_REQUIRED_STATUS):
        raise SystemExit(2)


if __name__ == "__main__":
    main()

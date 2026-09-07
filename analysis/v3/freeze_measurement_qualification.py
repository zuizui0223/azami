"""Freeze Chapter 1 v3 image-only measurement qualification before ecology.

The freeze consumes the predeclared qualification contract, the completed
128-photo large-vs-original technical report, and an explicit endpoint decision.
It does not read environmental values, trait-environment coefficients or P values.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


FREEZE_STATUS = "V3_MEASUREMENT_QUALIFICATION_FROZEN_BEFORE_TRAIT_ENVIRONMENT_JOIN"
DECISION_STATUS = "phenotype_blind_measurement_qualification_decision"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def endpoint_universe(contract: dict) -> set[str]:
    endpoints: list[str] = []
    for module in contract["module_rules"].values():
        endpoints.extend(map(str, module["endpoints"]))
    if len(endpoints) != int(contract.get("endpoint_partition_count", -1)):
        raise ValueError("Measurement contract endpoint count differs from partition")
    if len(endpoints) != len(set(endpoints)):
        raise ValueError("Measurement contract assigns an endpoint to multiple modules")
    return set(endpoints)


def expected_route(status: str) -> set[str]:
    if status == "NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE":
        return {"large_allowed_with_qc_and_uncertainty", "blocked_by_other_measurement_rule"}
    if status == "ORIGINAL_MAY_ADD_INFORMATION_BEYOND_LARGE":
        return {"stream_original_required", "blocked_by_other_measurement_rule"}
    if status == "INSUFFICIENT_FOLLOWUP_INFORMATION":
        return {"blocked_pending_additional_gate"}
    raise ValueError(f"Unknown resolution followup status: {status}")


def freeze(contract_path: Path, followup_report_path: Path, decision_path: Path, out_path: Path) -> dict:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    report = json.loads(followup_report_path.read_text(encoding="utf-8"))
    decision = json.loads(decision_path.read_text(encoding="utf-8"))

    if contract.get("status") != "fixed_before_128_photo_followup_results":
        raise ValueError("Measurement qualification contract is not the pre-result frozen contract")
    if contract.get("trait_environment_results_inspected_for_this_decision") is not False:
        raise ValueError("Measurement contract is not phenotype/environment-outcome blind")
    if report.get("status") != "LARGE_ORIGINAL_STREAM_FOLLOWUP_COMPLETE":
        raise ValueError("128-photo large-original technical followup is incomplete")
    if int(report.get("photos_requested", -1)) != 128 or int(report.get("photos_selected", -1)) != 128:
        raise ValueError("Measurement freeze requires the declared 128-photo followup")
    if int(report.get("source_images_persisted", -1)) != 0:
        raise ValueError("Followup persisted source images unexpectedly")
    if report.get("full_original_archive_required_by_design") is not False:
        raise ValueError("Followup unexpectedly requires a full image archive")
    if report.get("ecological_models_executed") is not False:
        raise ValueError("Followup is not ecology-blind")

    if decision.get("status") != DECISION_STATUS:
        raise ValueError("Measurement decision has wrong status")
    if decision.get("trait_environment_results_inspected") is not False:
        raise ValueError("Measurement decision must state trait_environment_results_inspected=false")
    if decision.get("contract_id") != contract.get("contract_id"):
        raise ValueError("Measurement decision and contract IDs differ")

    endpoints = endpoint_universe(contract)
    followup = report.get("endpoint_summary") or {}
    if set(followup) != endpoints:
        missing = sorted(endpoints - set(followup))
        extra = sorted(set(followup) - endpoints)
        raise ValueError(f"128-photo report endpoint universe differs; missing={missing}, extra={extra}")

    entries = decision.get("endpoints")
    if not isinstance(entries, list) or len(entries) != len(endpoints):
        raise ValueError("Measurement decision must contain one row per endpoint")
    by_endpoint: dict[str, dict] = {}
    for entry in entries:
        endpoint = str(entry.get("endpoint_id") or "")
        if endpoint in by_endpoint or endpoint not in endpoints:
            raise ValueError(f"Duplicate or unknown measurement decision endpoint: {endpoint}")
        status = str(entry.get("resolution_status") or "")
        if status != str(followup[endpoint].get("status")):
            raise ValueError(f"Decision resolution status disagrees for {endpoint}")
        route = str(entry.get("ecological_route") or "")
        if route not in expected_route(status):
            raise ValueError(f"Ecological route {route} is not allowed for {endpoint} with status {status}")
        rationale = str(entry.get("rationale") or "").strip()
        if not rationale:
            raise ValueError(f"Missing measurement rationale for {endpoint}")
        by_endpoint[endpoint] = entry
    if set(by_endpoint) != endpoints:
        raise ValueError("Measurement decision does not exactly partition endpoints")

    orientation = by_endpoint["orientation_image_vertical_angle"]
    if orientation.get("bbox_uncertainty_required") is not True:
        raise ValueError("Orientation qualification must retain upstream bbox uncertainty")
    if decision.get("raw_values_retained_when_inferentially_ineligible") is not True:
        raise ValueError("Measurement decision must retain raw values when inference is blocked")

    route_counts: dict[str, int] = {}
    resolution_counts: dict[str, int] = {}
    for entry in entries:
        route_counts[entry["ecological_route"]] = route_counts.get(entry["ecological_route"], 0) + 1
        resolution_counts[entry["resolution_status"]] = resolution_counts.get(entry["resolution_status"], 0) + 1

    receipt = {
        "status": FREEZE_STATUS,
        "contract_id": contract["contract_id"],
        "followup": {
            "status": report["status"],
            "photos_requested": report["photos_requested"],
            "photos_selected": report["photos_selected"],
            "successful_pairs": report.get("successful_pairs"),
            "matched_head_pairs": report.get("matched_head_pairs"),
            "endpoint_status_counts": report.get("endpoint_status_counts"),
            "report_sha256": sha256_file(followup_report_path),
        },
        "decision_sha256": sha256_file(decision_path),
        "contract_sha256": sha256_file(contract_path),
        "endpoint_count": len(endpoints),
        "resolution_status_counts": resolution_counts,
        "ecological_route_counts": route_counts,
        "orientation": {
            "resolution_status": orientation["resolution_status"],
            "ecological_route": orientation["ecological_route"],
            "bbox_uncertainty_required": True,
            "bbox_uncertainty_source": "predeclared four +/-5% bounding-box shifts from cached 14-condition perturbation design"
        },
        "raw_values_retained_when_inferentially_ineligible": True,
        "trait_environment_results_inspected": False,
        "ecological_models_executed": 0,
        "trait_environment_join_authorized_by_measurement_gate": True,
        "claim_boundary": "Qualification addresses image-pipeline information retention and synthetic technical sensitivity, not detector truth or physical botanical accuracy."
    }
    if out_path.exists():
        raise ValueError("Measurement freeze output already exists")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=Path("analysis/v3/measurement_qualification_contract.json"))
    parser.add_argument("--followup-report", type=Path, required=True)
    parser.add_argument("--decision", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    receipt = freeze(args.contract, args.followup_report, args.decision, args.out)
    print(json.dumps({"status": receipt["status"], "routes": receipt["ecological_route_counts"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

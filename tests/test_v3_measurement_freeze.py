import json
from pathlib import Path

import pytest

from analysis.v3.freeze_measurement_qualification import FREEZE_STATUS, freeze


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "measurement_qualification_contract.json"


def endpoints():
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    return [e for module in c["module_rules"].values() for e in module["endpoints"]]


def write_followup(path: Path, status="NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE"):
    report = {
        "status": "LARGE_ORIGINAL_STREAM_FOLLOWUP_COMPLETE",
        "photos_requested": 128,
        "photos_selected": 128,
        "successful_pairs": 128,
        "matched_head_pairs": 160,
        "source_images_persisted": 0,
        "full_original_archive_required_by_design": False,
        "ecological_models_executed": False,
        "endpoint_status_counts": {status: 27},
        "endpoint_summary": {e: {"status": status} for e in endpoints()},
    }
    path.write_text(json.dumps(report), encoding="utf-8")


def write_decision(path: Path, status="NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE"):
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    decision = {
        "status": "phenotype_blind_measurement_qualification_decision",
        "contract_id": contract["contract_id"],
        "trait_environment_results_inspected": False,
        "raw_values_retained_when_inferentially_ineligible": True,
        "endpoints": [
            {
                "endpoint_id": e,
                "resolution_status": status,
                "ecological_route": "large_allowed_with_qc_and_uncertainty",
                "bbox_uncertainty_required": e == "orientation_image_vertical_angle",
                "rationale": "Synthetic stable-route test.",
            }
            for e in endpoints()
        ],
    }
    path.write_text(json.dumps(decision), encoding="utf-8")


def test_measurement_freeze_requires_exact_27_endpoint_partition(tmp_path: Path):
    report = tmp_path / "report.json"
    decision = tmp_path / "decision.json"
    out = tmp_path / "freeze.json"
    write_followup(report)
    write_decision(decision)
    receipt = freeze(CONTRACT, report, decision, out)
    assert receipt["status"] == FREEZE_STATUS
    assert receipt["endpoint_count"] == 27
    assert receipt["orientation"]["bbox_uncertainty_required"] is True
    assert receipt["trait_environment_results_inspected"] is False
    assert out.is_file()


def test_original_required_endpoint_cannot_use_large_route(tmp_path: Path):
    status = "ORIGINAL_MAY_ADD_INFORMATION_BEYOND_LARGE"
    report = tmp_path / "report.json"
    decision = tmp_path / "decision.json"
    write_followup(report, status=status)
    write_decision(decision, status=status)
    data = json.loads(decision.read_text())
    for row in data["endpoints"]:
        row["ecological_route"] = "large_allowed_with_qc_and_uncertainty"
    decision.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="not allowed"):
        freeze(CONTRACT, report, decision, tmp_path / "freeze.json")


def test_insufficient_endpoint_must_remain_blocked(tmp_path: Path):
    status = "INSUFFICIENT_FOLLOWUP_INFORMATION"
    report = tmp_path / "report.json"
    decision = tmp_path / "decision.json"
    write_followup(report, status=status)
    write_decision(decision, status=status)
    data = json.loads(decision.read_text())
    for row in data["endpoints"]:
        row["ecological_route"] = "blocked_pending_additional_gate"
    decision.write_text(json.dumps(data), encoding="utf-8")
    receipt = freeze(CONTRACT, report, decision, tmp_path / "freeze.json")
    assert receipt["ecological_route_counts"] == {"blocked_pending_additional_gate": 27}

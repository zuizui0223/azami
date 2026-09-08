import hashlib
import json
from pathlib import Path

import pytest

from analysis.v3.freeze_measurement_qualification import FREEZE_STATUS, followup_thresholds, freeze


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "measurement_qualification_contract.json"


def endpoints():
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    return [e for module in c["module_rules"].values() for e in module["endpoints"]]


def write_followup(path: Path, status="NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE"):
    common = 3 if status == "INSUFFICIENT_FOLLOWUP_INFORMATION" else 40
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
        "decision_thresholds": followup_thresholds(),
        "endpoint_summary": {e: {"status": status, "both_usable": common,
                                  "numeric_decidable": common >= 30} for e in endpoints()},
    }
    path.write_text(json.dumps(report), encoding="utf-8")


def write_decision(path: Path, status="NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE"):
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    decision = {
        "status": "phenotype_blind_measurement_qualification_decision",
        "contract_id": contract["contract_id"],
        "trait_environment_results_inspected": False,
        "raw_values_retained_when_inferentially_ineligible": True,
        "followup_provenance": {"report_sha256": hashlib.sha256(path.with_name("report.json").read_bytes()).hexdigest()},
        "endpoints": [
            {
                "endpoint_id": e,
                "module": name,
                "resolution_status": status,
                "ecological_route": "large_allowed_with_qc_and_uncertainty",
                "bbox_uncertainty_required": "bbox_rule" in module,
                "rationale": "Synthetic stable-route test.",
            }
            for name, module in contract["module_rules"].items() for e in module["endpoints"]
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


def synthetic_inputs(tmp_path):
    report, decision, out = (tmp_path / name for name in ("report.json", "decision.json", "freeze.json"))
    write_followup(report)
    write_decision(decision)
    return report, decision, out


def replace_report_and_rebind(report, decision, data):
    """Declare an intentional synthetic evidence change for structural tests."""
    report.write_text(json.dumps(data), encoding="utf-8")
    chosen = json.loads(decision.read_text(encoding="utf-8"))
    chosen["followup_provenance"]["report_sha256"] = hashlib.sha256(report.read_bytes()).hexdigest()
    decision.write_text(json.dumps(chosen), encoding="utf-8")


def test_followup_thresholds_come_from_existing_predeclared_analysis():
    assert followup_thresholds() == {
        "minimum_common_usable_head_pairs": 30,
        "minimum_discordant_eligibility_pairs": 15,
        "original_gain_share_among_discordant": .80,
        "eligibility_exact_p_max": .01,
        "rank_correlation_min": .90,
        "median_abs_delta_over_original_iqr_max": .15,
    }


def test_same_status_substituted_report_cannot_be_frozen(tmp_path):
    report, decision, out = synthetic_inputs(tmp_path)
    report.write_text(report.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 differs"):
        freeze(CONTRACT, report, decision, out)
    assert not out.exists()


@pytest.mark.parametrize("declared_sha", [None, "", "g" * 64, "a" * 63])
def test_missing_or_invalid_evidence_hash_fails_closed(tmp_path, declared_sha):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(decision.read_text(encoding="utf-8"))
    data["followup_provenance"]["report_sha256"] = declared_sha
    decision.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="exact followup report SHA-256"):
        freeze(CONTRACT, report, decision, out)
    assert not out.exists()


@pytest.mark.parametrize("module_name", ["orientation", "gross_shape"])
def test_contract_bbox_rule_cannot_be_disabled(tmp_path, module_name):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(decision.read_text(encoding="utf-8"))
    next(row for row in data["endpoints"] if row["module"] == module_name)["bbox_uncertainty_required"] = False
    decision.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="requires bbox uncertainty"):
        freeze(CONTRACT, report, decision, out)
    assert not out.exists()


def test_relabelling_endpoint_cannot_bypass_its_module(tmp_path):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(decision.read_text(encoding="utf-8"))
    data["endpoints"][0]["module"] = "visible_colour"
    decision.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="module differs"):
        freeze(CONTRACT, report, decision, out)


@pytest.mark.parametrize("common", [9, 29])
def test_fine_geometry_cannot_be_promoted_from_insufficient_common_pairs(tmp_path, common):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(report.read_text(encoding="utf-8"))
    endpoint = "bract_projection_asymmetry"
    data["endpoint_summary"][endpoint].update(both_usable=common, numeric_decidable=False)
    replace_report_and_rebind(report, decision, data)
    chosen = json.loads(decision.read_text(encoding="utf-8"))
    next(row for row in chosen["endpoints"] if row["endpoint_id"] == endpoint)["followup_common_usable"] = 999
    decision.write_text(json.dumps(chosen), encoding="utf-8")
    with pytest.raises(ValueError, match="requires a blocked route"):
        freeze(CONTRACT, report, decision, out)
    assert not out.exists()


def test_fine_geometry_retains_blocked_route_below_existing_boundary(tmp_path):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(report.read_text(encoding="utf-8"))
    chosen = json.loads(decision.read_text(encoding="utf-8"))
    for row in chosen["endpoints"]:
        if row["module"] == "fine_architecture_and_surface":
            data["endpoint_summary"][row["endpoint_id"]].update(both_usable=9, numeric_decidable=False)
            row["ecological_route"] = "blocked_by_other_measurement_rule"
    decision.write_text(json.dumps(chosen), encoding="utf-8")
    replace_report_and_rebind(report, decision, data)
    receipt = freeze(CONTRACT, report, decision, out)
    assert receipt["ecological_route_counts"] == {
        "large_allowed_with_qc_and_uncertainty": 14, "blocked_by_other_measurement_rule": 13}


@pytest.mark.parametrize("common,numeric_decidable", [(None, False), (-1, False), (True, False), (9, True), (30, False)])
def test_fine_geometry_requires_consistent_numeric_support(tmp_path, common, numeric_decidable):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(report.read_text(encoding="utf-8"))
    data["endpoint_summary"]["bract_projection_asymmetry"].update(both_usable=common, numeric_decidable=numeric_decidable)
    replace_report_and_rebind(report, decision, data)
    with pytest.raises(ValueError, match="common-usable count|Common-usable count"):
        freeze(CONTRACT, report, decision, out)


def test_followup_numeric_threshold_cannot_be_lowered_after_results(tmp_path):
    report, decision, out = synthetic_inputs(tmp_path)
    data = json.loads(report.read_text(encoding="utf-8"))
    data["decision_thresholds"]["minimum_common_usable_head_pairs"] = 9
    replace_report_and_rebind(report, decision, data)
    with pytest.raises(ValueError, match="thresholds differ"):
        freeze(CONTRACT, report, decision, out)

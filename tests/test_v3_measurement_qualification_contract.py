import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "measurement_qualification_contract.json"


def load_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_contract_was_fixed_before_followup_results_and_environment_outcomes():
    c = load_contract()
    assert c["status"] == "fixed_before_128_photo_followup_results"
    assert c["trait_environment_results_inspected_for_this_decision"] is False
    assert c["environment_representation_already_frozen"] is True
    assert c["measurement_freeze_requirement"]["required_before_trait_environment_join"] is True


def test_all_27_endpoints_are_partitioned_exactly_once():
    c = load_contract()
    endpoints = []
    for module in c["module_rules"].values():
        endpoints.extend(module["endpoints"])
    assert len(endpoints) == 27
    assert len(set(endpoints)) == 27
    assert c["endpoint_partition_count"] == 27
    assert "visible_floret_fraction" in endpoints
    assert "orientation_image_vertical_angle" in endpoints


def test_resolution_statuses_have_predeclared_actions():
    rules = load_contract()["resolution_followup_rules"]
    assert set(rules) == {
        "NO_FOLLOWUP_EVIDENCE_ORIGINAL_REQUIRED_BEYOND_LARGE",
        "ORIGINAL_MAY_ADD_INFORMATION_BEYOND_LARGE",
        "INSUFFICIENT_FOLLOWUP_INFORMATION",
    }
    assert "streamed original-resolution processing" in rules["ORIGINAL_MAY_ADD_INFORMATION_BEYOND_LARGE"]
    assert "outside primary ecological inference" in rules["INSUFFICIENT_FOLLOWUP_INFORMATION"]


def test_orientation_bbox_uncertainty_is_upstream():
    orientation = load_contract()["module_rules"]["orientation"]
    assert orientation["ecological_priority"] == "Tier_1"
    assert "+/-5% shifts" in orientation["bbox_rule"]
    assert "before environmental join" in orientation["bbox_rule"]

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "environment_exposure_contract.json"
PRIORITY = ROOT / "analysis" / "v3" / "environment_selection_priority_20260907.json"


def load_priority():
    return json.loads(PRIORITY.read_text(encoding="utf-8"))


def test_priority_was_fixed_without_complete_diagnostic_or_trait_outcomes():
    p = load_priority()
    assert p["status"] == "phenotype_blind_selection_priority_fixed_before_complete_15_candidate_diagnostics"
    assert p["complete_15_candidate_diagnostic_results_inspected"] is False
    assert p["partial_13_candidate_run_used_for_selection"] is False
    assert set(p["selection_basis"]) == {
        "biological_proximity",
        "coverage",
        "environment_only_redundancy",
    }
    assert "Trait coefficients, P values" in p["general_rule"]


def test_priority_covers_exactly_all_contract_candidates():
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    priority = load_priority()
    candidates = {row["id"] for row in contract["monthly_candidates"]}
    candidates |= {row["id"] for row in contract["broader_climate_representations"]}
    assert len(candidates) == 15
    assert set(priority["candidate_roles_before_diagnostics"]) == candidates
    assert priority["contract_id"] == contract["contract_id"]


def test_direct_process_priority_matches_frozen_hypotheses():
    roles = load_priority()["candidate_roles_before_diagnostics"]
    assert roles["pr_month"]["priority"] == "primary_direct"
    assert roles["rsds_month"]["priority"] == "primary_direct"
    assert roles["vpd_month"]["priority"] == "primary_direct"
    assert roles["tasmax_month"]["priority"] == "primary_direct_thermal_candidate"
    assert roles["sfcWind_month"]["priority"] == "secondary_direct"
    assert roles["pet_month"]["priority"] == "derived_alternative"
    assert roles["cmi_month"]["priority"] == "derived_alternative"
    assert roles["NPP"]["priority"] == "distal_context"


def test_redundancy_priority_is_predeclared_not_p_value_selected():
    p = load_priority()
    rules = " ".join(p["decision_rules"])
    assert "absolute correlation >= 0.80" in rules
    bio12 = p["candidate_roles_before_diagnostics"]["BIO12"]
    assert "pr_month" in bio12["redundancy_preference"] or "observation-month precipitation" in bio12["redundancy_preference"]
    assert "trait outcomes" in rules
    # Semantic boundary: a broad or distal candidate is not promoted merely
    # because it escapes the redundancy flag; biological proximity remains required.
    assert "Non-redundancy alone does not promote" in rules or "non-redundancy alone" in rules.casefold()

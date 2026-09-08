import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "ecological_analysis_contract.json"


def load_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_ecological_fitting_is_paused_until_implementation_is_verified():
    c = load_contract()
    assert c["status"] == "implementation_corrections_pending_before_trait_join"
    assert hold_review(c)["ecological_fitting_authorized"] is False
    hold = c["design_hold"]
    assert hold["ecological_fitting_paused"] is True
    assert hold["environment_gate_completed"] is True
    assert hold["measurement_freeze_completed"] is True
    assert hold["trait_join_requires_source_cohort_freeze_receipt"] is True
    assert hold["trait_results_must_not_be_inspected_during_source_cohort_or_support_selection"] is True
    assert "analysis/v3/ecological_source_cohort_contract.json" in hold["next_required_artifacts"]
    assert any("phenotype-blind" in item for item in hold["next_required_artifacts"])
    assert any("taxon-support" in item for item in hold["next_required_artifacts"])


def hold_review(contract):
    return json.loads((ROOT / contract["design_hold"]["implementation_review"]).read_text(encoding="utf-8"))


def test_numerical_review_does_not_silently_promote_ecological_execution():
    review = hold_review(load_contract())
    assert review["full_original_stream_authorized"] is False
    assert review["v3_trait_environment_models_executed_in_this_review"] == 0
    assert review["historical_results_rewritten"] is False
    assert "cross-taxon" in review["numerical_correction"]["solver"]
    assert len(review["required_before_ecological_execution"]) >= 4


def test_completed_environment_and_measurement_gates_are_recorded():
    c = load_contract()
    done = c["completed_design_gates"]
    assert done["biological_hypotheses"] == "analysis/v3/capitulum_abiotic_hypotheses_v3.md"
    assert done["full_source_native_environment_upstream"] == "reproducibility/v3_environment_design_upstream_20260907.json"
    assert done["environment_decision"] == "analysis/v3/environment_representation_decision_20260907.json"
    assert done["environment_freeze_receipt"] == "reproducibility/v3_environment_representation_freeze_20260907.json"
    assert done["environment_freeze_status"] == "V3_ENVIRONMENT_REPRESENTATION_FROZEN_BEFORE_TRAIT_JOIN"
    assert done["measurement_contract"] == "analysis/v3/measurement_qualification_contract.json"
    assert done["measurement_decision"] == "analysis/v3/measurement_qualification_decision_20260908.json"
    assert done["measurement_freeze_receipt"] == "reproducibility/v3_measurement_qualification_freeze_20260908.json"
    assert done["measurement_freeze_status"] == "V3_MEASUREMENT_QUALIFICATION_FROZEN_BEFORE_ECOLOGY"


def test_preserved_scientific_requirements_survive_redesign():
    c = load_contract()
    p = c["preserved_requirements"]
    assert p["primary_ecological_range_scope"] == "native_only"
    assert p["taxonomic_join_required"] is True
    assert p["environment_representation_freeze_required"] is True
    assert p["hierarchy"] == ["head", "photo", "observation", "taxon"]
    assert p["dominant_taxa_solution"] == "hierarchical_partial_pooling_random_slopes_with_taxon_level_summaries"
    assert p["measurement_uncertainty_upstream"] is True
    assert p["orientation_stability_required_before_environment_join"] is True
    assert p["measurement_values_are_retained_even_when_inferentially_ineligible"] is True
    assert p["no_taxon_row_cap_in_final_ecological_source_cohort"] is True
    assert p["no_spatial_thinning_in_master_ecological_source_cohort"] is True
    assert set(p["environment_collinearity_diagnostics_required"]) == {"correlation_matrix", "matrix_rank", "condition_number", "VIF"}
    assert set(p["environment_selection_basis"]) == {"biological_proximity", "coverage", "environment_only_redundancy"}
    assert p["legacy_v2_candidates_need_not_survive"] is True


def test_source_cohort_gate_is_upstream_not_result_defense():
    c = load_contract()
    reason = c["design_hold"]["reason"]
    assert "Environment and measurement qualification are frozen" in reason
    assert "phenotype-blind ecological source cohort" in reason
    assert "before original-image trait execution" in reason
    assert "environmental exposures" in reason


def test_abiotic_only_scope_is_retained():
    c = load_contract()
    assert "outside the current estimand" in c["scope"]["biotic_scope"]
    assert "adaptation" in c["scope"]["causal_ceiling"]


def test_design_ledgers_and_freeze_artifacts_exist():
    assert (ROOT / "analysis" / "v3" / "capitulum_abiotic_hypotheses_v3.md").is_file()
    assert (ROOT / "analysis" / "v3" / "design_rationale.md").is_file()
    assert (ROOT / "analysis" / "v3" / "environment_exposure_contract.json").is_file()
    assert (ROOT / "analysis" / "v3" / "environment_representation_decision_20260907.json").is_file()
    assert (ROOT / "reproducibility" / "v3_environment_representation_freeze_20260907.json").is_file()
    assert (ROOT / "analysis" / "v3" / "measurement_qualification_contract.json").is_file()
    assert (ROOT / "analysis" / "v3" / "measurement_qualification_decision_20260908.json").is_file()
    assert (ROOT / "reproducibility" / "v3_measurement_qualification_freeze_20260908.json").is_file()
    assert (ROOT / "analysis" / "v3" / "ecological_source_cohort_contract.json").is_file()

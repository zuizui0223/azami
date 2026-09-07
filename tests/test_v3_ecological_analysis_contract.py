import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "ecological_analysis_contract.json"


def load_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_ecological_fitting_is_paused_until_environment_representation_is_frozen():
    c = load_contract()
    assert c["status"] == "environment_representation_freeze_pending_before_v3_ecological_fitting"
    hold = c["design_hold"]
    assert hold["ecological_fitting_paused"] is True
    assert hold["completed_upstream_artifact"] == "analysis/v3/capitulum_abiotic_hypotheses_v3.md"
    assert hold["trait_join_requires_environment_freeze_receipt"] is True
    assert hold["environment_freeze_receipt_status"] == "V3_ENVIRONMENT_REPRESENTATION_FROZEN_BEFORE_TRAIT_JOIN"
    assert hold["trait_results_must_not_be_inspected_during_environment_selection"] is True
    assert "analysis/v3/environment_representation_decision_20260907.json" in hold["next_required_artifacts"]
    assert "reproducibility/v3_environment_representation_freeze_20260907.json" in hold["next_required_artifacts"]


def test_preserved_scientific_requirements_survive_redesign():
    c = load_contract()
    p = c["preserved_requirements"]
    assert p["primary_ecological_range_scope"] == "native_only"
    assert p["taxonomic_join_required"] is True
    assert p["environment_representation_freeze_required"] is True
    assert p["hierarchy"] == ["head", "photo", "observation", "taxon"]
    assert p["dominant_taxa_solution"] == "hierarchical_partial_pooling_random_slopes_with_taxon_level_summaries"
    assert p["measurement_uncertainty_upstream"] is True
    assert set(p["environment_collinearity_diagnostics_required"]) == {"correlation_matrix", "matrix_rank", "condition_number", "VIF"}
    assert set(p["environment_selection_basis"]) == {"biological_proximity", "coverage", "environment_only_redundancy"}
    assert p["legacy_v2_candidates_need_not_survive"] is True


def test_no_climate_core_is_privileged_while_environment_freeze_is_open():
    c = load_contract()
    reason = c["design_hold"]["reason"]
    assert "phenotype-blind environment candidate matrix" in reason
    assert "abiotic_environment" not in c


def test_abiotic_only_scope_is_retained():
    c = load_contract()
    assert "outside the current estimand" in c["scope"]["biotic_scope"]
    assert "adaptation" in c["scope"]["causal_ceiling"]


def test_hypothesis_environment_and_supervisor_ledgers_exist():
    assert (ROOT / "analysis" / "v3" / "capitulum_abiotic_hypotheses_v3.md").is_file()
    assert (ROOT / "analysis" / "v3" / "supervisor_comment_integration_20260907.md").is_file()
    assert (ROOT / "analysis" / "v3" / "environment_exposure_contract.json").is_file()
    assert (ROOT / "analysis" / "v3" / "freeze_environment_representation.py").is_file()

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "ecological_analysis_contract.json"


def load_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_primary_ecology_is_native_range_and_abiotic_only():
    c = load_contract()
    assert c["sampling_and_replication"]["primary_ecological_range_scope"] == "native_only"
    assert "outside the current estimand" in c["scope"]["biotic_scope"]
    assert c["execution_gate"]["legacy_v2_candidates_need_not_survive"] is True


def test_dominant_taxa_are_handled_structurally_not_by_posthoc_deletion():
    c = load_contract()
    assert c["biological_scale"]["within_taxon"]["model"] == "hierarchical_partial_pooling_random_slopes"
    assert "raw pooled observation slope" in c["sampling_and_replication"]["dominant_taxa_rule"]
    assert c["biological_scale"]["among_taxon"]["taxon_weighting_rule"].startswith("One taxon summary")


def test_environmental_blocks_are_integrated_before_endpoint_decomposition():
    c = load_contract()
    order = c["abiotic_environment"]["primary_test_order"]
    assert order[0] == "fit_core_environment_representation"
    assert order[1] == "test_each_predeclared_additional_block_for_information_beyond_core"
    assert "correlation_matrix" in c["abiotic_environment"]["collinearity_before_fit"]["diagnostics"]
    assert c["abiotic_environment"]["warm_season_precipitation_alternative"]["predictor"] == "BIO18"


def test_measurement_uncertainty_precedes_ecological_interpretation():
    c = load_contract()
    orientation = c["measurement_uncertainty"]["orientation"]
    assert orientation["stable_subset_required"] is True
    assert "before joining environmental outcomes" in orientation["stable_subset_thresholds"]
    assert c["visible_variation"]["hierarchical_decomposition_required"] is True


def test_phylogeny_and_claim_language_are_bounded():
    c = load_contract()
    assert c["phylogenetic_sensitivity"]["candidate_promotion_gate"] is False
    assert c["claim_language"]["disallowed_primary_label"] == "adaptive-pattern candidate"
    assert c["claim_language"]["allowed_candidate_label"] == "candidate_association_for_functional_validation"


def test_comment_integration_ledger_exists():
    assert (ROOT / "analysis" / "v3" / "supervisor_comment_integration_20260907.md").is_file()

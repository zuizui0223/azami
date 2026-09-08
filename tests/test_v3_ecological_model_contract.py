import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "v3" / "ecological_model_contract.json"


def load_contract():
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_model_is_fixed_before_v3_trait_outcomes():
    c = load_contract()
    assert c["status"] == "fixed_before_original_stream_trait_outcomes_and_trait_environment_coefficients"
    assert c["separation"]["v3_original_stream_trait_values_inspected_when_fixed"] is False
    assert c["separation"]["v3_trait_environment_coefficients_inspected_when_fixed"] is False
    assert c["separation"]["v2_candidate_survival_used_to_choose_model"] is False


def test_hypothesis_environment_mapping_is_bounded():
    mapping = load_contract()["environment_mapping"]
    assert mapping["orientation_H1_primary"]["conditional_predictors"] == ["pr_month"]
    assert mapping["orientation_H1_primary"]["wind_secondary"] == ["pr_month", "sfcWind_month"]
    assert mapping["visible_colour_H2_primary"]["drying_formulation"] == ["pr_month", "rsds_month", "vpd_month"]
    assert mapping["visible_colour_H2_primary"]["thermal_formulation"] == ["pr_month", "rsds_month", "tasmax_month"]
    assert mapping["gross_shape_H3_primary"]["drying_formulation"] == ["pr_month", "rsds_month", "vpd_month"]
    assert mapping["fine_architecture_and_surface"]["primary_hypothesis_test"] is False


def test_hierarchy_dominance_and_spatial_controls_are_primary():
    c = load_contract()
    assert c["analysis_unit"]["nesting"] == ["head", "photo", "observation", "taxon"]
    assert "random-effects" in c["within_taxon_estimand"]["primary_estimand"]
    assert c["spatial_primary_specification"]["included_in_every_primary_within_and_among_fit"] is True
    assert len(c["spatial_primary_specification"]["terms"]) == 8
    assert "Moran" in c["spatial_primary_specification"]["residual_diagnostic"]


def test_bbox_stability_and_colour_context_are_upstream():
    m = load_contract()["measurement_formulation"]
    assert len(m["bbox_sensitive_endpoints"]) == 5
    assert "all four" in m["bbox_primary_stability_rule"]
    assert "five predeclared measurement formulations" in m["bbox_uncertainty_reporting"]
    assert "negative-control" in m["colour_context_negative_control"]


def test_joint_colour_and_shape_not_endpoint_fishing():
    units = load_contract()["response_inferential_units"]
    assert units["visible_colour"]["joint_first"] is True
    assert set(units["visible_colour"]["units"]["circular_hue"]["endpoints"]) == {"corolla_hue_sin", "corolla_hue_cos"}
    assert len(units["visible_colour"]["units"]["closed_colour_composition"]["endpoints"]) == 4
    assert units["gross_shape"]["joint_first"] is True
    mult = load_contract()["multiplicity_and_decision"]
    assert mult["no_global_endpoint_by_raster_screen"] is True
    assert len(mult["primary_hypothesis_families"]) == 3
    assert mult["v2_survival"] == "not a criterion"

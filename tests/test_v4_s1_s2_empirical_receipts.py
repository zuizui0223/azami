"""Regression guards for actual (not synthetic) S1 and S2 ecological receipts.

Only aggregate published values, no original observation/photographer/coords.
"""
from pathlib import Path
import json
import pytest

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"analysis/v4/results/20261021_unified_abcde"


def load(name):
    return json.loads((DIR/name).read_text(encoding="utf-8"))


def test_original_source_is_fully_partitioned_and_not_introduced_complement():
    a=load("S2_INTRODUCED_SOURCE_AUDIT_20261021.json")
    p=a["original_population"]
    assert p["observations"]==665115
    assert (p["native_wild"],p["introduced_wild_explicit_WCVP"],
            p["other_unresolved_captive_unmapped_or_unusable"])==(
                313350,156750,195015)
    assert sum((p["native_wild"],p["introduced_wild_explicit_WCVP"],
                p["other_unresolved_captive_unmapped_or_unusable"]))==665115
    assert p["introduced_orientation_eligible"]==91106
    assert a["source"]["source_json_sha256"]==(
        "bbfc9a459a1a5c230c31a405671633e3032491f92f0123634e5cbb8f4db33562")


def test_observer_sensitivity_keeps_exact_full_species_rain_signal_with_no_same_head_claim():
    a=load("S1_OBSERVER_RAIN_SENSITIVITY_20261021.json")
    o=a["primary_between_rainfall_degrees_per100mm"]
    assert o["original"]["observations"]==170066
    assert o["original"]["taxa"]==244
    assert o["original"]["beta"]==pytest.approx(.4381991922885207)
    assert o["exclude_7_WCVP_collisions"]["observations"]==169612
    assert o["one_original_observation_per_taxon_cell_observer"]["observations"]==91597
    assert o["remove_top_one_percent_source_observers"]["observers_removed"]==929
    for r in o.values():
        assert r["ci"][0]>0 and r["ci"][1]>0
        assert r["within_ci"][0]<0<r["within_ci"][1]
    assert a["same_observation_original_photo_multiplicity"][
        "confirmed_same_head_repeated_measurements"]==0
    assert a["source"]["member_sha256"]==(
        "6d060a09a35e83d4281609942674a68bc005522f5ba147dd1d2abac918c12b1e")


def test_introduced_cline_uses_source_verified_original_pixels_and_no_evolutionary_equality():
    source=load("S2_INTRODUCED_SOURCE_AUDIT_20261021.json")
    result=load("S2_NATIVE_INTRODUCED_CLIMATE_TRANSFER_20261021.json")
    a=result["original_pixel_focal_coverage"]
    assert a["target_orientation_observations"]==173737
    assert a["all_five_climate_complete"]==163687
    assert a["all_five_climate_missing"]==10050
    assert result["source"]["native_introduced_original_wcvp_snapshot_sha256"]==(
        source["source"]["source_json_sha256"])
    assert result["source"]["source_climate_metadata_sha256"]==(
        "905dd41d987a742761d9bb5c8de23e908b730ea70691491483840b998e1ca342")
    assert result["method"]["fixed_predictors"]==[
        "chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15",
        "chelsa_rsds_mean"]
    assert result["method"]["null_equivalence_test"] is False
    for name,info in result["species"].items():
        for state in ("native","introduced"):
            o=info[state]
            assert o["cluster_ci95"][0]<0<o["cluster_ci95"][1]
            assert o["balanced_observer_cell_n"]<o["environment_complete"]
            assert o["cells"]>=70
        assert info["introduced_minus_native"]["approximate_cluster_ci95"][0]<0
        assert info["introduced_minus_native"]["approximate_cluster_ci95"][1]>0
        assert info["predictive_transfer"][
            "native_trained_skill_relative_to_native_constant"]<0
        assert not info["predictive_transfer"]["climate_overlap_gated"]
    assert result["species"]["Cirsium arvense"]["introduced"][
        "environment_complete"]==41995
    assert result["species"]["Cirsium vulgare"]["introduced"][
        "environment_complete"]==38983
    assert result["scientific_boundary"]["not_supports"]


def test_primary_research_register_and_readme_reference_S1_S2_sources():
    register=json.loads((
        ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"
    ).read_text(encoding="utf-8"))
    s1=register["S1"]["observer_source_robustness"]
    s2=register["S2"]["climate_transfer"]
    assert s1["run"]==39537702876
    assert s2["run"]==39539209624
    assert not s2["evolutionary_conclusion_established"]
    assert not register["quality_controls"]["photographer_random_effect_fitted"]
    assert not register["quality_controls"][
        "introduced_native_photo_angle_equal_error_verified"]
    assert "S2_NATIVE_INTRODUCED_CLIMATE_TRANSFER_20261021.json" in s2["receipt"]
    assert (DIR/"S2_NATIVE_INTRODUCED_TRANSFER_RESULT_JA.md").is_file()

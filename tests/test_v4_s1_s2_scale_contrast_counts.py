"""Numerical contract: S1 between-taxon BIO12 and two-species wetter introduced S2.

Tests source-borne observation DENOMINATORS across all native/introduction
sampling levels. A photo-eligible image is NOT an independent new observation
from the same S1 source, nor identical to climate-occurrence sample counts.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"analysis/v4/results/20261022_S1_S2_BETWEEN_WITHIN_NATIVE_INTRO_COUNTS_VERIFIED.json"
REG=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"

def load():
    return json.loads(R.read_text(encoding="utf-8"))


def test_source_original_population_exhausts_exactly_original_665115_records():
    d=load()
    all_=d["original_population_all350"]
    assert all_["source_original_observations"]==665115
    assert all_["source_original_photo_links"]==1122854
    assert all_["S1_native_wild_original_observations"]==313350
    assert all_["S2_explicit_WCVP_introduced_wild_observations"]==156750
    assert all_["unresolved_other_original_observations"]==195015
    assert 313350+156750+195015==665115
    assert all_["S1_native_orientation_eligible"]==170066
    assert all_["S2_explicit_WCVP_introduced_orientation_eligible"]==91106


def test_same_species_native_intro_raw_counts_and_two_distinct_sampling_levels():
    d=load()
    s=d["S2"]
    ar=s["same_species_comparisons"]["Cirsium arvense"]
    vu=s["same_species_comparisons"]["Cirsium vulgare"]
    assert [ar["native"]["original_explicit_WCVP_wild_observations"],
            ar["introduced"]["original_explicit_WCVP_wild_observations"],
            vu["native"]["original_explicit_WCVP_wild_observations"],
            vu["introduced"]["original_explicit_WCVP_wild_observations"]]==[
                78646,70964,62368,85644]
    assert [ar["native"]["original_orientation_eligible_observations"],
            ar["introduced"]["original_orientation_eligible_observations"],
            vu["native"]["original_orientation_eligible_observations"],
            vu["introduced"]["original_orientation_eligible_observations"]]==[
                49072,43195,33637,47833]
    n=s["focal2_totals"]["native"]
    i=s["focal2_totals"]["introduced"]
    assert n["original_explicit_WCVP_wild_observations"]==141014
    assert i["original_explicit_WCVP_wild_observations"]==156608
    assert n["original_orientation_eligible_observations"]==82709
    assert i["original_orientation_eligible_observations"]==91028
    assert n["observer_by_5deg_cell_unique_climate_observations"]==56684
    assert i["observer_by_5deg_cell_unique_climate_observations"]==87494
    assert n["orientation_balanced_observations"]==38116
    assert i["orientation_balanced_observations"]==55944
    assert n["five_CHELSA_cached_original_pixel_complete_oriented_observations"]==82709
    assert i["five_CHELSA_cached_original_pixel_complete_oriented_observations"]==80978
    assert n["five_CHELSA_cached_balanced_model_observations"]==38116
    assert i["five_CHELSA_cached_balanced_model_observations"]==50385
    assert (n["orientation_balanced_observations"]+
            i["orientation_balanced_observations"])==94060
    assert (n["observer_by_5deg_cell_unique_climate_observations"]+
            i["observer_by_5deg_cell_unique_climate_observations"])==144178
    assert s["introduced_wild_share_focal2_of_all_introduced"]>.999
    assert s["introduced_orientation_share_focal2_of_all_introduced"]>.999


def test_between_species_S1_is_positive_but_within_species_gradient_not_confirmed():
    d=load()
    s=d["S1"]
    assert s["n_observations"]==170066 and s["n_accepted_keys"]==244
    assert s["native_244_between_species_BIO12_degrees_per_100mm"]==0.4381991922885207
    assert s["between_95pct_interval"][0]>0
    assert s["within_95pct_interval"][0]<0<s["within_95pct_interval"][1]
    assert s["within_species_annual_precipitation_significant"] is False


def test_both_species_wetter_occurrence_medians_not_consistent_head_orientation():
    d=load()
    for name in ("Cirsium arvense","Cirsium vulgare"):
        r=d["S2"]["same_species_comparisons"][name]
        assert r["introduced"]["median_annual_BIO12_mm"]>r["native"]["median_annual_BIO12_mm"]
        assert r["differences"]["introduced_minus_native_photo_angle_median_degrees"]==-1
        lo,hi=r["differences"]["photo_angle_ge90deg_spatial_cell_bootstrap_95pct_interval_points"]
        assert lo<0<hi
    ar=d["S2"]["same_species_comparisons"]["Cirsium arvense"]
    vu=d["S2"]["same_species_comparisons"]["Cirsium vulgare"]
    assert round(ar["differences"]["introduced_minus_native_median_BIO12_mm"],1)==217.4
    assert round(vu["differences"]["introduced_minus_native_median_BIO12_mm"],1)==293.0
    assert ar["differences"]["introduced_minus_native_photo_angle_ge90deg_percentage_points"]<0
    assert vu["differences"]["introduced_minus_native_photo_angle_ge90deg_percentage_points"]>0
    full=d["S2"]["full_ECDF_descriptive_support"]
    assert full["original_399_resamples_completed_for_both_species"]==399
    assert full["both_introduced_BIO12_distributions_right_shifted_over_broad_source_ranges"] is True
    assert full["both_full_photo_angle_ECDF_difference_simultaneous_source_cell_bands_contain_zero_at_all_display_points"] is True
    assert full["not_a_proof_of_no_phenotypic_response_or_no_genetic_evolution"] is True


def test_primary_narrative_keeps_source_overlap_and_noncausal_boundaries():
    d=load()
    reg=json.loads(REG.read_text(encoding="utf-8"))
    assert d["target_narrative"]["not_reframed_as_geography_adjusted_S1_primary"] is True
    assert d["interpretation_stop_rules"]["S1_and_S2_native_observations_are_overlapping_source_not_two_independent_datasets"] is True
    assert d["interpretation_stop_rules"]["non_significant_photo_angle_movement_is_not_proof_of_no_angle_change"] is True
    assert d["interpretation_stop_rules"]["original_V3_and_S1_frozen_unchanged"] is True
    assert "170066" in reg["contracts"]["epistemic_status"]
    assert "S1_S2_primary_positive_between_species_native_introduced_contrast_20261022" in reg
    assert reg["S1_S2_primary_positive_between_species_native_introduced_contrast_20261022"]["source_S2"]["head_angle_exact_zero_shift_established"] is False
    assert reg["contracts"]["S1_S2_main_scale_contrast_count_ledger"]==str(R.relative_to(ROOT))

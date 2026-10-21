"""Locked empirical 194/191 continuous S1 history: disallow categorical origin promotion."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"analysis/v4/results/20261021_unified_abcde/S1_CONTINUOUS_194_191_52_VERIFIED_20261021.json"

def test_all_194_measured_continuous_and_191_no_taxon_collision():
    r=json.loads(R.read_text(encoding="utf-8"))
    d=r["verified_numerical_summary"]
    clean=d["collision_free_191_taxon_sensitivity"]
    assert r["source_run_id"]==39591120377
    assert r["artifact"]["name"]=="v4-s1-continuous-194-taxa-52-placement-history"
    assert r["artifact"]["sha256"].startswith("sha256:")
    assert d["status"]=="EMPIRICAL_194_S1_CONTINUOUS_IMAGE_ORIENTATION_PHYLOGENY"
    assert (d["n_original_S1_eligible_taxon_keys"],
            d["n_original_images_in_S1"],
            d["n_tree_taxon_keys_with_original_continuous_orientation"],
            d["n_PGLS_historical_placement_scenarios"])==(244,170066,194,52)
    assert clean["n_tree_taxon_keys_with_original_continuous_orientation"]==191
    assert d["source_name_collisions_within_194"]==3
    assert clean["source_name_collisions_within_194"]==0
    assert set(d["source_name_collision_tips"])=={
        "Cirsium_parryi","Cirsium_setidens","Cirsium_yezoense"}
    assert d["S1_original_five_exposures"]==[
       "chelsa_bio01","chelsa_bio04","chelsa_bio12",
       "chelsa_bio15","chelsa_rsds_mean"]
    assert d["measurement_is_continuous_not_binary"] is True
    assert d["evolutionary_nodding_gain_count_claimed"] is False
    assert clean["evolutionary_nodding_gain_count_claimed"] is False
    assert d["reused_data_no_frozen_S1_refit"] is True

def test_rooted_continuous_branch_expected_number_is_not_independent_origin_count():
    r=json.loads(R.read_text(encoding="utf-8"))
    d=r["verified_numerical_summary"]
    clean=d["collision_free_191_taxon_sensitivity"]
    assert d["species_mean_all_edges_increase_ge_10deg_conditional_expected_edge_scenario_range"]==[
        43.518257014727126,48.6842861925713]
    assert clean["species_mean_all_edges_increase_ge_10deg_conditional_expected_edge_scenario_range"]==[
        42.876073028232994,47.9350385995627]
    assert d["species_mean_internal_edges_only_increase_ge_10deg_high_probability_branch_count_range"]==[0.,0.]
    assert clean["species_mean_internal_edges_only_increase_ge_10deg_high_probability_branch_count_range"]==[0.,0.]
    assert d["species_median_internal_edges_only_increase_ge_10deg_high_probability_branch_count_range"]==[0.,0.]
    assert d["n_BM_nonindependent_branch_shifts_cannot_equal_discrete_origins"] is True
    assert r["synopsis"]["independent_nodding_origins_proven"] is False

def test_five_covar_GLS_direction_and_not_specific_bio15():
    d=json.loads(R.read_text(encoding="utf-8"))["verified_numerical_summary"]
    strict=d["collision_free_191_taxon_sensitivity"]
    assert d["diagnostic_same194_GLS_chelsa_bio12_estimate_range"][0]>0
    assert strict["diagnostic_same191_GLS_chelsa_bio12_estimate_range"][0]>0
    assert d["diagnostic_same194_GLS_chelsa_bio12_lower95_range"][0]>0
    assert strict["diagnostic_same191_GLS_chelsa_bio12_lower95_range"][0]<0
    for panel,prefix in ((d,"same194"),(strict,"same191")):
        a=panel[f"diagnostic_{prefix}_GLS_chelsa_bio15_estimate_range"]
        assert a[0]<0<a[1]

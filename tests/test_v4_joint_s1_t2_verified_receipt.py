"""Regression gate for source-authenticated S1→T2 244/194/46 joint receipt."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"analysis/v4/results/20261021_unified_abcde/S1_T2_JOINT_244_194_VERIFIED_20261021.json"
REGISTRY=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"

def test_source_verified_joint_cohorts_and_climate_contrasts():
    r=json.loads(RECEIPT.read_text())
    e=r["results"]
    s=e["model_scope"]
    assert r["workflow_run_id"]==39584744194
    assert r["workflow_artifact_digest"].startswith("sha256:")
    assert s["orientation_eligible_accepted_keys"]==244
    assert s["orientation_eligible_observations"]==170066
    assert s["placement_tree_covered_keys"]==194
    assert s["placement_tree_covered_observations"]==169912
    assert s["grafted_placement_scenarios"]==52
    assert r["source_provenance"]["S1_fixed_exposures"]==[
      "chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean"]
    b12=e["BIO12_between"]
    b15=e["BIO15_between"]
    assert b12["S1_effect_image_degrees"]>0
    assert b12["placement_pointwise_positive_interval_count"]==52
    assert b15["S1_interval"][0]<0<b15["S1_interval"][1]
    assert b15["placement_pointwise_positive_interval_count"]==0

def test_no_source_name_collision_or_pseudo_boolean_promotion():
    r=json.loads(RECEIPT.read_text())
    inv=r["results"]["botanical_state_inventory"]
    assert inv["EAzami_resolved_states"]==17
    assert inv["NA_original_author_states"]==63
    assert inv["EAzami_source_concept_state_rows"]==20
    assert inv["EastAsia_exact_S1_244_matches"]==5
    assert inv["EastAsia_exact_S1_194_tree_matches"]==4
    assert inv["NA_exact_S1_244_matches"]==42
    assert inv["NA_exact_S1_194_tree_matches"]==42
    assert inv["s1_join_exclusion_status_counts"]=={
      "blocked_WCVP_source_name_collision":5,
      "exact_S1_and_tree":46,
      "exact_S1_only_tree_uncovered":1,
      "not_in_S1_orientation_cohort":28,
    }
    assert inv["same_photo_derived_state_or_inferred_threshold_used"] is False
    h=r["results"]["T2_PGLS_52_exploratory"]
    assert h["taxa_with_independent_botanical_U_D_and_actual_S1_angle_and_tree"]==46
    assert h["n_erect"]==40 and h["n_nodding"]==6
    assert h["unknown_root_minimum_steps_scenario_range"]==[5,6]
    assert h["Mk_erect_to_nodding_median_scenario_range"]==[5,7]
    assert h["number_of_independent_real_biological_origins_estimated"] is False
    assert h["rain_climate_dependent_transition_rate_estimated"] is False
    assert r["old_S1_and_frozen_V3_not_refitted"] is True

def test_registered_scope_does_not_confuse_image_angle_with_gravity_or_t3():
    reg=json.loads(REGISTRY.read_text())
    x=reg["S1_T2_shared_phylogenetic_evidence_20261021"]
    assert x["S1_image_vertical_orientation_eligible_accepted_keys"]==244
    assert x["exact_join_S1_image_and_52_tree_and_independent_states"]["total"]==46
    assert x["EAzami_v10_BIO15_axis_specificity_confirmed"] is False
    assert x["independent_Cirsium_origin_events_proven"] is False
    assert x["gravitational_angle_inferred_from_S1_image_axis"] is False
    assert x["not_treat_EAzami_as_independent_out_of_sample_replication"] is True
    assert reg["T3"]["cross_tribe_convergence_inferred"] is False

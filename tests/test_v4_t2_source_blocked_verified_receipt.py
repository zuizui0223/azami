"""Never promote the verified 2026-10-21 S1/T2 regional sharp-null to evolutionary causation."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/"analysis/v4/results/20261021_unified_abcde/T2_SOURCE_BLOCKED_FIVE_CLIMATE_AND_MK_VERIFIED_20261021.json"
REG=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"

def test_real_five_climate_outcome_and_scientific_inference_ceiling():
    r=json.loads(RECEIPT.read_text(encoding="utf-8"))
    x=r["summary"]
    assert r["run_id"]==39587904816
    assert r["artifact"]["name"]=="v4-t2-source-blocked-five-climate"
    assert r["artifact"]["digest"].startswith("sha256:")
    assert (x["S1_original_image_angle_native_taxa"],
            x["S1_placement_tree_original_cohort"],
            x["independent_binary_states_matched_to_S1_and_tree"])==(244,194,46)
    assert (x["matched_NA"],x["matched_EA"],x["nodding_NA"],x["nodding_EA"])==(42,4,3,3)
    assert x["NA_EA_source_blocked_null_maps"]==45920
    assert x["NA_only_null_maps"]==11480
    assert x["S1_five_fixed_exposures"]==[
      "chelsa_bio01","chelsa_bio04","chelsa_bio12","chelsa_bio15","chelsa_rsds_mean"]
    assert x["PGLS_scenario_count"]==52
    mk=x["NA42_independent_nuclear_Mk_prevalence_null"]
    assert mk["nuclear_tip_cohort"]==42
    assert mk["NA_EA_full_46_phylogenetic_null_not_computed"] is True
    for root in ("erect","stationary"):
        d=mk["root_sensitivity"][root]
        assert d["status"]=="conditional_nuclear_Mk_null_computed"
        assert d["null_draws_used"]==1500
        assert all(p>0.8 for p in d["two_axis_maxT_Mk_p"])
    assert x["S1_frozen_model_refitted"] is False
    assert x["climate_dependent_transition_rates_estimated"] is False
    assert x["independent_evolutionary_origins_or_dates_estimated"] is False
    assert r["interpretation_boundaries"]["source_blocked_45920_exact_permutations_are_not_phylogenetic_null"] is True

def test_joint_register_does_not_claim_bio15_specificity_or_t3_events():
    reg=json.loads(REG.read_text(encoding="utf-8"))
    joint=reg["S1_T2_shared_phylogenetic_evidence_20261021"]
    new=joint["source_stratified_five_climate_and_nuclear_null"]
    assert new["region_blocked_null_is_not_phylogeny_correlated"] is True
    assert new["markov_nuclear_null_phylogeny_correlated"] is True
    assert new["no_true_independent_evolutionary_origins_estimated"] is True
    assert new["no_historical_climate_dependent_rate_estimated"] is True
    assert joint["EAzami_v10_BIO15_axis_specificity_confirmed"] is False
    assert reg["S1"]["phylogeny_scenarios_are_independent_evolutionary_events"] is False
    assert reg["T3"]["cross_tribe_convergence_inferred"] is False

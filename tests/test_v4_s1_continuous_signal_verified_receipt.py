"""Evidence-integrity gate for *measured* S1 image angle phylogenetic-signal results.

A historical placement does not become an independent evolutionary origin
when an ancestral Brownian conditional branch change is generated.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"analysis/v4/results/20261021_unified_abcde/S1_CONTINUOUS_PHYLOGENETIC_SIGNAL_VERIFIED_20261021.json"


def test_real_phylogenetic_lambda_frozen_observations_and_nonindependence():
    r=json.loads(RESULT.read_text(encoding="utf-8"))
    s=r["verified_summary"]
    assert r["run_id"]==39594257539
    assert r["archive"]["id"]==12343166776
    assert r["archive"]["sha256"].startswith("sha256:")
    assert s["status"]=="EMPIRICAL_S1_IMAGE_ANGLE_PHYLOGENETIC_SIGNAL_SENSITIVITY"
    assert s["source_cohorts"]=={
       "original194":194,
       "collision_free191":191,
       "collision_free_at_least5":171,
       "collision_free_at_least10":143,
    }
    assert r["hypothesis_climate_covariates_fixed"]==[
       "chelsa_bio01","chelsa_bio04","chelsa_bio12",
       "chelsa_bio15","chelsa_rsds_mean"]
    assert r["original_52_GBOTB_VPhyloMaker2_trees_are_not_independent_species_histories"]
    assert r["measured_angle_is_not_gravity_referenced"]
    assert r["no_proven_independent_downward_capitulum_origins"]
    assert r["main_observational_BIO12_effect_244_not_reestimated"]


def test_climate_corrected_s1_phylogenetic_signal_is_not_supported_in_52_scenarios():
    s=json.loads(RESULT.read_text(encoding="utf-8"))["verified_summary"]
    p=s["result_panels"]
    expected={
       "collision_free191|species_mean|intercept_only|nugget0.05":(0.,.377046632468365,3),
       "collision_free191|species_mean|S1_five_climate_adjusted|nugget0.05":(0.,.31901088460861526,4),
       "collision_free_at_least10|species_mean|S1_five_climate_adjusted|nugget0.05":(0.,.2919903302330273,1),
    }
    for name,(lam_min,lam_max,n_high) in expected.items():
        result=p[name]
        assert result["n_nonindependent_placement_scenarios"]==52
        assert result["lambda_ML_range"]==[lam_min,lam_max]
        assert result["lambda_ML_median"]==0
        assert result["n_lambda_ML_at_least_0_1"]==n_high
        assert result["n_deltaAIC_favours_lambda_more_than_2"]==0
        assert result["lambda_values_not_independent_evolutionary_replicates"]


def test_first_scenario_null_p_not_misreported_as_52_independent_tests():
    r=json.loads(RESULT.read_text(encoding="utf-8"))
    b=r["verified_summary"]["first_scenario_bootstrap"]
    expected={
      "collision_free191__intercept_only__species_mean":(.377046632468365,.0625),
      "collision_free191__S1_five_climate_adjusted__species_mean":(.31901088460861526,.1775),
      "collision_free_at_least10__S1_five_climate_adjusted__species_mean":(0.,1.),
    }
    for name,(lam,p) in expected.items():
        row=b[name]
        assert row["requested_draws"]==399
        assert row["observed_lambda_ML"]==lam
        assert row["boundary_calibrated_bootstrap_p"]==p
        assert row["null"]=="parametric_bootstrap_phylogenetically_independent_tips_given_original_five_climates"
        assert "one fixed" in row["scope"]

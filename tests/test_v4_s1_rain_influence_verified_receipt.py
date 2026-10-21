"""Regression gates for original source-hashed 2026-10-21 BIO12 > 0 S1 audit."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/"analysis/v4/results/20261021_unified_abcde/S1_CONTINUOUS_BIO12_BIO15_TAXON_INFLUENCE_VERIFIED_20261021.json"
REG=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"

def test_verified_191_and_143_taxon_annual_rain_is_positive_without_any_single_taxon():
    doc=json.loads(RESULT.read_text(encoding="utf-8"))
    assert doc["status"].startswith("REAL_ORIGINAL_SOURCE_S1")
    assert doc["source_GitHub_Actions_run_id"]==39598389991
    assert doc["source_artifact"]["id"]==12345610456
    assert doc["source_artifact"]["digest"].startswith("sha256:")
    s=doc["source_verified_stdout_numerical_summary"]
    assert s["status"]=="EMPIRICAL_S1_ANNUAL_RAINFALL_CONTINUOUS_ANGLE_SINGLE_TAXON_INFLUENCE_AUDIT"
    for cohort,n in (("collision_free191",191),("collision_free_at_least10",143)):
        r=s["rain_audit"][cohort]
        b12=r["chelsa_bio12"]
        b15=r["chelsa_bio15"]
        assert b12["source_original_5_joint_predictors"]==[
           "chelsa_bio01","chelsa_bio04","chelsa_bio12",
           "chelsa_bio15","chelsa_rsds_mean"]
        assert b12["increment_original_units"]==100
        assert b12["species_leave_one_out_n"]==n
        assert b12["species_leave_one_out_effect_all_positive"] is True
        assert b12["species_leave_one_out_pointwise_HC3_lower95_positive_count"]==n
        assert b12["HC3_conditional_t95_interval"][0]>0
        assert b12["species_pairs_bootstrap_completed"]==1999
        assert b12["species_pairs_bootstrap_percentile95"][0]>0
        assert b15["HC3_conditional_t95_interval"][0]<0<b15["HC3_conditional_t95_interval"][1]
        assert b15["species_pairs_bootstrap_percentile95"][0]<0<b15["species_pairs_bootstrap_percentile95"][1]
        assert b15["species_leave_one_out_pointwise_HC3_lower95_positive_count"]==0
    assert doc["scientific_claim_bounds"]["original_S1_observation_level_between_within_model_not_refitted"] is True
    assert doc["scientific_claim_bounds"]["precipitation_association_not_evolutionary_causation"] is True


def test_scientific_registry_does_not_claim_seasonality_or_evolved_nodding_origins():
    reg=json.loads(REG.read_text(encoding="utf-8"))
    diag=reg["S1"]["source_verified_BIO12_BIO15_taxon_influence_20261021"]
    assert diag["taxon_191_BIO12_all_leave_one_out_HC3_lower95_positive"] is True
    assert diag["taxon_143_BIO12_all_leave_one_out_143_HC3_lower95_positive"] is True
    assert diag["BIO15_independent_effect_confirmed"] is False
    assert diag["within_species_BIO12_signal_significant"] is False
    assert diag["independent_botanical_downward_origins_inferred"] is False
    assert reg["S1"]["independent_rain_functional_effect_established"] is False

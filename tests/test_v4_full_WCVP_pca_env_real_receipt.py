"""Lock source-verified complete WCVP availability-background PCA-env evidence.

All observations and original TIFFs were private in Actions. The source
receipt records aggregate species-region quantities and SHA-verified complete
CHELSA original backgrounds; do not infer causality or head morphology.
"""
from pathlib import Path
import json
import pytest

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"analysis/v4/results/20261021_unified_abcde"
def data():
    return json.loads((R/"S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json").read_text(encoding="utf-8"))


def test_full_background_85k_WCVP_and_original_SHA_are_source_verified():
    r=data()
    assert r["source"]["full_background_actions_run"]==39557273395
    assert r["source"]["artifact_id"]==12326676034
    assert r["source"]["source_json_sha256"]==(
        "3809d84d478b1178737d5518d045bedda5917e755a4725a7db7749cdb50a3776"
    )
    assert r["source"]["all_five_CHELSA_1981_2010_original_full_TIFFs_SHA256_checked_against_frozen_source_metadata"]
    assert r["source"]["original_photo_angle_values_not_used_in_environment_niche_selection"]
    assert r["global_background"]["TDWG_0p5_degree_mesh_points"]==85414
    assert r["global_background"]["CHELSA_five_predictors_finite_mesh_points"]==85414
    assert r["global_background"]["WCVP_available_native_area_coverage_both_species"]==1
    assert r["global_background"]["WCVP_available_introduced_area_coverage_both_species"]==1
    assert not r["global_background"]["ecospat_exact_equivalence_similarity_test_executed"]
    assert not r["global_background"]["MaxEnt_executed"]


def test_native_to_introduced_BIO12_increases_in_both_species_but_BIO15_differs():
    j=data()["species"]
    ar=j["Cirsium arvense"]["climate_median_native_to_introduced"]
    vu=j["Cirsium vulgare"]["climate_median_native_to_introduced"]
    assert ar["BIO12_mm"][1]>ar["BIO12_mm"][0]
    assert vu["BIO12_mm"][1]>vu["BIO12_mm"][0]
    assert ar["BIO15"][1]<ar["BIO15"][0]
    assert vu["BIO15"][1]>vu["BIO15"][0]
    assert vu["BIO1_C"][1]-vu["BIO1_C"][0]==pytest.approx(4.2)
    assert sum(j["Cirsium arvense"]["PCA_env_2axis_variance_fraction"])>.75
    assert sum(j["Cirsium vulgare"]["PCA_env_2axis_variance_fraction"])>.80


def test_PCA_env_background_correction_changes_D_and_does_not_prove_expansion():
    j=data()["species"]
    for name,s in j.items():
        assert s["WCVP_area_weighted_background_coverage"]=={"native":1,"introduced":1}
        assert s["PCA_env_background_corrected_Schoener_D"]!=s["PCA_env_uncorrected_occurrence_density_D"]
        assert 0<=s["PCA_env_background_corrected_Schoener_D"]<=1
        assert 0<=s["introduced_density_outside_native_95percent_probability_content"]<=1
    assert j["Cirsium arvense"]["PCA_env_background_corrected_Schoener_D"]==pytest.approx(.47393796854546466)
    assert j["Cirsium vulgare"]["PCA_env_background_corrected_Schoener_D"]==pytest.approx(.4424543915282124)
    assert j["Cirsium vulgare"]["introduced_density_outside_native_95percent_probability_content"]==pytest.approx(.43979374526144926)
    assert data()["global_background"]["ecospat_exact_equivalence_similarity_test_executed"] is False


def test_vulgare_Africa_and_South_America_have_more_rain_and_seasonality():
    r=data()["species"]["Cirsium vulgare"]["introduced_regions"]
    assert r["Africa"]["BIO12_median_minus_native_mm"]>200
    assert r["Africa"]["BIO15_median_minus_native"]>20
    assert r["South America"]["BIO12_median_minus_native_mm"]>500
    assert r["South America"]["BIO15_median_minus_native"]>5
    assert r["Africa"]["available_background"]==1
    assert r["South America"]["available_background"]==1
    assert r["Africa"]["D"]<r["South America"]["D"]<r["North America"]["D"]
    assert data()["source"]["private_original_observers_and_locations_not_exported"]


def test_prior_native_occupied_316_cache_D_not_swapped_into_full_global_niche_claim():
    r=data()["prior_partial_cache_comparison"]
    assert r["limited_to_316_original_native_occupied_storage_tiles"]
    for name in ("Cirsium arvense","Cirsium vulgare"):
        assert r["partial_cached_D"][name]<r["full_source_D"][name]
    assert r["refrain_interpreting_cached_D_as_true_global_niche_difference"]
    doc=(R/"S2_FULL_WCVP_TWO_SPECIES_CLIMATE_NICHE_RESULT_JA.md").read_text(encoding="utf-8")
    assert "85,414" in doc
    assert "背景" in doc
    assert "アフリカ" in doc and "南米" in doc

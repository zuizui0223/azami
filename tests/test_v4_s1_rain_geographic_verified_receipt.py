"""Immutable source-backed S1 geographic/rainfall evidence gate.

A robust *non-spatial* BIO12 species coefficient is not a statistically
isolated global rain effect once native taxon spatial distribution is modeled.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
R=ROOT/"analysis/v4/results/20261022_S1_RAIN_SPATIAL_GEOGRAPHY_SOURCE_VERIFIED.json"
REG=ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json"

def test_exact_original_photo_and_five_unchanged_climates():
    r=json.loads(R.read_text(encoding="utf-8"))
    assert r["source_GitHub_run_id"]==39659502828
    assert r["source_artifact"]["id"]==12373246235
    assert r["source_artifact"]["digest"].startswith("sha256:")
    assert r["observation_source_audit"]["n_S1_five_climate_eligible"]==170066
    assert r["observation_source_audit"]["n_S1_accepted_keys"]==244
    assert r["observation_source_audit"]["n_source_invalid_coordinates_excluded"]==0
    assert r["observation_source_audit"]["five_selected_exposures"]==[
       "chelsa_bio01","chelsa_bio04","chelsa_bio12",
       "chelsa_bio15","chelsa_rsds_mean"]
    assert r["original_native_S1_decoded_sha256"]=="342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
    assert r["inference_boundaries"]["original_S1_between_within_244_model_not_refitted"] is True
    assert r["inference_boundaries"]["source_equal_5deg_cells_are_not_equal_area_geodesic_cells"] is True
    assert r["inference_boundaries"]["photos_continue_as_image_vertical_not_gravity"] is True

def test_photo_and_5deg_cells_geographic_adjustments_no_longer_support_global_BIO12():
    r=json.loads(R.read_text(encoding="utf-8"))
    panels=r["models_compact"]
    assert len(panels)==6
    for mode in ("photo_weighted","equal_5deg_spatial_cell"):
        for cohort,n in (("clean_native",237),("clean_tree",191),("clean_tree_min10",143)):
            row=panels[f"{mode}|{cohort}"]
            assert row["n_taxa"]==n
            assert row["no_geo"]["CI"][0]>0
            assert row["spherical_second_order"]["CI"][0]<0<row["spherical_second_order"]["CI"][1]
            assert row["longitude_composition"]["CI"][1]>0
    strict=panels["photo_weighted|clean_tree"]
    assert abs(strict["no_geo"]["coef"]-.5075715791859658)<1e-10
    assert abs(strict["spherical_second_order"]["coef"]-.044477493217758735)<1e-10
    assert strict["longitude_composition"]["CI"][0]<0
    assert panels["equal_5deg_spatial_cell|clean_tree"]["longitude_composition"]["CI"][0]>0

def test_rain_angle_east_longitude_composition_and_within_band_uncertainty():
    r=json.loads(R.read_text(encoding="utf-8"))["models_compact"]
    full=r["photo_weighted|clean_tree"]
    west=full["within_sector"]["W180_W030"]
    east=full["within_sector"]["E060_E180"]
    assert (east["n_taxa"],west["n_taxa"])==(40,93)
    assert east["mean_native_source_BIO12_mm_per_year"]>west["mean_native_source_BIO12_mm_per_year"]
    assert east["mean_original_image_angle_deg"]>west["mean_original_image_angle_deg"]
    assert full["leave_one_longitude_sector_out"]["E060_E180"]["CI"][0]<0
    assert full["leave_one_longitude_sector_out"]["E060_E180"]["CI"][1]>0
    assert east["within_sector_5CHELSA_HC3"]["CI"][0]<0
    # A single nominal positive east-band result after 5-degree cell balance
    # cannot be called universal, independent or multiple-testing-confirmed.
    eastcell=r["equal_5deg_spatial_cell|clean_tree"]["within_sector"]["E060_E180"]
    assert eastcell["within_sector_5CHELSA_HC3"]["CI"][0]>0

def test_scientific_registry_corrected_claim_does_not_forge_causal_rain_evolution():
    reg=json.loads(REG.read_text(encoding="utf-8"))
    diag=reg["S1"]["source_verified_native_geographical_confounds_20261022"]
    assert diag["run_id"]==39659502828
    assert diag["all_six_spherical_adjusted_cohort_HC3_BIO12_interval_includes_zero"] is True
    assert diag["independent_annual_precipitation_effect_after_spatial_controls_confirmed"] is False
    assert diag["original_S1_observation_level_within_between_fit_reestimated"] if "original_S1_observation_level_within_between_fit_reestimated" in diag else diag["source_S1_observation_within_between_model_refitted"] is False
    assert diag["location_coordinates_shared_in_output"] is False
    assert reg["S1"]["independent_rain_functional_effect_established"] is False
    assert reg["T3"]["cross_tribe_convergence_inferred"] is False

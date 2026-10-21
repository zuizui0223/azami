"""Source-hashed S1 annual-rain primary and S2 within-introduction photo-angle evidence.

The source-numerical 173737 archived original image-angle observations are
NOT genetic, head-identity or physical-gravity calibrated measurements.
The BIO15 between-species S1 tendency has a zero-crossing 95% interval.
"""
from pathlib import Path
import json
import pytest

ROOT=Path(__file__).resolve().parents[1]
A=ROOT/"analysis/v4/results/20261021_unified_abcde"
def load():
    return json.loads((A/"S2_BIO12_PRIMARY_BIO15_SECONDARY_VERIFIED_20261021.json").read_text())


def test_annual_precip_S1_positive_and_seasonality_only_tendency():
    result=load()
    s1=result["S1_cross_taxon_reference"]
    assert s1["n"]==170066 and s1["taxa"]==244
    assert s1["main_BIO12_between_deg_per_plus100mm"]==pytest.approx(.4381991922885207)
    assert s1["main_BIO12_between_95pct"][0]>0
    assert s1["main_BIO12_Holm_p"]<.01
    assert s1["secondary_BIO15_between_deg_per_plus10"]>0
    low,hi=s1["secondary_BIO15_between_95pct"]
    assert low<0<hi
    wlow,whi=s1["secondary_BIO15_within_95pct"]
    assert wlow<0<whi
    assert not s1["seasonality_independently_supported"]
    assert s1["S1_and_S2_estimands_differ"]


def test_full_verified_original_rain_and_seasonality_photo_source_and_observer_counts():
    j=load()
    s=j["source"]
    assert s["github_action_run"]==39559452415
    assert s["artifact_id"]==12326369656
    assert s["member_sha256"]==(
        "084b28ec90c20bf08826fca412e4b4f685257334affc3dcd14f20899f66ee0f7"
    )
    assert s["focal_original_image_orientation_eligible"]==173737
    assert s["focal_original_observer_5degree_cell_balanced"]==94060
    assert s["new_original_raster_1981_2010_checksum_verified"]


def test_North_America_frozen_two_species_primary_BIO12_wet_bins_do_not_prove_nodding():
    a=load()["primary_BIO12_North_America"]
    for species in ("Cirsium arvense","Cirsium vulgare"):
        r=a[species]
        assert r["above_native_q95"]["cells"]>=15
        assert r["above_native_q95"]["observations"]>=500
        q=r["cell_bootstrap_q025_q50_q975"]
        assert q[0]<0<q[2]
        local=r["within_intro_top_third_minus_bottom_third"]
        assert local["high_cells"]>=15 and local["low_cells"]>=15
        ci=local["cell_bootstrap_q025_q50_q975"]
        assert ci[0]<0<ci[2]
    assert a["all_reported_interval_bounds_include_zero"]
    assert not a["direction_consistently_positive"]
    assert a["no_supported_rain_shelter_photo_angle_cline"]


def test_small_geographic_cell_support_and_BIO15_uncertainty_are_not_dropped():
    d=load()
    reg=d["other_introduced_BIO12_regions"]["Cirsium vulgare"]
    assert reg["Africa"]["uncertainty_interval_not_identifiable"]
    assert reg["South America"]["uncertainty_interval_not_identifiable"]
    assert reg["Australasia"]["local_high_low_cell_bootstrap_q025_q50_q975"][0]<0
    s=d["secondary_BIO15_introduced_results"]
    assert s["Cirsium arvense_North_America"]["high_BIO15_local_tercile_minus_low_tercile_photo_angle_deg"]<0
    assert s["Cirsium vulgare_North_America"]["cell_bootstrap_q025_q50_q975"][0]<0
    assert s["Cirsium vulgare_Africa"]["bootstrap_interval_not_identifiable"]
    assert s["Cirsium vulgare_South_America"]["bootstrap_interval_not_identifiable"]
    assert "Exploratory" in s["secondary_status"]
    assert any("1981–2010" in line for line in d["technical_and_causal_limits"])


def test_registry_and_Main_result_doc_are_consistent():
    registry=json.loads((ROOT/"analysis/v4/RAIN_ORIENTATION_RESEARCH_DECISIONS_20261021.json").read_text())
    s1=registry["S1"]["precipitation_priority"]
    assert s1["primary"]=="chelsa_bio12"
    assert s1["secondary"]=="chelsa_bio15"
    assert "not an independently resolved" in s1["secondary_status"] or (
       "do not call" in s1["secondary_status"])
    assert registry["S2"]["climate_transfer"]["primary_rain_axis"]=="chelsa_bio12"
    doc=(A/"S1_S2_BIO12_MAIN_BIO15_SECONDARY_RESULTS_JA.md").read_text()
    assert "降水季節性" in doc
    assert "年降水BIO12" in doc
    assert "未検証" in doc or "できない" in doc

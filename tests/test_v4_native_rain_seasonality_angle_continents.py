"""Synthetic source tests for rainfall and seasonal niche vs photo head angle.

The real geographic photo-orientation and CHELSA estimates must come from
private source-checksummed GitHub Actions, never fabricated from these toys.
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import native_rain_seasonality_angle_continents as m


def frames():
    rng=np.random.default_rng(20261021)
    records=[]
    for species in ("Cirsium arvense","Cirsium vulgare"):
        for state,region in [
            ("native",1),("introduced",7),("introduced",2),
            ("introduced",8),("introduced",5)
        ]:
            for cell in range(90):
                lat=-40+(cell//25)*16
                lon=-175+(cell%25)*14
                # The introduced regions vary annually wet or seasonal climate,
                # without assigning groups using any photo-angle values.
                for j in range(3):
                    rain=530+cell*3 + rng.normal(scale=15)
                    season=40+(cell%35)*.9+rng.normal(scale=.8)
                    if state=="introduced":
                        rain+=45
                        if region in (2,8):
                            season+=8
                    angle=55+.025*rain+.25*season+rng.normal(scale=2)
                    records.append({
                        "obs_id":f"{species}-{state}-{region}-{cell}-{j}",
                        "user_id":f"simulated-human-{species}-{state}-{region}-{cell}-{j}",
                        "taxon_name":species,"range_state":state,
                        "tdwg1":region,"latitude":lat,"longitude":lon,
                        "observed_on":f"202{j+1}-{(cell%12)+1:02d}-15",
                        "quality_grade":"research",
                        "orientation_angle_degrees_median":float(angle),
                        "chelsa_bio12":float(rain),
                        "chelsa_bio15":float(season)
                    })
    return pd.DataFrame(records)


def test_native_q95_and_continent_terciles_are_direct_same_photo_measurements(monkeypatch):
    frame=frames()
    monkeypatch.setattr(m,"SOURCE_ORIGINAL_ALL_HEADS",len(frame))
    r=m.run_angle_seasonality(frame,{"chelsa_bio12":{"source_SHA256":"a"*64},
                                       "chelsa_bio15":{"source_SHA256":"b"*64}})
    assert r["status"]=="EXECUTED_ORIGINAL_TWO_CHELSA_RAIN_VARIABLES_VS_ORIGINAL_IMAGE_HEAD_ANGLES"
    assert len(r["distinct_focal_species"])==2
    assert r["source_original_eligible_image_head_angles"]==len(frame)
    assert not r["private_observation_user_geo_or_image_data_published"]
    for name,s in r["distinct_focal_species"].items():
        assert set(s["introduced_continents"])=={
            "North America","Africa","South America","Australasia"}
        for region,analysis in s["introduced_continents"].items():
            assert set(analysis["axes"])==set(m.AXES)
            for variable,case in analysis["axes"].items():
                ter=case["upper_vs_lower_introduced_climate_tercile"]
                assert ter["source_cells"]>=15
                assert ter["climate_terciles_are_outcome_blind"]
                assert ter["upper_climate_third"]["source_photo_orientation_measurements"]>=30
                assert ter["lower_climate_third"]["source_photo_orientation_measurements"]>=30
                assert ter["photo_angle_tercile_contrast"]["not_gravity_angle_or_within_genet_change"]
                q=case["above_native_q95_vs_other_within_same_introduced_continent"]
                assert q["angle_reference"]=="image vertical, not gravity"
    encoded=json.dumps(r)
    assert "simulated-human" not in encoded
    assert '"obs_id"' not in encoded


def test_climate_axis_selection_does_not_depend_on_source_image_angles(monkeypatch):
    a=frames()
    monkeypatch.setattr(m,"SOURCE_ORIGINAL_ALL_HEADS",len(a))
    before=m.run_angle_seasonality(a,{})
    b=a.copy()
    b["orientation_angle_degrees_median"]=180-b["orientation_angle_degrees_median"]
    after=m.run_angle_seasonality(b,{})
    for species in before["distinct_focal_species"]:
        for continent in before["distinct_focal_species"][species]["introduced_continents"]:
            for climate in m.AXES:
                x=before["distinct_focal_species"][species]["introduced_continents"][continent]["axes"][climate]
                y=after["distinct_focal_species"][species]["introduced_continents"][continent]["axes"][climate]
                t=x["upper_vs_lower_introduced_climate_tercile"]
                z=y["upper_vs_lower_introduced_climate_tercile"]
                assert t["introduced_climate_tercile_lower_threshold"]==pytest.approx(
                    z["introduced_climate_tercile_lower_threshold"])
                assert t["introduced_climate_tercile_upper_threshold"]==pytest.approx(
                    z["introduced_climate_tercile_upper_threshold"])
                assert t["upper_climate_third"]["source_photo_orientation_measurements"]==(
                    z["upper_climate_third"]["source_photo_orientation_measurements"])


def test_source_absence_cannot_count_as_zero_head_photo_response(monkeypatch):
    a=frames()
    monkeypatch.setattr(m,"SOURCE_ORIGINAL_ALL_HEADS",len(a))
    a.loc[a.index[0],"chelsa_bio15"]=np.nan
    with pytest.raises(ValueError,match="Original BIO12 and BIO15 pixels missing"):
        m.run_angle_seasonality(a,{})
    a=frames()
    a.loc[a.index[0],"obs_id"]=a.loc[a.index[1],"obs_id"]
    with pytest.raises(ValueError,match="photographer validation"):
        m.run_angle_seasonality(a,{})


def test_tiny_region_and_no_native_reference_does_not_emit_fake_difference():
    x=frames()
    native=x.loc[x.range_state.eq("native")&x.taxon_name.eq("Cirsium vulgare")]
    tiny=x.loc[x.range_state.eq("introduced")&x.tdwg1.eq(2)&
               x.taxon_name.eq("Cirsium vulgare")].head(8).copy()
    from analysis.v4.native_intro_continent_orientation import hypothesis_within_continent
    from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
    native=dated_and_balanced(native)
    tiny=dated_and_balanced(tiny)
    case=hypothesis_within_continent(native,tiny,
          climate_variable="chelsa_bio15",scenario="strict-seasonality")
    assert case["status"].startswith("STOP_")
    assert case["within_continent_angle_change"] is None

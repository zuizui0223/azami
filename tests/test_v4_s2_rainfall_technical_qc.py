"""Focused S2 wet-band observation-quality/season/photographer QA (synthetic)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.v4.native_spatial import cells
from analysis.v4.s2_rainfall_occupancy_trait import (
    _source_rain_bands, rainfall_wet_band_technical_audit
)


def source_and_photo():
    rng=np.random.default_rng(20261021)
    rows=[]
    for state in ("native","introduced"):
        for i in range(320):
            rain=500.+(i%100)*10 if state=="native" else 710.+(i%100)*10
            lat=-60+(i%15)*7.
            lon=-168+((i//15)%20)*14.
            rows.append({
                "obs_id":f"{state}-{i}",
                "taxon_name":"Cirsium vulgare",
                "range_state":state,
                "user_id":f"{state}-u{i}",
                "_cell":int(cells([lat],[lon])[0]),
                "tdwg1":1 if state=="native" else 7,
                "observed_on":f"2024-0{6+(i%2)}-15",
                "quality_grade":"research",
                "source_photo_count":1+(i%3==0),
                "chelsa_bio12":rain,
                "chelsa_bio01":8+rng.normal(),
                "chelsa_bio04":680+rng.normal(0,40),
                "chelsa_bio15":35+rng.normal(),
                "chelsa_rsds_mean":12+rng.normal(),
                "orientation_angle_degrees_median":
                    float(np.clip(60+rng.normal(0,12)+(state=="introduced")*3,0,180)),
            })
    photo=pd.DataFrame(rows)
    source=photo.drop(columns=["orientation_angle_degrees_median"])
    return source,photo


def test_month_quality_reporter_qc_and_no_pose_illusion():
    source,photo=source_and_photo()
    native=source.loc[source.range_state.eq("native")]
    q=_source_rain_bands(native)
    ids=set(photo.obs_id)
    obs=rainfall_wet_band_technical_audit(
        source,photo,ids,q,seed=12345)
    assert obs["status"]=="ORIGINAL_SOURCE_WET_BAND_TECHNICAL_BIAS_AUDIT"
    assert not obs["verified_camera_roll_pitch_or_actual_gravity"]
    assert not obs["verified_flowering_stage"]
    for field in ("native_source_metadata","introduced_source_metadata"):
        m=obs[field]
        assert m["same_selected_source_observation_has_eligible_archived_angle"]==1.0
        assert m["original_multiple_linked_photographs_fraction_NOT_same_head"]>0
        assert m["original_quality_grade_research_fraction_NOT_photo_pose"]==1.0
    std=obs["same_wet_band_month_x_record_grade_standardised_angle"]
    assert std["status"]=="EXPLORATORY_COMMON_MONTH_GRADE_STANDARDISATION"
    assert std["supported_month_grade_strata"]>=1
    assert len(std["introduced_minus_native_Q10_Q50_Q90"])==3
    assert std["standardised_native_cells"]>=15
    assert std["standardised_introduced_cells"]>=15
    assert std["spatial_interval_status"] in (
        "EXPLORATORY_FIXED_SEASON_GRADE_CELL_BOOTSTRAP",
        "INSUFFICIENT_VALID_MONTH_GRADE_SPATIAL_RESAMPLES",
    )
    trim=obs["same_wet_band_excluding_top_1pct_source_observers"]
    assert trim["status"]=="EXPLORATORY_REMOVE_TOP_1PCT_ORIGINAL_REPORTERS"
    assert "obs_id" not in str(obs)
    assert "user_id" not in str(obs)


def test_original_source_angle_eligibility_changes_but_climate_reference_not():
    source,photo=source_and_photo()
    q=_source_rain_bands(source.loc[source.range_state.eq("native")])
    altered=photo.loc[
        ~((photo.range_state=="introduced")&(photo.obs_id.str.endswith("1")))
    ].copy()
    a=rainfall_wet_band_technical_audit(
        source,photo,set(photo.obs_id),q,seed=77)
    b=rainfall_wet_band_technical_audit(
        source,altered,set(altered.obs_id),q,seed=77)
    assert a["native_source_metadata"]==b["native_source_metadata"]
    assert (b["introduced_source_metadata"]
             ["same_selected_source_observation_has_eligible_archived_angle"]
            <a["introduced_source_metadata"]
              ["same_selected_source_observation_has_eligible_archived_angle"])
    assert not b["verified_camera_roll_pitch_or_actual_gravity"]


def test_does_not_fabricate_camera_pose_or_obs_metadata():
    source,photo=source_and_photo()
    q=_source_rain_bands(source.loc[source.range_state.eq("native")])
    invalid=rainfall_wet_band_technical_audit(
        source.drop(columns=["source_photo_count"]),photo,
        set(photo.obs_id),q,seed=0)
    assert invalid["status"]=="MISSING_ORIGINAL_GRADE_OR_PHOTO_LINK_METADATA"
    assert invalid["no_camera_pose_inferred"]

def test_zero_balancing_weights_are_exclusions_not_negative_weights():
    import numpy as np
    import pytest
    from analysis.v4.native_climate_trait_coupling import weighted_quantiles
    a=np.array([-999.,10.,20.,1000.])
    w=np.array([0.,1.,1.,0.])
    q=weighted_quantiles(a,w,(.1,.5,.9))
    assert all(10.<=v<=20. for v in q)
    assert q[0]<=q[1]<=q[2]
    with pytest.raises(ValueError,match="nonnegative"):
        weighted_quantiles(a,np.array([-1.,1.,1.,1.]))
    with pytest.raises(ValueError,match="positive total"):
        weighted_quantiles(a,np.zeros(4))

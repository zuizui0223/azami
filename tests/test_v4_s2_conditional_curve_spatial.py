"""Same BIO12 versus full head-angle CDF; cell-bootstrap and provenance tests."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4 import s2_conditional_curve_spatial as s2
from analysis.v4 import s2_joint_rain_angle_landscape as original


def example(monkeypatch):
    rng=np.random.default_rng(20261022)
    rows=[]
    for k,name in enumerate(s2.SPECIES):
        for q,state in enumerate(("native","introduced")):
            for cell in range(40):
                for rain_level in range(15):
                    for rep in range(2):
                        # Overlapping native/introduced rain climates are
                        # necessary to inspect a conditional phenotype curve.
                        rain=300+70*rain_level+q*30+rng.normal(0,45)
                        angle=np.clip(54+0.017*rain+3*q+
                            rng.normal(0,22),0,180)
                        j=len(rows)
                        rows.append({
                            "obs_id":f"obs{j}",
                            "user_id":f"user{j}",
                            "taxon_name":name,"range_state":state,
                            "_cell":1000*k+500*q+cell,
                            "latitude":float(-50+cell*2.5),
                            "longitude":float(-160+cell*7),
                            "tdwg1":1 if state=="native" else 7,
                            "observed_on":"2024-07-20",
                            "chelsa_bio01":11+rng.normal(),
                            "chelsa_bio04":770+rng.normal(0,55),
                            "chelsa_bio12":rain,
                            "chelsa_bio15":34+rng.normal(0,4),
                            "chelsa_rsds_mean":14+rng.normal(0,.7),
                            s2.ANGLE:float(angle)
                        })
    photo=pd.DataFrame(rows)
    wild=photo.drop(columns=[s2.ANGLE])
    monkeypatch.setattr(s2,"EXPECTED_PHOTO",len(photo))
    monkeypatch.setattr(s2,"EXPECTED_BALANCED",len(photo))
    monkeypatch.setattr(original,"EXPECTED_PHOTO",len(photo))
    monkeypatch.setattr(original,"EXPECTED_SOURCE",len(wild))
    # Test sample is tiny compared with 94,060 real observations; keep the
    # plotting grid coarser rather than requiring fake original support.
    monkeypatch.setattr(original,"RAIN_DISPLAY_KNOTS",21)
    receipt=original.evaluate(wild,photo)
    return photo,receipt


def test_source_photo_joint_full_conditional_ecdf_geographic_sensitivity(monkeypatch):
    photo,receipt=example(monkeypatch)
    result=s2.evaluate(photo,receipt,draws=90)
    assert result["status"]=="SOURCE_VERIFIED_SAME_RAIN_UNSELECTED_ANGLE_FULL_CDF_CELL_SENSITIVITY"
    assert result["no_chosen_width_quantile_extreme_tail_effect_metric"]
    assert not any(result["source_privacy"].values())
    for name,case in result["species"].items():
        assert case["status"]=="EXPLORATORY_SIMULTANEOUS_SAME_RAIN_FULL_HEAD_CDF_CELL_ENVELOPE"
        assert case["native_vs_ALL_introduced_only"]
        assert case["full_image_angle_axis_not_preselected_threshold"]
        assert case["common_rainfall_bins_originally_supported"]>=5
        diag=case["simultaneous_interval_diagnostic"]
        assert diag["draws_valid_for_all_common_rain_bins"]>=81
        assert diag["absolute_supremum_deviation_95percent"]>=0
        obs=case["observed_conditional_head_angle_CDF_introduced_minus_native"]
        low=case["simultaneous_95pct_spatial_sensitivity_lower"]
        high=case["simultaneous_95pct_spatial_sensitivity_upper"]
        assert len(obs)==len(case["rainfall_display_bins_originally_supported_both_ranges"])
        for accepted,a,b,c in zip(case["rainfall_display_bins_originally_supported_both_ranges"],
                                  obs,low,high):
            if accepted:
                aa=np.asarray(a,float);bb=np.asarray(b,float);cc=np.asarray(c,float)
                assert len(aa)==len(case["image_angle_CDF_upper_edges_degrees"])
                assert np.all(bb<=aa+1e-10)
                assert np.all(cc>=aa-1e-10)
                assert abs(aa[-1])<1e-8
            else:
                assert all(v is None for v in a+b+c)


def test_unverified_photo_joint_refuses_bootstrap_and_privacy(monkeypatch,tmp_path):
    photo,receipt=example(monkeypatch)
    tampered=json.loads(json.dumps(receipt))
    tampered["species"]["Cirsium vulgare"]["native"][
        "original_photo_joint_mass"][0][0]+=0.04
    with pytest.raises(ValueError,match="inconsistent"):
        s2.evaluate(photo,tampered,draws=60)
    altered=photo.drop(columns=["obs_id"])
    with pytest.raises(ValueError,match="photo identity"):
        s2.evaluate(altered,receipt,draws=60)
    saved=s2.save(photo,receipt,tmp_path,draws=60)
    output=(tmp_path/"S2_full_conditional_rain_angle_ECDF_cell_sensitivity.json").read_text()
    assert json.loads(output)["status"]==saved["status"]
    for protected in ('"obs_id"','"user_id"','"latitude"','"longitude"'):
        assert protected not in output

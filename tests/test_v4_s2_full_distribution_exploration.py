"""Distribution-first, metric-free S2 original-source comparison."""
import json
import numpy as np
import pandas as pd
import pytest
from analysis.v4.native_spatial import cells
from analysis.v4.s2_full_distribution_exploration import explore,save


def sample():
    rng=np.random.default_rng(20261021)
    rows=[]
    for species in ("Cirsium arvense","Cirsium vulgare"):
        for status in ("native","introduced"):
            for i in range(300):
                lat=-54+(i%17)*6
                lon=-170+((i//17)%12)*23
                rain=400+(i%100)*11+(status=="introduced")*170
                angle=56+rng.normal(0,14)+(i%9==0)*40
                rows.append({
                    "obs_id":f"{species}:{status}:{i}",
                    "taxon_name":species,"range_state":status,
                    "latitude":float(lat),"longitude":float(lon),
                    "user_id":f"observer:{species}:{status}:{i}",
                    "_cell":int(cells([lat],[lon])[0]),
                    "observed_on":"2024-07-18",
                    "chelsa_bio01":10+rng.normal(),
                    "chelsa_bio04":750+rng.normal(0,90),
                    "chelsa_bio12":float(rain),
                    "chelsa_bio15":30+rng.normal(0,7),
                    "chelsa_rsds_mean":13+rng.normal(),
                    "orientation_angle_degrees_median":float(np.clip(angle,0,180)),
                })
    p=pd.DataFrame(rows)
    return p.drop(columns=["orientation_angle_degrees_median"]),p


def test_no_preselected_width_or_tail_metrics():
    wild,photos=sample()
    o=explore(wild,photos)
    assert o["status"]=="ORIGINAL_NATIVE_ALL_INTRODUCED_FULL_EMPIRICAL_CURVES"
    assert o["no_selected_outcome_threshold_or_effect_score"]
    for sp in o["species"].values():
        assert set(sp)=={"rainfall_environment_full_distribution",
                        "image_angle_full_distribution"}
        for case in sp.values():
            x=case["display_only_unselected_common_knots"]
            a=np.asarray(case["native_empirical_cdf"])
            b=np.asarray(case["all_introduced_empirical_cdf"])
            assert len(x)==len(a)==len(b)
            assert len(x)>100
            assert np.diff(x).min()>0
            assert np.diff(a).min()>=-1e-10
            assert np.diff(b).min()>=-1e-10
            assert a[-1]==b[-1]==1
            # No preselected quantile or threshold: uncertainty envelopes
            # are reported for EVERY original empirical-CDF display knot.
            band=case["whole_curve_spatial_uncertainty"]
            assert band["status"]=="EXPLORATORY_FULL_ECDF_GEOGRAPHIC_CELL_UNCERTAINTY"
            assert band["no_biological_tail_or_width_metric_selected"]
            assert band["resamples_valid"]>=.90*band["resamples_requested"]
            delta=np.asarray(band["introduced_minus_native_ECDF"])
            low=np.asarray(band["simultaneous_95pct_spatial_sensitivity_lower"])
            high=np.asarray(band["simultaneous_95pct_spatial_sensitivity_upper"])
            assert np.allclose(delta,b-a)
            assert len(low)==len(delta)==len(high)==len(x)
            assert np.all(low<=delta+1e-10)
            assert np.all(high>=delta-1e-10)
            assert np.all(low>=-1) and np.all(high<=1)
            assert not {"Q10","Q50","Q90","width","IQR","score"} & set(case)


def test_inverting_image_angles_does_not_change_source_climate():
    wild,photos=sample()
    a=explore(wild,photos)
    changed=photos.copy()
    changed["orientation_angle_degrees_median"]=180-changed["orientation_angle_degrees_median"]
    b=explore(wild,changed)
    name="Cirsium arvense"
    assert a["species"][name]["rainfall_environment_full_distribution"]==(
        b["species"][name]["rainfall_environment_full_distribution"])
    assert a["species"][name]["image_angle_full_distribution"]!=(
        b["species"][name]["image_angle_full_distribution"])


def test_original_identity_required_and_only_aggregate_exported(tmp_path):
    wild,photos=sample()
    with pytest.raises(ValueError,match="source identity"):
        explore(wild.drop(columns=["user_id"]),photos)
    o=save(wild,photos,tmp_path)
    raw=(tmp_path/"S2_unselected_full_rain_and_head_angle_distributions.json").read_text()
    assert json.loads(raw)["status"]==o["status"]
    for term in ('"obs_id"','"user_id"','"latitude"','"longitude"'):
        assert term not in raw
    assert not any(o["privacy"].values())

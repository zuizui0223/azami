"""Synthetic fixed-five CHELSA and native-vs-introduced transfer tests."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from analysis.v4.native_cached_climate import sample_cached_pixels,VARIABLES
from analysis.v4.native_introduced_climate_transfer import (
    dated_and_balanced,standardization,cell_weighted_fit,
    climate_common_support,transfer_against_native_baseline,evaluate_focal
)


def toy_focal(per_state=120):
    rng=np.random.default_rng(20261021)
    records=[]
    for species in ("Cirsium arvense","Cirsium vulgare"):
        for state in ("native","introduced"):
            for i in range(per_state):
                cell=i%30
                latitude=-40.+(cell//10)*14.
                longitude=-106.+(cell%10)*9.
                climate={name:float(rng.normal()) for name in VARIABLES}
                climate["chelsa_bio12"]=float(500+200*rng.normal())
                y=45+climate["chelsa_bio12"]*.035+0.3*climate["chelsa_bio01"]+rng.normal(scale=.2)
                item={
                    "taxon_name":species,"range_state":state,"obs_id":f"{species}-{state}-{i}",
                    "user_id":f"person-{state}-{i%47}",
                    "latitude":latitude,"longitude":longitude,
                    "observed_on":f"{2017+(i%7)}-{(i%12+1):02d}-05",
                    "orientation_angle_degrees_median":float(y),
                }
                item.update(climate)
                records.append(item)
    return pd.DataFrame(records)


def test_true_pixel_value_not_broad_five_degree_average(tmp_path):
    grid={
        "transform":[1.,0.,-180.,0.,-1.,90.],
        "tiles":{"18_36":{"file":"tiles-000.npz","row_off":89,"col_off":180}},
        "nodata":65535
    }
    np.savez_compressed(tmp_path/"tiles-000.npz",**{"18_36":np.array([
        [13,14,15],[16,17,18],[19,20,21]
    ],dtype=np.uint16)})
    # Two longitudes differ by one original raster pixel inside the same 5-degree storage tile.
    a=sample_cached_pixels(tmp_path,grid,np.array([0.1,0.1]),
                          np.array([0.1,1.1]))
    assert np.array_equal(a,np.array([13,14]))
    # Outside cached tile returns missing, never 0 or a coarse-cell mean.
    b=sample_cached_pixels(tmp_path,grid,np.array([20.]),np.array([30.]))
    assert np.isnan(b[0])


def test_one_observation_per_original_observer_x_cell_and_time_index():
    df=toy_focal(120)
    a=dated_and_balanced(df)
    assert len(a)>0
    assert not a.duplicated(["taxon_name","range_state","_cell","user_id"]).any()
    assert a._month.notna().all()
    assert set(a._year)=={2023}


def test_two_focal_native_introduced_bio12_slopes_are_separately_identifiable():
    f=toy_focal(240)
    a=evaluate_focal(f,{"type":"toy-cache-not-a-real-climate-extraction"})
    assert a["assertions"]["evolution_rate_or_plasticity_causally_identified"] is False
    assert len(a["focal_species"])==2
    for label,result in a["focal_species"].items():
        assert result["status"]=="estimated_provisional_S2_native_introduced_climate_contrast"
        assert result["estimated_within_range_slope"]["native"]["source_distinct_cells"]>=15
        assert result["estimated_within_range_slope"]["introduced"]["source_distinct_cells"]>=15
        beta=result["estimated_within_range_slope"]["native"]["annual_precip_degrees_per_100mm"]
        assert 1<beta<5
        assert result["original_native_trained_introduced_transfer"]["status"].startswith("cross_range")
        assert result["range_BIO12_difference"]["contrasts_are_not_equivalence_tests"]
        assert 0<=result["joint_climatic_support"]["chelsa_bio12"]["native_fraction_in_overlap"]<=1


def test_climate_rank_and_rare_cells_not_falsely_fit_as_absence_of_introduced_slope():
    d=toy_focal(60)
    p=dated_and_balanced(d.loc[(d.taxon_name=="Cirsium arvense") & (d.range_state=="native")])
    for variable in VARIABLES:p[variable]=10.
    sc={name:{"mean":10.,"sd":1.} for name in VARIABLES}
    v=cell_weighted_fit(p,sc)
    assert v["status"]=="five_climate_and_month_design_not_identifiable"
    small=p.head(10)
    assert cell_weighted_fit(small,sc)["status"]=="insufficient_original_pixels_cells_or_focal_orientation"


def test_source_clip_ranges_are_not_implicit_introduction_ages():
    d=toy_focal(240)
    nat=dated_and_balanced(d.loc[(d.taxon_name=="Cirsium arvense") & (d.range_state=="native")])
    intro=dated_and_balanced(d.loc[(d.taxon_name=="Cirsium arvense") & (d.range_state=="introduced")])
    sc=standardization(pd.concat([nat,intro],ignore_index=True))
    v=transfer_against_native_baseline(nat,intro,sc)
    assert "introduced_weighted_rmse" in v
    assert v["climatic_overlap_gated"] is False
    assert "evolution" in v["scope"]


def test_cache_sampler_rejects_wrong_coordinate_grain(tmp_path):
    with pytest.raises(ValueError,match="coordinate mismatch"):
        sample_cached_pixels(tmp_path,{"transform":[1,0,0,0,-1,0]},
                             np.ones(3),np.ones(2))

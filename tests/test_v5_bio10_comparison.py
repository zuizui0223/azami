"""Synthetic tests for BIO05→BIO10 same-time-window V5 I + II analysis."""
import numpy as np
import pandas as pd

from analysis.v5.bio10_s1 import (
    EXTRA, CLIMATE, CLIMATE_NAMES, THRESH, _bootstrap_moments
)
from analysis.v5.bio10_s2 import MODELS
from analysis.v5.bio10_s2_angle_chroma import AXES, _fit, RESPONSES
from analysis.v5.bio10_source_audit import URLS


def test_three_five_axis_models_and_only_one_temperature_replaced():
    assert set(CLIMATE)==set(MODELS)==set(AXES)
    assert len(MODELS)==3
    assert len(EXTRA)==5 and EXTRA[-1]=="chelsa_bio10"
    assert MODELS["annual_plus_seasonality"]==(
        "chelsa_bio01","chelsa_bio04","chelsa_bio12",
        "chelsa_bio15","chelsa_rsds_mean")
    old=MODELS["warm_plus_seasonality"]
    new=MODELS["warm_bio10_plus_seasonality"]
    assert len(old)==len(new)==5
    assert [o for o,n in zip(old,new) if o!=n]==["chelsa_bio05"]
    assert [n for o,n in zip(old,new) if o!=n]==["chelsa_bio10"]
    assert ("chelsa_bio18" in old and "chelsa_bio18" in new)
    assert CLIMATE["warm_bio10_plus_seasonality"]==(
        "_z8","_z10","_z2","_z6","_z7")
    assert set(AXES["warm_bio10_plus_seasonality"])==set(new)
    assert CLIMATE_NAMES["warm_bio10_plus_seasonality"]==(
        "BIO18","BIO10","RSDS","BIO04","BIO15")


def test_chelsa_bio10_official_original_1981_2010_path():
    assert len(URLS)==2
    assert URLS[0].endswith("/bio10/1981-2010/CHELSA_bio10_1981-2010_V.2.1.tif")
    assert "CHELSA_bio10_1981-2010_V.2.1.tif" in URLS[1]


def test_same_s1_taxon_block_bootstrap_dimensions():
    rng=np.random.default_rng(21)
    x=rng.normal(size=(3,8,5))
    y=rng.normal(size=(3,8,3))
    i,j=np.array([0,0]),np.array([1,2])
    mult=np.ones((6,3))
    moments=_bootstrap_moments(x,y,i,j,mult)
    assert set(moments)=={"between_pair","within_pair"}
    for xx,xy in moments.values():
        assert xx.shape==(6,25) and xy.shape==(6,15)
    assert THRESH>0


def test_s2_within_region_angle_chroma_outcomes_both_are_identifiable():
    rng=np.random.default_rng(2026)
    records=[]
    for state in ("native","introduced"):
        for block in range(5):
            for k in range(15):
                rain,temp=rng.normal(size=2)
                seasonality=rng.normal(size=2)
                rsds=rng.normal()
                records.append(dict(
                    range_state=state,block=block,cluster=state+"_"+str(block),
                    chelsa_bio01=temp+rng.normal()*.2,
                    chelsa_bio12=rain+rng.normal()*.2,
                    chelsa_bio05=temp+rng.normal()*.2,
                    chelsa_bio10=temp+rng.normal()*.15,
                    chelsa_bio18=rain+rng.normal()*.1,
                    chelsa_rsds_mean=rsds,
                    chelsa_bio04=seasonality[0],
                    chelsa_bio15=seasonality[1],
                    orientation_angle_degrees_median=65+2*rain+rng.normal(),
                    corolla_lab_chroma_median=30+3*rsds+rng.normal(),
                ))
    d=pd.DataFrame(records)
    for name,axes in AXES.items():
        for _,column in RESPONSES.items():
            fit=_fit(d,axes,demeaning=True,response=column)
            assert fit["n_geographic_region_clusters"]==10
            assert np.isfinite(fit["rain_partial_slope"]["native"]["beta"])
            assert np.isfinite(fit["radiation_partial_slope"]["introduced"]["beta"])
            assert fit["model"]=="within_20degree_region"


def test_s2_bio10_chroma_has_radiation_focal_without_conflating_rarity():
    assert RESPONSES["chroma"]=="corolla_lab_chroma_median"
    from analysis.v5.bio10_s2_angle_chroma import FOCAL
    assert FOCAL["chroma"]=="radiation"
    assert FOCAL["angle"]=="rain"

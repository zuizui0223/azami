"""Source-free tests for native vs introduced climate–head angle slopes."""
import numpy as np
import pandas as pd

from analysis.v5.annual_warm_s2_rain_angle import (
    ANGLE, AXES, _fit, _supported_cells
)
from analysis.v5.annual_warm_s2 import MODELS


def make_synthetic(seed=9):
    rng=np.random.default_rng(seed)
    rows=[]
    for state in ("native","introduced"):
        for block in range(10 if state=="native" else 100,
                           15 if state=="native" else 105):
            for i in range(18):
                x=rng.normal(size=7)
                # True within-block rain effect positive in both ranges, but
                # stronger when introduced. Geographic intercept is a confound.
                response=75+3*x[0]+2*(state=="introduced")*x[0]+(block%5)*4
                response+=rng.normal(scale=.35)
                rows.append({
                    "taxon_name":"Cirsium arvense",
                    "range_state":state,"block":block,"cluster":state+"_"+str(block),
                    "i":block,"j":i,
                    "chelsa_bio12":x[0],
                    "chelsa_bio01":x[1],
                    "chelsa_bio18":x[2],
                    "chelsa_bio05":x[3],
                    "chelsa_rsds_mean":x[4],
                    "chelsa_bio04":x[5],
                    "chelsa_bio15":x[6],
                    ANGLE:response,
                })
    return pd.DataFrame(rows)


def test_same_seasonality_and_radiation_in_both_rain_angle_models():
    assert set(AXES)==set(MODELS)
    for name in AXES:
        assert len(AXES[name])==5
        assert set(AXES[name])==set(MODELS[name])
    annual,warm=AXES.values()
    assert annual[2:]==warm[2:]
    assert annual[:2]==("chelsa_bio12","chelsa_bio01")
    assert warm[:2]==("chelsa_bio18","chelsa_bio05")


def test_supported_20degree_regions_and_one_degree_cells_are_kept():
    c=make_synthetic()
    out=_supported_cells(c,"Cirsium arvense")
    assert len(out)==180
    assert out.groupby("range_state").block.nunique().to_dict()=={
        "native":5,"introduced":5}
    assert out[["range_state","block","i","j"]].duplicated().sum()==0


def test_within_region_native_introduced_direct_rain_angle_slopes():
    c=_supported_cells(make_synthetic(),"Cirsium arvense")
    within=_fit(c,AXES["annual_plus_seasonality"],demeaning=True)
    p=within["rain_partial_slope"]
    assert within["n_geographic_region_clusters"]==10
    assert within["model"]=="within_20degree_region"
    assert p["native"]["beta"]>0.1
    assert p["introduced"]["beta"]>p["native"]["beta"]
    assert p["introduced_minus_native"]["beta"]>0.1
    assert p["introduced_minus_native"]["p_pointwise"]<.05
    pooled=_fit(c,AXES["annual_plus_seasonality"],demeaning=False)
    assert pooled["n_1degree_cells"]==within["n_1degree_cells"]


def test_warm_climate_does_not_relabel_annual_climate_field():
    c=_supported_cells(make_synthetic(),"Cirsium arvense")
    warm=_fit(c,AXES["warm_plus_seasonality"],demeaning=True)
    annual=_fit(c,AXES["annual_plus_seasonality"],demeaning=True)
    assert warm["rain_partial_slope"]["native"]["beta"] != annual["rain_partial_slope"]["native"]["beta"]
    assert warm["n_1degree_cells"]==annual["n_1degree_cells"]

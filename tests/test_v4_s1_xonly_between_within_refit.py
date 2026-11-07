"""Synthetic invariants for original-S1 outcome-blind gradient fitting."""
import numpy as np
import pandas as pd
import pytest
from analysis.v4.s1_xonly_between_within_refit import _sandwich,_inference,selected_axes,_prepare

def test_original_blind_selected_axes_only():
    x={"status":"original_S1_authenticated_X_only_blind_selection",
       "uses_original_or_previous_flower_trait_results":False,
       "uses_photographic_trait_columns":False,
       "original_native_observations":313350,
       "section31":{"between_taxa_X_only":{"selected_variables":["chelsa_bio12"]},
                    "within_taxa_X_only":{"selected_variables":["chelsa_bio18"]}}}
    b,w=selected_axes(x)
    assert b==["chelsa_bio12"] and w==["chelsa_bio18"]
    x["uses_photographic_trait_columns"]=True
    with pytest.raises(ValueError,match="X_ONLY"):selected_axes(x)

def test_sandwich_recover_linear_between_within():
    g=np.repeat(np.arange(40),8)
    w=np.repeat(1/8,320)
    rng=np.random.default_rng(5)
    X=np.column_stack([np.ones(320),rng.normal(size=320),
                       rng.normal(size=320)])
    Y=np.column_stack([1+X[:,1]*.4-X[:,2]*.8+rng.normal(0,.1,320)])
    b,s,c,cond=_sandwich(X,Y,w,g)
    rows=_inference(b,s,c,["intercept","between:one","within:two"],
                    {"between:one","within:two"})
    assert len(rows)==2 and all(0<=r["cluster_robust_joint_p"]<=1 for r in rows)
    assert abs(b[1,0]-.4)<.05 and abs(b[2,0]+.8)<.05

def test_photo_qc_keeps_same_rows():
    n=36*6
    groups=np.repeat(np.arange(36),6)
    d=pd.DataFrame({"accepted_key":groups,
                    "chelsa_bio12":np.linspace(1,4,n),
                    "chelsa_bio18":np.linspace(1,7,n),
                    "corolla_lab_chroma":np.arange(n,dtype=float),
                    "corolla_lab_chroma__eligible":True})
    selected,y,w,meta=_prepare(d,"chroma",["chelsa_bio12"],["chelsa_bio18"])
    assert len(selected)==n and np.isclose(w.sum(),36)

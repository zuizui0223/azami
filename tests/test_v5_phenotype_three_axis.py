import numpy as np
import pandas as pd
from analysis.v5.phenotype_three_axis import (
    AXES,SPACES,load_three_axis,codes,gains,evaluate
)

def test_frozen_axes_and_qc(tmp_path):
    from analysis.v5.novel_head_overlap import T
    raw={"obs_id":["a","b"]}
    for key,col in AXES.values():
        raw[col]=[10.,20.]
        raw[col.removesuffix("_median")+"_n_usable_heads"]=[1,1]
    raw[T[0]]=[40.,181.]
    raw[T[1].removesuffix("_median")+"_n_usable_heads"]=[1,0]
    raw["corolla_lab_lightness_median"]=[20.,30.]
    path=tmp_path/"source.csv"
    pd.DataFrame(raw).to_csv(path,index=False)
    d=load_three_axis(path)
    assert list(AXES)==["chroma","angle","outline"]
    assert list(SPACES)==["chroma","angle","outline","joint_3axis"]
    assert len(SPACES["joint_3axis"])==3
    assert not any("lightness" in c or "hue" in c for c in d)
    assert d.loc[1,T[0]]==181 and d.loc[1,T[1]]==20
    assert np.isnan(d.loc[1,"scan_orientation_image_vertical_angle"])
    assert np.isnan(d.loc[1,"scan_corolla_lab_chroma"])
    assert d.loc[0,"scan_capitulum_outline_aspect_ratio"]==10

def test_novel_joint_combination_within_all_marginal_ranges():
    r=np.array([[i,i,i] for i in range(5)],float)
    x=np.array([[0,4,0],[4,4,4]],float)
    row,score=evaluate(r,np.ones(5),None,x,np.array([0,1]),2,np.array([1]),np.array([0]),False)
    assert row["possible_bins"]==343
    assert row["novel_unobserved"]==1 and row["analogue_unobserved"]==0
    assert row["existing_gain"]==0 and row["new_gain"]==1

def test_existing_gain_partition():
    a=np.array([.6,.4,0.]);b=np.array([.3,.6,.1])
    np.testing.assert_allclose(gains(a,b,np.array([True,True,False])),[.2,.1])

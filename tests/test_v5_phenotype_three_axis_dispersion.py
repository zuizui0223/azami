import numpy as np
from analysis.v5.phenotype_three_axis_dispersion import (
    cell_moments,summarize,group_metrics,reference_scale,block_test,
)

def test_equal_cell_mixture_and_decomposition():
    v=np.array([[0.,0.,0.],[2.,2.,2.],[4.,4.,4.]])
    m,count=cell_moments(v,np.array([0,0,1]),2)
    a,total,between,within=summarize(m.mean(0))
    np.testing.assert_allclose(a,[2.5]*3)
    np.testing.assert_allclose(total,[2.75]*3)
    np.testing.assert_allclose(between,[2.25]*3)
    np.testing.assert_allclose(within,[.5]*3)
    np.testing.assert_array_equal(count,[2,1])

def test_translation_is_not_dispersion():
    v=np.array([[0.,1.,2.],[2.,3.,4.],[10.,11.,12.],[12.,13.,14.]])
    m,_=cell_moments(v,np.arange(4),4)
    np.testing.assert_allclose(group_metrics(m[:2].mean(0)),group_metrics(m[2:].mean(0)))

def test_native_scaling():
    v=np.array([[0.,2.,4.],[2.,4.,8.],[4.,6.,12.]])
    ix=np.array([0,0,1])
    mu,sd=reference_scale(v,ix,2)
    m,_=cell_moments((v-mu)/sd,ix,2)
    np.testing.assert_allclose(summarize(m.mean(0))[1],1.)

def test_permutation_recomputes_group_centroid():
    a=np.repeat(np.linspace(-1,1,8)[:,None],3,axis=1)
    v=np.vstack([a,2*a+3])
    m,_=cell_moments(v,np.arange(16),16)
    pairs=[(np.array([j]),np.array([j+8])) for j in range(8)]
    rows=block_test(m,pairs)
    assert len(rows)==2 and rows[0]["permutations"]==256
    np.testing.assert_allclose(rows[0]["ratio"],4.)
    np.testing.assert_allclose(rows[0]["difference"],rows[1]["difference"])
    assert 0<rows[0]["p"]<=1

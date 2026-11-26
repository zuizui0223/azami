import numpy as np
from analysis.v5.phenotype_reservoir import hist, convex_fit, supported_reweight, tv
def test_hist_joint_equal_cell():
    v=np.array([[0.,2.],[2.,0.],[2.,0.]])
    h=hist(v,np.array([0,0,1]),3,[[.5,1.5],[.5,1.5]])
    np.testing.assert_allclose(h[0,[2,6]],[.5,.5])
    assert h[1,6]==1 and np.isnan(h[2]).all()
def test_mixture():
    a=np.array([[1.,0.,0.],[0.,1.,0.]])
    q=convex_fit(a,np.array([.2,.8,0.]),np.array([.5,.5]))
    np.testing.assert_allclose(q,[.2,.8,0.],atol=1e-6)
    assert tv(q,[.2,.8,0.])<1e-6
def test_no_native_support():
    q,m=supported_reweight(np.array([.6,.4,0.]),np.array([.2,.4,.4]))
    np.testing.assert_allclose(q,[1/3,2/3,0])
    np.testing.assert_allclose(m,.4)

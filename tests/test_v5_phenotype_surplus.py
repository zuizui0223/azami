import numpy as np
from analysis.v5.phenotype_surplus import targets, slope_row
def test_decomposition():
    n=np.array([[.4,.1],[.3,.2]])
    i=np.array([[.2,.3],[.1,.4]])
    common,excess,mass,w=targets(n,i,np.array([0,0,1]),np.array([0,0,1]))
    np.testing.assert_allclose([common.sum(),excess.sum()],[.6,.4])
    np.testing.assert_allclose(mass,[.4,.2])
    np.testing.assert_allclose(w,[[.1,.1,.2],[0,0,.2]])
def test_within_block():
    blocks=np.repeat(np.arange(9),4)
    x=np.tile([-.75,-.25,.25,.75],9)
    y=2*x+blocks*10
    r=slope_row("synthetic","median","angle","mean",y,x,blocks)
    assert r["status"]=="estimated"
    np.testing.assert_allclose(r["beta"],2)
    np.testing.assert_allclose([r["loo_min"],r["loo_max"]],[2,2])

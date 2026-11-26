import numpy as np
import pandas as pd
from analysis.v5.phenotype_observation_audit import center_within,shared_weights,estimate
def test_weighted_demeaning():
    x=np.array([[1.,2.],[3.,4.],[7.,3.],[9.,5.]])
    w=np.array([1.,3.,2.,1.]);g=np.array([0,0,1,1])
    y=center_within(x,w,g)
    for k in [0,1]:
        np.testing.assert_allclose(np.average(y[g==k],axis=0,weights=w[g==k]),0,atol=1e-12)
def test_balanced_bins():
    n=pd.DataFrame(dict(u=[0,1,8],v=[0,1,0],common=True))
    i=pd.DataFrame(dict(u=[0,8,8,16],v=[0,0,1,0],common=True))
    d=shared_weights(n,i)
    np.testing.assert_allclose(d.groupby("introduced").w.sum(),[1,1])
    t=d.assign(bin=d.u//4).groupby(["introduced","bin"]).w.sum().unstack(0)
    np.testing.assert_allclose(t[0],t[1])
def test_regression():
    b=np.repeat(np.arange(9),6);x=np.tile(np.arange(6),9)
    y=2*x+np.sin(np.arange(len(x)))*.1
    r=estimate(y,np.column_stack([np.ones(len(x)),x]),np.ones(len(x)),b,1)
    assert r["status"]=="estimated" and abs(r["beta"]-2)<.02 and r["loo_min"]>1.9

import numpy as np
from analysis.v5.phenotype_flowering_exposure import flow,local_weights

def test_fractional_transport_unequal_cell_counts():
    a=np.array([[0.],[1.]])
    b=np.array([[.1],[.2],[1.1]])
    l,r,m=flow(a,b,np.array([.5,.5]),np.full(3,1/3),(0,))
    np.testing.assert_allclose(m.sum(),5/6,atol=1e-8)
    assert np.all(np.bincount(l,weights=m,minlength=2)<=.5+1e-8)
    assert np.all(np.bincount(r,weights=m,minlength=3)<=1/3+1e-8)
    assert np.all(abs(a[l,0]-b[r,0])<=.52)

def test_no_matching_missing_and_far():
    a=np.array([[np.nan],[0.]])
    b=np.array([[4.]])
    l,r,m=flow(a,b,np.array([.5,.5]),np.ones(1),(0,))
    assert len(m)==0

def test_local_weights_equalize_original_pair_mixture():
    e=np.array([0.,.1,.3,0.,3.,np.nan])
    p=np.array([0,0,1,0,1,1]);y=np.array([0,1,0,0,1,1])
    w=np.array([.5,1.,1.,.5,.5,.5])
    ww,ess,gess,mass=local_weights(e,p,y,w,np.array([0,1]),0.)
    np.testing.assert_allclose(ww[y==0].sum(),1.)
    np.testing.assert_allclose(ww[y==1].sum(),1.)
    assert ww[2]==ww[4]==ww[5]==0
    assert ess==gess==1.

"""Synthetic unit tests only; never used as empirical results."""
import numpy as np
from analysis.v5.phenotype_redistribution import cell_hist,encode,paired_features,weighted_quantile,tv

def test_cell_weights():
    h=cell_hist(np.array([0,0,1,1]),np.array([0,0,0,1]),2,2)
    np.testing.assert_allclose(h,[[2/3,1/3],[0,1]])
    np.testing.assert_allclose(h.mean(0),[1/3,2/3])

def test_joint_categories():
    v=np.array([[0,0],[0,2],[2,0],[2,2]])
    h=encode(v,np.arange(4),4,[np.array([1]),np.array([1])],[0,1])
    np.testing.assert_allclose(h,np.eye(4))

def test_paired_exact():
    mat=np.array([[1.,0.],[0.,1.]]*3)
    pairs=[(np.array([i]),np.array([i+1])) for i in [0,2,4]]
    a,b,d,p,lo,hi,ls,rs,exact,n=paired_features(mat,mat,np.array([.5,.5]),pairs,5)
    assert exact and n==8
    np.testing.assert_allclose(d,[-1,1])
    np.testing.assert_allclose(p,[.25,.25])
    assert tv(a,b)==1

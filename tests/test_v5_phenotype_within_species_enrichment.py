import numpy as np
from analysis.v5.phenotype_within_species_enrichment import (
    K,X,RANKS,native_cuts,rank_codes,cell_hist,fit,fit_coefficients,score_omnibus
)

def test_fixed_native_quintile_encoding_and_outside_gate():
    v=RANKS.copy()
    cuts,lo,hi=native_cuts(v,np.ones(len(v)))
    # Explicit cuts avoid ties in the discrete synthetic native sample.
    cuts=np.tile([-1.5,-.5,.5,1.5],(3,1))
    code,valid=rank_codes(v,cuts,lo,hi)
    np.testing.assert_array_equal(code,np.arange(K))
    assert valid.all()
    h,t,k=cell_hist(np.vstack([v,[10.,0.,0.]]),np.zeros(K+1,int),1,cuts,lo,hi,np.ones(K,bool))
    assert t[0]==126 and k[0]==125
    np.testing.assert_allclose(h[0],np.full(K,1/K))

def test_recover_known_direction_and_curvature_with_pair_intercepts():
    rng=np.random.default_rng(22)
    a=rng.uniform(.2,1,(16,K))
    a/=a.sum(1)[:,None]
    truth=np.array([.22,-.12,.05,-.2,.1,.03])
    b=a*np.exp(X@truth)[None,:]
    b/=b.sum(1)[:,None]
    beta,_,_=fit(a,b,X,np.arange(6))
    np.testing.assert_allclose(beta,truth,atol=1e-6)
    coef,rows=fit_coefficients(a,b,X,np.repeat(np.arange(8),2))
    np.testing.assert_allclose(coef,truth,atol=1e-6)
    assert len(rows)==6
    np.testing.assert_allclose(rows[0]["contrast_ratio"],np.exp(4*truth[0]),atol=1e-6)

def test_identical_groups_have_no_enrichment():
    a=np.tile(np.full(K,1/K),(8,1))
    g=score_omnibus(a,a,X,np.arange(8))
    assert g["p"]==1 and g["score_rank"]==0
    beta,_,_=fit(a,a,X,np.arange(6))
    np.testing.assert_allclose(beta,0)

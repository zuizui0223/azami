import numpy as np
from analysis.v5.phenotype_continuous_enrichment import (
    raw_pairs,state,fit_raw,features,relative_curve,kde_mass,env_profile,X,K,NAMES,RAIN
)
from analysis.v5.phenotype_within_species_enrichment import fit,profile_state

def test_raw_objective_equals_histogram_objective():
    rng=np.random.default_rng(16)
    n=10
    aa=rng.multinomial(80,np.ones(K)/K,size=n)
    bb=rng.multinomial(100,np.ones(K)/K,size=n)
    a,b=aa/80.,bb/100.
    rank=X[:,:3]
    v=np.concatenate([np.repeat(rank,k,axis=0) for k in np.vstack([aa,bb])])
    kept=np.r_[aa.sum(1),bb.sum(1)]
    ix=np.repeat(np.arange(2*n),kept)
    d=raw_pairs(v,ix,np.ones(len(v),bool),np.arange(n),np.arange(n,2*n),
                kept,np.zeros(3),np.ones(3),np.arange(n))
    beta=rng.normal(0,.1,6)
    loss,s,h=state(d,beta)
    l2,s2,h2=profile_state(a,b,X,beta)
    np.testing.assert_allclose(loss,l2,atol=1e-10)
    np.testing.assert_allclose(s,s2,atol=1e-10)
    np.testing.assert_allclose(h,h2,atol=1e-10)
    b1=fit_raw(d,np.arange(6))[0]
    b2=fit(a,b,X,np.arange(6))[0]
    np.testing.assert_allclose(b1,b2,atol=1e-7)

def test_recover_continuous_quadratic():
    rng=np.random.default_rng(72)
    n=25
    z=rng.normal(size=(n*2*180,3))
    p=np.repeat(np.arange(n),360)
    f=features(z)
    w=np.full(len(z),1/180.)
    truth=np.array([.2,-.1,.05,-.16,.12,0.])
    # Generate fractional response using the pair-balanced profile.
    from scipy.special import expit
    eta=f@truth
    for k in range(n):
        mask=p==k;lo,hi=-30.,30.
        for _ in range(52):
            m=(lo+hi)/2
            if (expit(eta[mask]+m)*w[mask]).sum()>1:hi=m
            else:lo=m
        eta[mask]+=(lo+hi)/2
    y=expit(eta)
    d=dict(f=f,p=p,w=w,y=y,blocks=np.arange(n))
    b,_,_=fit_raw(d,np.arange(6))
    np.testing.assert_allclose(b,truth,atol=1e-6)
    r=relative_curve(0,0,np.zeros(3),np.ones(3),b,np.eye(6),1.96)
    assert r["relative_enrichment"]==r["lo"]==r["hi"]==1.

def test_kde_equal_cell_weight():
    v=np.array([0.,0.,0.,1.]);ix=np.array([0,0,0,1])
    z=kde_mass(v,ix,np.array([3,1]),np.array([0,1]),np.array([0.,1.]),.25)
    np.testing.assert_allclose(z[0,0],1.)
    np.testing.assert_allclose(z[1,1],1.)
    np.testing.assert_allclose(z[0,1],z[1,0])


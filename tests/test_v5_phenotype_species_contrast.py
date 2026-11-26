import numpy as np
from analysis.v5.phenotype_species_contrast import X,K,common_support,direct_contrast,balanced_diagnostic

def synthetic(beta,blocks):
    rng=np.random.default_rng(8)
    a=rng.uniform(.5,1.,(len(blocks),K))
    a/=a.sum(1)[:,None]
    b=a*np.exp(X@np.asarray(beta))[None,:]
    b/=b.sum(1)[:,None]
    return a,b,np.asarray(blocks)

def test_difference():
    blocks=np.repeat(np.arange(8),3)
    ca=synthetic(np.zeros(6),blocks)
    cv=synthetic([.2,-.1,.03,-.1,.04,.02],blocks)
    r=direct_contrast(ca,cv)
    np.testing.assert_allclose([x["cv_minus_ca"] for x in r],[.2,-.1,.03,-.1,.04,.02],atol=1e-6)
    assert r[0]["shared_blocks"]==8

def test_shared_counts_and_gate():
    ca=synthetic(np.zeros(6),np.repeat(np.arange(4),2))
    cv=synthetic(np.zeros(6),np.repeat(np.arange(2,10),3))
    ids,n=common_support(ca,cv)
    np.testing.assert_array_equal(ids,[2,3])
    np.testing.assert_array_equal(n,[2,2])
    info,rows=balanced_diagnostic(ca,cv)
    assert info["balanced_pairs_per_species"]==4
    assert all(r["status"]=="insufficient_support" for r in rows)

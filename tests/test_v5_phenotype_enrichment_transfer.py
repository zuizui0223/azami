import numpy as np
from analysis.v5.phenotype_enrichment_transfer import cross_validate,summarize_transfer,X,K

def test_transfer_of_shared_enrichment():
    rng=np.random.default_rng(8)
    a=rng.uniform(.5,1.,(12,K));a/=a.sum(1)[:,None]
    truth=np.array([.2,-.1,.03,-.1,.04,.02])
    b=a*np.exp(X@truth)[None,:];b/=b.sum(1)[:,None]
    base,betas,scores=cross_validate(a,b,np.repeat(np.arange(6),2))
    np.testing.assert_allclose(base,truth,atol=1e-6)
    np.testing.assert_allclose(betas,np.tile(truth,(6,1)),atol=1e-6)
    r=summarize_transfer(scores)
    assert r["full"]["positive_blocks"]==6
    assert r["full"]["loss_reduction_pct"]>0
    assert r["drop_chroma"]["full_advantage_pct_of_null"]>0

def test_no_enrichment():
    a=np.tile(np.full(K,1/K),(6,1))
    _,betas,scores=cross_validate(a,a,np.arange(6))
    np.testing.assert_allclose(betas,0,atol=1e-10)
    r=summarize_transfer(scores)
    np.testing.assert_allclose(r["full"]["loss_reduction_pct"],0,atol=1e-9)
    assert r["full"]["positive_blocks"]==0

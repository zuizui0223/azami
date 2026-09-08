import numpy as np
import pytest
from analysis.v3.joint_partial_pooling import JointSlopeLikelihood


@pytest.mark.parametrize('variance',[0.,.1,1.,10.])
def test_nine_predictor_contractions_match_literal_reference(variance):
    rng = np.random.default_rng(42112)
    groups = np.repeat(np.arange(16),12)
    x = rng.normal(size=(len(groups),9))
    b = rng.normal(size=(len(groups),10))
    y = rng.normal(size=len(groups))
    model = JointSlopeLikelihood(y,x,groups,b)
    out = model.evaluate(np.full(9,variance))
    m = out['conditional_random_covariance_ratio']
    info = model.ff - np.einsum('tpi,tpq,tqj->ij',model.xf,m,model.xf)
    score = model.fy - np.einsum('tpi,tpq,tq->i',model.xf,m,model.xy)
    q = model.yy - np.einsum('tp,tpq,tq->',model.xy,m,model.xy)
    inverse = np.linalg.inv(info)
    alpha = inverse@score
    rss = q-score@alpha
    np.testing.assert_allclose(out['alpha'],alpha,rtol=1e-9,atol=1e-9)
    np.testing.assert_allclose(out['rss'],rss,rtol=1e-11,atol=1e-9)
    vxf = model.xf-model.gram@m@model.xf
    residual_score = model.xy-np.einsum('tpi,i->tp',model.xf,alpha)
    modes = np.einsum('tpq,tq->tp',m,residual_score)
    vxresid = residual_score-np.einsum('tpq,tq->tp',model.gram,modes)
    trace = np.diagonal(model.gram-model.gram@m@model.gram,axis1=1,axis2=2).sum(axis=0)
    trace -= np.einsum('tpi,ij,tpj->p',vxf,inverse,vxf)
    np.testing.assert_allclose(out['gradient'],trace-model.df/rss*(vxresid**2).sum(axis=0),atol=1e-8,rtol=1e-8)

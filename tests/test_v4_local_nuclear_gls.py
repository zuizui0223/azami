import numpy as np
import pytest
import statsmodels.api as sm
from statsmodels.multivariate.manova import MANOVA
from analysis.v4 import local_nuclear_gls as m

def inputs():
    rng=np.random.default_rng(51);e=rng.normal(size=(24,4))
    v=.2*np.ones((24,24))+.8*np.eye(24)
    y=e[:,0]+rng.normal(size=24)
    return e,v,y

def test_scalar_gls_matches_statsmodels_with_all_four_predictors():
    e,v,y=inputs();r=m.fit(y[:,None],e,v)
    x=sm.add_constant((e-e.mean(0))/e.std(0,ddof=1))
    reference=sm.GLS(y,x,sigma=v).fit()
    assert np.allclose(r['coefficients'],reference.params[:,None])
    assert np.allclose([t['p_value'] for t in r['tests']],reference.pvalues[1:])
    assert r['residual_df']==19

def test_paired_hue_joint_test_matches_whitened_manova_and_rotation():
    e,v,y=inputs();z=np.column_stack([np.sin(y),np.cos(y)])
    r=m.fit(z,e,v);l=np.linalg.cholesky(v)
    x=sm.add_constant((e-e.mean(0))/e.std(0,ddof=1))
    for j in range(1,5):
        result=MANOVA(np.linalg.solve(l,z),np.linalg.solve(l,x)).mv_test([('axis',np.eye(5)[j:j+1])])
        assert r['tests'][j-1]['p_value']==pytest.approx(result.results['axis']['stat'].loc["Hotelling-Lawley trace",'Pr > F'])
    q=np.array([[.6,-.8],[.8,.6]])
    rotated=m.fit(z@q,e,v)
    assert np.allclose([t['p_value'] for t in r['tests']],[t['p_value'] for t in rotated['tests']])

def test_singular_climate_design_rejected():
    e,v,y=inputs();e[:,1]=e[:,0]
    with pytest.raises(ValueError):m.fit(y[:,None],e,v)

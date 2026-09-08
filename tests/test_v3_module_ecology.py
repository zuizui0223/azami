import numpy as np
import pytest
from analysis.v3.dependence_resampling import source_partition
from analysis.v3.module_ecology import fit_matched_module,quadratic_geometry,quadratic_values,summarize_draws,bootstrap_matched_module


def synthetic(seed=1928):
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(18),24)
    x = rng.normal(size=(len(g),2))+rng.normal(size=(18,2))[g]
    b = rng.normal(size=(len(g),2))
    y = np.empty(len(g))
    beta_among = np.array([.4,-.3])
    for i in range(18):
        ix = np.flatnonzero(g==i)
        centered = x[ix]-x[ix].mean(axis=0)
        e = rng.normal(size=len(ix)); e -= e.mean()
        y[ix] = centered@(np.array([.2,.7])+rng.normal(0,.2,2)) + x[ix].mean(axis=0)@beta_among + b[ix]@np.array([.1,.5])+e
    return np.column_stack([y,2*y+3]),x,g,b,beta_among


def test_within_and_among_use_identical_taxa_units_and_response_rows():
    y,x,g,b,beta = synthetic()
    result = fit_matched_module(y,x,g,b)
    np.testing.assert_allclose(result['coefficients'][1],np.array([beta,2*beta]),atol=1e-11)
    np.testing.assert_allclose(result['coefficients'][0,1],2*result['coefficients'][0,0],atol=1e-7)
    np.testing.assert_allclose(result['coefficients'][2],result['coefficients'][1]-result['coefficients'][0])
    assert result['taxa']==tuple(sorted(map(str,range(18))))
    assert all(d['within_retained_taxa']==18 for d in result['diagnostics'])


def test_coordinate_matched_nuisance_is_used_at_both_scales():
    y,x,g,b,_ = synthetic()
    repeated = np.stack([b,b],axis=1)
    one = fit_matched_module(y,x,g,b)
    two = fit_matched_module(y,x,g,repeated)
    np.testing.assert_allclose(one['coefficients'],two['coefficients'],atol=1e-12)


def test_among_uses_mean_observation_nuisance_not_function_of_mean_location():
    y,x,g,b,_ = synthetic()
    # A nonlinear supplied coordinate is averaged exactly as supplied.
    b = np.column_stack([b,b[:,0]**2])
    result = fit_matched_module(y,x,g,b)
    labels = sorted(set(g))
    xm = np.array([x[g==i].mean(axis=0) for i in labels])
    ym = np.array([y[g==i].mean(axis=0) for i in labels])
    bm = np.array([b[g==i].mean(axis=0) for i in labels])
    dense = np.linalg.lstsq(np.column_stack([np.ones(len(labels)),xm,bm]),ym,rcond=None)[0][1:3].T
    np.testing.assert_allclose(result['coefficients'][1],dense,atol=1e-10)


def test_rank_aware_statistic_is_invariant_to_units_and_exact_redundancy():
    c = np.array([[2.,.3],[.3,1.]])
    beta = np.array([.2,-.7])
    expected = beta@np.linalg.solve(c,beta)
    assert quadratic_values(beta,quadratic_geometry(c))[0] == pytest.approx(expected)
    units = np.array([1000.,.001])
    assert quadratic_values(beta*units,quadratic_geometry(c*units[:,None]*units[None,:]))[0] == pytest.approx(expected)
    mapping = np.array([[1,0],[0,1],[2,-3]])
    assert quadratic_values(mapping@beta,quadratic_geometry(mapping@c@mapping.T))[0] == pytest.approx(expected)
    with pytest.raises(ValueError,match='outside'):
        quadratic_values(mapping@beta+[0,0,.1],quadratic_geometry(mapping@c@mapping.T))


def test_negative_or_rank_zero_covariance_is_not_a_probability():
    for c in (np.zeros((2,2)),np.array([[1,2],[2,1]])):
        with pytest.raises(ValueError):
            quadratic_geometry(c)


def test_same_draw_covariance_preserves_cross_scale_and_coordinate_dependence():
    rng = np.random.default_rng(667)
    within = rng.normal(size=(120,2,3))
    among = .8*within+rng.normal(0,.1,size=within.shape)
    draws = np.stack([within,among,among-within],axis=1)
    point = np.zeros((3,2,3))
    out = summarize_draws(point,draws,{'a':[0,1],'b':[2]})
    covariance = out['covariance']
    w,a,d = slice(0,6),slice(6,12),slice(12,18)
    np.testing.assert_allclose(covariance[d,d],covariance[a,a]+covariance[w,w]-covariance[a,w]-covariance[w,a],atol=1e-14)
    assert np.trace(covariance[d,d]) < .1*np.trace(covariance[a,a]+covariance[w,w])
    assert len(out['candidate_tests'])==6
    assert out['ecological_fitting_authorized'] is False


def test_failed_bootstrap_draw_is_retained_without_redraw_or_success_conditioning():
    n = 36
    source = source_partition(range(n),np.repeat(np.arange(12),3),range(n),
                              np.tile([1,11,21],12),np.tile([2,12,22],12),grid_degrees=2)
    calls = []
    def fake_fit(y,x,t,b):
        calls.append(1)
        if len(calls)==3:
            raise ValueError('Deliberate nonestimable draw')
        return {'coefficients':np.array([[[y.mean()]],[[2*y.mean()]],[[y.mean()]]]),
                'taxa':tuple(set(t)),'diagnostics':[]}
    out = bootstrap_matched_module(np.arange(n)[:,None],np.ones((n,1)),np.empty((n,0)),
                                   source,np.arange(n),seed=77,replicates=5,blocks={'a':[0]},fit=fake_fit)
    assert len(calls)==6
    assert out['planned_replicates']==5 and out['estimable_replicates']==4
    assert out['records'][1]['status']=='not_estimable'
    assert out['summary'] is None
    assert np.isnan(out['draws'][1]).all()


def test_nonfinite_module_support_and_among_aliases_fail_without_dropping_columns():
    y,x,g,b,_ = synthetic()
    damaged = y.copy(); damaged[0,0] = np.nan
    with pytest.raises(ValueError,match='aligned'):
        fit_matched_module(damaged,x,g,b)
    with pytest.raises(ValueError,match='Among exposure'):
        fit_matched_module(y,x,g,np.column_stack([b,x[:,0]]))

import numpy as np
import pandas as pd
import pytest
import json
from pathlib import Path
from analysis.v4.observation_models import reuse_inventory, spherical_basis
from analysis.v4.observation_models import ENV, MODULES, fit_response, fit_module, design_diagnostics


def test_crossproduct_fit_matches_dense_gls_without_averaging():
    rng=np.random.default_rng(942)
    groups=np.repeat(np.arange(25),7)
    x=np.column_stack([np.ones(len(groups)),rng.normal(size=len(groups))])
    y=x@np.array([.3,.7])+np.repeat(rng.normal(size=25),7)+rng.normal(size=len(groups))*.4
    fit=fit_response(x,y,groups)
    assert fit['status']=='estimated'
    v=np.eye(len(y))+fit['variance_ratio']*(groups[:,None]==groups[None,:])
    info=x.T@np.linalg.solve(v,x)
    beta=np.linalg.solve(info,x.T@np.linalg.solve(v,y))
    assert fit['beta']==pytest.approx(beta,abs=1e-9)
    assert fit['n']==len(y)
    assert np.sum(fit['influence'],axis=0)==pytest.approx(np.zeros(2),abs=1e-8)


def test_scalar_reml_agrees_with_statsmodels():
    from statsmodels.regression.mixed_linear_model import MixedLM
    rng=np.random.default_rng(20261020)
    g=np.repeat(np.arange(30),9)
    x=np.column_stack([np.ones(len(g)),rng.normal(size=len(g))])
    y=.4*x[:,1]+rng.normal(size=30)[g]*.8+rng.normal(size=len(g))*.4
    actual=fit_response(x,y,g)
    reference=MixedLM(y,x,groups=g).fit(reml=True,method='lbfgs')
    assert reference.converged
    assert actual['beta']==pytest.approx(reference.fe_params,abs=1e-4)
    assert actual['variance_ratio']==pytest.approx(float(reference.cov_re[0,0]/reference.scale),rel=2e-3)


def sample():
    rng=np.random.default_rng(300)
    g=np.repeat(np.arange(55),12)
    frame=pd.DataFrame({'accepted_key':g})
    scales={}
    for v in ENV:
        b=rng.normal(size=55)[g];w=rng.normal(size=len(g))
        w-=pd.Series(w).groupby(g).transform('mean').to_numpy()
        frame[v+'__between']=b;frame[v+'__within']=w;frame[v]=b+w
        scales[v]={'mean':0.,'sd':1.}
    signal=1.5*frame[ENV[0]+'__between']-1.7*frame[ENV[0]+'__within']
    for outcome in ['y1','y2']:
        frame[outcome]=signal+rng.normal(size=55)[g]*.4+rng.normal(size=len(g))*.5
        frame[outcome+'__eligible']=True
    return frame,scales


def test_partial_colour_coordinates_retained_and_scales_distinguished():
    frame,scales=sample()
    frame.loc[:24,'y2']=np.nan
    fit=fit_module(frame,['y1','y2'],scales)
    assert fit['status']=='estimated'
    assert fit['components']['y1']['n']==660
    assert fit['components']['y2']['n']==635
    effects=fit['components']['y1']['effects']
    assert effects[ENV[0]+'__between']['raw_units_per_environment_sd']>1.
    assert effects[ENV[0]+'__within']['raw_units_per_environment_sd']<-1.
    assert fit['tests'][ENV[0]+'__within']['df1']==2
    assert fit['tests'][ENV[0]+'__within']['p']<.001


def test_joint_test_unchanged_by_response_units():
    frame,scales=sample()
    before=fit_module(frame,['y1','y2'],scales)
    frame['y2']=frame.y2*100+40
    after=fit_module(frame,['y1','y2'],scales)
    for term,test in before['tests'].items():
        assert after['tests'][term]['F']==pytest.approx(test['F'],rel=1e-6)
    for term,test in before['joint_temperature_vpd'].items():
        assert test['df1']==4
        assert after['joint_temperature_vpd'][term]['F']==pytest.approx(test['F'],rel=1e-6)


def test_singular_design_is_not_promoted_to_inference():
    x=np.ones((80,3))
    assert fit_response(x,np.arange(80),np.repeat(np.arange(20),4))['status']=='not_identifiable'


def test_actual_design_diagnostics_preserve_endpoint_missingness_and_detect_collinearity():
    frame,scales=sample()
    for outcomes in MODULES.values():
        for outcome in outcomes:
            frame[outcome]=frame.y1
            frame[outcome+'__eligible']=True
    first=MODULES['colour'][0]
    frame.loc[:19,first+'__eligible']=False
    a=design_diagnostics(frame)['endpoints']
    assert a[first]['full_eight_terms']['n']==640
    assert a[MODULES['orientation'][0]]['full_eight_terms']['n']==660
    assert a[first]['full_eight_terms']['status']=='estimated'
    frame[ENV[2]+'__between']=frame[ENV[0]+'__between']
    b=design_diagnostics(frame)['endpoints'][first]
    assert b['full_eight_terms']['status']=='singular'
    assert b['between_observation_weighted']['status']=='singular'
    assert b['within']['status']=='estimated'


def test_omission_sensitivity_does_not_recover_missing_predictor_rows():
    frame,scales=sample()
    frame.loc[:15,ENV[2]]=np.nan
    full=fit_module(frame,['y1','y2'],scales)
    reduced=fit_module(frame,['y1','y2'],scales,[v for v in ENV if v!=ENV[2]])
    assert full['components']['y1']['n']==644
    assert reduced['components']['y1']['n']==644
    assert reduced['joint_temperature_vpd']=={}
    assert ENV[2]+'__between' not in reduced['tests']


def test_inventory_keeps_unmeasured_and_ineligible_values_out_of_inference():
    frame,_=sample()
    crosswalk=json.loads((Path(__file__).parents[1]/'analysis/v4/trait27_functional_crosswalk_v1.json').read_text())
    report={'endpoints':{}}
    for entry in crosswalk['endpoints']:
        name=entry['id']
        frame[name]=np.nan
        frame[name+'__eligible']=False
        report['endpoints'][name]={'finite':0,'eligible':0,
                                  'eligible_with_four_climate':0,'source':'fixture'}
    main='corolla_lab_lightness'
    secondary='capitulum_outline_aspect_ratio'
    technical='involucre_surface_edge_density'
    for name in (main,secondary,technical):
        frame[name]=frame.y1
        report['endpoints'][name]['finite']=len(frame)
        if name!=technical:
            frame[name+'__eligible']=True
            report['endpoints'][name].update(eligible=len(frame),eligible_with_four_climate=len(frame))
    rows={r['endpoint']:r for r in reuse_inventory(frame,report,crosswalk)}
    assert len(rows)==27
    assert rows[main]['current_use']=='joint_primary_module'
    assert rows[secondary]['current_use']=='secondary_phenotype_description'
    assert rows[technical]['current_use']=='technical_validation_only'
    assert rows[technical]['eligible_median'] is None
    assert rows['corolla_hue_sin']['current_use']=='unmeasured'
    report['endpoints'][main]['eligible']=1
    with pytest.raises(ValueError,match='coverage differs'):
        reuse_inventory(frame,report,crosswalk)


def test_spherical_basis_is_periodic_across_longitude_and_has_eight_terms():
    lat=np.array([-40.,0.,35.,65.])
    lon=np.array([-179.,-40.,80.,179.])
    assert spherical_basis(lat,lon).shape==(4,8)
    assert spherical_basis(lat,lon)==pytest.approx(spherical_basis(lat,lon+360),abs=1e-12)


def test_fixed_adjustment_removes_known_confounding_without_averaging_rows():
    frame,scales=sample()
    rng=np.random.default_rng(822)
    frame['adjustment']=frame[ENV[0]+'__within']+rng.normal(size=len(frame))
    frame['y1']=2*frame.adjustment+rng.normal(size=len(frame))*.02
    reference=fit_module(frame,['y1'],scales)
    adjusted=fit_module(frame,['y1'],scales,adjustments=['adjustment'])
    term=ENV[0]+'__within'
    assert reference['components']['y1']['effects'][term]['raw_units_per_environment_sd']>1.8
    assert abs(adjusted['components']['y1']['effects'][term]['raw_units_per_environment_sd'])<.01
    assert adjusted['components']['y1']['n']==reference['components']['y1']['n']==660
    assert len(adjusted['tests'])==8
    assert 'adjustment' not in adjusted['tests']


def test_added_environment_uses_finite_matched_rows():
    frame,scales=sample()
    rng=np.random.default_rng(542)
    name='seasonality'
    frame[name]=rng.normal(size=len(frame))
    frame[name+'__between']=frame.groupby('accepted_key')[name].transform('mean')
    frame[name+'__within']=frame[name]-frame[name+'__between']
    scales[name]={'mean':frame[name].mean(),'sd':frame[name].std(ddof=0)}
    frame.loc[0,name]=np.nan
    expanded=fit_module(frame,['y1'],scales,ENV+[name])
    reference=fit_module(frame.loc[frame[name].notna()],['y1'],scales)
    assert expanded['status']=='estimated'
    assert expanded['components']['y1']['n']==reference['components']['y1']['n']==659
    assert len(expanded['tests'])==10


def test_selected_design_does_not_require_excluded_environment():
    frame,scales=sample()
    frame.loc[0,ENV[2]]=np.nan
    selected=[v for v in ENV if v!=ENV[2]]
    model=fit_module(frame,['y1'],scales,selected,required_environment=selected)
    assert model['components']['y1']['n']==660
    assert ENV[2]+'__within' not in model['tests']


def test_common_slope_uses_raw_predictors_with_species_random_intercept():
    frame,scales=sample()
    model=fit_module(frame,['y1'],scales,ENV,required_environment=ENV,decompose=False)
    x=np.column_stack([np.ones(len(frame))]+[
        (frame[v]-scales[v]['mean'])/scales[v]['sd'] for v in ENV])
    y=frame.y1.to_numpy()
    reference=fit_response(x,(y-y.mean())/y.std(),frame.accepted_key.to_numpy())
    effects=model['components']['y1']['effects']
    assert set(effects)==set(ENV)
    assert [effects[v]['standardized_beta'] for v in ENV]==pytest.approx(reference['beta'][1:])
    assert model['components']['y1']['n']==660

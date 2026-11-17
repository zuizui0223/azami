"""All data generated here are synthetic; no original V4 S1 individual records."""
import json
import numpy as np
import pandas as pd
import pytest
from sklearn.cross_decomposition import PLSRegression

from analysis.v4 import species_environment_spls2 as s


def rows(n=80,obs=5,seed=718):
    rng=np.random.default_rng(seed)
    xs=rng.normal(size=(n,len(s.CLIMATE)))
    X=np.repeat(xs,obs,axis=0)
    X+=rng.normal(scale=.025,size=X.shape)
    y=np.column_stack([
        xs[:,0]*2+xs[:,4]*.4+rng.normal(size=n)*.1,
        xs[:,0]+rng.normal(size=n)*.1,
        xs[:,3]+rng.normal(size=n)*.1,
        np.sin(xs[:,0]),np.cos(xs[:,0]),
        xs[:,6]+rng.normal(size=n)*.25,
        xs[:,5]+rng.normal(size=n)*.25,
        xs[:,7]+rng.normal(size=n)*.25,
        xs[:,8]+rng.normal(size=n)*.25,
    ])
    Y=np.repeat(y,obs,axis=0)+rng.normal(scale=.05,size=(n*obs,len(s.YCOLS)))
    df=pd.DataFrame(X,columns=s.CLIMATE)
    for j,c in enumerate(s.YCOLS):
        df[c]=Y[:,j]
        df[c+'__eligible']=True
    df['accepted_key']=np.repeat([f'ACCEPTED_{i:03}' for i in range(n)],obs)
    df['obs_id']=[f'photo_{i:06}' for i in range(len(df))]
    df['longitude']=np.repeat(np.linspace(-170,150,n),obs)
    return df


def test_joint_original_source_photo_qc_and_hue_unit():
    d=rows(n=80,obs=6)
    d.loc[:2,'corolla_lab_chroma__eligible']=False
    d.loc[6:11,'capitulum_outline_circularity__eligible']=False
    a,c=s.aggregate_species_joint(d)
    assert c['n_joint_supported_taxa']==78
    assert c['all_eight_joint_photo_eligible']==len(d)-9
    assert c['source_photo_contributions_to_taxon_means']==len(d)-12
    assert np.allclose(np.linalg.norm(a[list(s.FEATURES['hue'])].to_numpy(float),axis=1),1)
    assert len(c['predictor_support'])==22
    assert 'topo_elevation' not in a


def test_order_independent_qc_aggregation():
    d=rows(n=58,obs=5)
    a,_=s.aggregate_species_joint(d)
    b,_=s.aggregate_species_joint(d.sample(frac=1,random_state=2).reset_index(drop=True))
    pd.testing.assert_frame_equal(a,b,check_exact=False,atol=1e-12)


def test_no_spatial_eight_cell_gate():
    d=rows(n=70,obs=5)
    d['longitude']=20.
    a,c=s.aggregate_species_joint(d)
    assert len(a)==70
    assert len(set(a.sector))==1


def test_qualified_phenotype_missing_X_still_keeps_other_exposures():
    d=rows(n=50,obs=5)
    d.loc[d.accepted_key=='ACCEPTED_000',s.CLIMATE[0]]=np.nan
    d.loc[d.accepted_key=='ACCEPTED_001',s.CLIMATE[1]]=np.nan
    a,c=s.aggregate_species_joint(d)
    assert len(a)==50
    assert np.isnan(a.loc[a._taxon_private=='ACCEPTED_000',s.CLIMATE[0]]).all()
    assert np.isfinite(a.loc[a._taxon_private=='ACCEPTED_000',s.CLIMATE[1]]).all()
    assert c['predictor_support'][0]['n_taxa_with_finite_species_mean']==49


def test_duplicate_ids_abort_before_model():
    d=rows(n=50,obs=5)
    d.loc[1,'obs_id']=d.loc[0,'obs_id']
    with pytest.raises(ValueError,match='IDENTITIES'):s.aggregate_species_joint(d)


def test_explicit_x_sparsity_and_soft_threshold():
    rng=np.random.default_rng(22)
    c=rng.normal(size=(22,9))
    for k in [2,4,7,11,16,22]:
        w=s._sparse_weight(c,k)
        assert np.isclose(np.linalg.norm(w),1)
        assert 1 <= np.count_nonzero(w) <= k


def test_fit_and_predict_train_scale_not_test_scale():
    rng=np.random.default_rng(2)
    X=rng.normal(size=(70,22));Y=np.column_stack([X[:,0],X[:,1],np.sin(X[:,0]),np.cos(X[:,0]),X[:,3],X[:,4],X[:,5],X[:,6],X[:,7]])
    m=s.SparsePLS2(4,2).fit(X,Y)
    v=m.predict(X[:10]);assert v.shape==(10,9)
    assert np.isfinite(v).all()
    assert len(m.selected)<=8
    shifted=X[:10].copy();shifted[:,2]=np.nan
    assert np.isfinite(m.predict(shifted)).all()


def test_first_component_all_predictors_equivalent_to_standard_pls1_score():
    rng=np.random.default_rng(123)
    X=rng.normal(size=(90,22));Y=rng.normal(size=(90,9))
    m=s.SparsePLS2(22,1).fit(X,Y)
    # Compare to PLSRegression after identical train scaling and multivariate Y weights.
    xx=(X-m.xmean)/m.xstd;yy=(Y-m.ymean)/m.yscale
    sk=PLSRegression(n_components=1,scale=False).fit(xx,yy)
    p1=(xx@m.beta);p2=sk.predict(xx)
    assert np.corrcoef(p1.ravel(),p2.ravel())[0,1]>.97


def test_synthetic_taxon_signal_nested_cv_and_observer_free_output():
    d=rows(n=80,obs=5)
    result=s.run_spls2(d,bootstrap=6,k_grid=(2,7,22),h_grid=(1,2))
    assert result['species_support']['n_joint_supported_taxa']==80
    assert result['taxon_holdout_nested_cv']['skill_all_eight']>.25
    assert result['taxon_holdout_nested_cv']['n_outer_partitions']==10
    assert result['three_sector_holdout_nested_cv']['status']=='estimated_three_dominant_longitude_sectors'
    assert len(result['environment_selection_stability'])==22
    assert result['final_fit']['hyperparameters']['keepX'] in (2,7,22)
    text=json.dumps(result)
    for secret in ('ACCEPTED_000','photo_000000','individual_prediction'):
        assert secret not in text


def test_no_confirmatory_q_and_labelled_exploratory():
    d=rows(n=55,obs=5)
    result=s.run_spls2(d,bootstrap=3,k_grid=(2,22),h_grid=(1,))
    assert result['causality_or_adaptation_established'] is False
    assert result['no_new_confirmatory_p_values'] is True
    assert result['uses_soil'] is False
    assert result['includes_elevation_or_terrain'] is False
    assert len(result['final_fit']['axes'])==1
    assert result['bootstrap_success_count']==3


def test_trait_hue_rotation_equivariance_of_regression_fitting():
    rng=np.random.default_rng(391)
    X=rng.normal(size=(100,22));Y=rng.normal(size=(100,9))
    Y[:,s.OUTCOMES['hue']]=rng.normal(size=(100,2))
    angle=.71
    rot=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    alt=Y.copy();alt[:,s.OUTCOMES['hue']]=Y[:,s.OUTCOMES['hue']]@rot
    a=s.SparsePLS2(7,2).fit(X,Y)
    b=s.SparsePLS2(7,2).fit(X,alt)
    assert np.isclose(a.response_loss(X,Y)[0],b.response_loss(X,alt)[0],rtol=1e-4)
    assert set(a.selected)==set(b.selected)


def test_no_source_variables_is_not_faked():
    d=rows(n=55,obs=5)
    d[s.CLIMATE[5]]=np.nan
    with pytest.raises(ValueError,match='CLIMATE_AXIS_MISSING'):s.run_spls2(d,bootstrap=1,k_grid=(2,),h_grid=(1,))

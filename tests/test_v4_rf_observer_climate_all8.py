"""Synthetic-only tests: no protected source data."""
import json
import hashlib
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

MODULE = Path(__file__).resolve().parents[1] / 'analysis/v4/rf_observer_climate_all8.py'

def get_module():
    assert MODULE.exists(), 'all-eight observer/climate implementation missing'
    spec=importlib.util.spec_from_file_location('observer_climate', MODULE)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def fixture_frame():
    # Globally reused observer IDs; trait is not allowed to determine assignment.
    m=get_module();rng=np.random.default_rng(9);rows=[]
    for t in range(6):
        for c in range(16):
            x=(c-7.5)/4
            for u in range(24):
                d={'obs_id':f'p{t}_{c}_{u}', 'user_id':f'u{u}',
                   'accepted_key':f't{t}', 'latitude':-42.1+c*5.1,
                   'longitude':12.2+t*.03}
                for j,e in enumerate(m.PREDICTORS):d[e]=x*(j+1)+rng.normal(0,.01)
                for name,cols in m.FEATURES.items():
                    for k,col in enumerate(cols):
                        d[col]=(np.sin(.2*x) if k==0 else np.cos(.2*x)) if name=='hue' else 3*t+4*x+rng.normal(0,.04)
                        d[col+'__eligible']=True
                rows.append(d)
    return pd.DataFrame(rows)


def test_implementation_exists():
    get_module()


def test_partition_global_and_one_observer_per_taxon_cell():
    m=get_module();f=fixture_frame()
    duplicate=f.iloc[[0]].copy();duplicate.obs_id='second_photo'
    pairs,audit=m.make_pairs(pd.concat([f,duplicate]),'angle',1,20261101)
    assert audit['observer_overlap']==0
    assert audit['selected_records']==len(f)
    assert audit['candidate_records']==len(f)+1
    assert pairs.shape[0]==96 and pairs.taxon.nunique()==6
    assert (pairs[['n_a','n_b']]>=3).all().all()
    assert set(pairs.columns)>=set(c+'_a' for c in m.PREDICTORS)


def test_order_invariant_and_qc_applied():
    m=get_module();f=fixture_frame()
    p,_=m.make_pairs(f,'angle',1,20261101)
    q,_=m.make_pairs(f.sample(frac=1,random_state=4),'angle',1,20261101)
    pd.testing.assert_frame_equal(p,q)
    f.loc[:,'orientation_image_vertical_angle__eligible']=False
    p,a=m.make_pairs(f,'angle',1,20261101)
    assert p.empty and a['candidate_records']==0


def test_hue_is_joint_unit_vector():
    m=get_module();p,_=m.make_pairs(fixture_frame(),'hue',1,20261101)
    for side in ('a','b'):
        y=p[[c+'_'+side for c in m.FEATURES['hue']]].to_numpy()
        assert np.allclose(np.linalg.norm(y,axis=1),1)


def test_under_supported_observers_and_unknown_endpoint():
    m=get_module();f=fixture_frame();f=f[f.user_id.isin(['u1','u2'])]
    p,_=m.make_pairs(f,'angle',1,20261101);assert p.empty
    with pytest.raises(ValueError):m.make_pairs(f,'bad_trait',1,20261101)
    with pytest.raises(ValueError):m.make_pairs(f,'angle',2,20261101)


def test_training_fit_does_not_see_evaluation_outcomes_or_coordinates():
    m=get_module();p,_=m.make_pairs(fixture_frame(),'angle',1,20261101)
    tr=p[p.lat_cell<88].copy();te=p[p.lat_cell>=88].copy()
    y=m.FEATURES['angle']
    a=m.fit_fold(tr,te,y,'a','b',False,trees=16,seed=14)
    altered=te.copy();altered[y[0]+'_b']+=10000
    b=m.fit_fold(tr,altered,y,'a','b',False,trees=16,seed=14)
    for model in ('baseline','ridge','rf'):
        assert np.array_equal(a[model],b[model])
    assert not np.array_equal(a['observed'],b['observed'])


def test_known_climate_signal_predicts_unseen_observers_and_blocks():
    m=get_module();p,_=m.make_pairs(fixture_frame(),'angle',1,20261101)
    r=m.evaluate(p,'angle',1,seed=20261101,trees=16)
    assert r['status']=='estimated'
    assert r['n_folds_per_direction']==5
    assert r['spatial_overlap_count']==0
    assert r['pooled']['ridge']['skill']>.8
    assert r['pooled']['rf']['skill']>.3
    assert r['n_prediction_cells_both_directions']==2*len(p)
    assert all(x not in json.dumps(r) for x in ['p0_','u0','"taxon"','"latitude"'])


def test_independent_noise_is_not_strongly_predictable():
    m=get_module();p,_=m.make_pairs(fixture_frame(),'angle',1,20261101)
    rng=np.random.default_rng(102)
    for side in ('a','b'):p[m.FEATURES['angle'][0]+'_'+side]=rng.normal(size=len(p))
    r=m.evaluate(p,'angle',1,seed=4,trees=16)
    assert r['pooled']['rf']['skill']<.25
    assert r['pooled']['ridge']['skill']<.25


def test_source_identity_digest_and_roster_join(tmp_path):
    m=get_module();f=fixture_frame();native=tmp_path/'n.csv.gz';roster=tmp_path/'r.csv.gz'
    f.drop(columns='user_id').to_csv(native,index=False);f[['obs_id','user_id']].to_csv(roster,index=False)
    receipt=tmp_path/'receipt.json'
    rec={'prepared_table_sha256':m.sha256(native),'native_observations':len(f),
         'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1',
         'environment_decoding':{e:{'source_sha256':'a'*64} for e in m.PREDICTORS}}
    receipt.write_text(json.dumps(rec))
    d,provenance=m.load_source(native,receipt,roster,expected_n=len(f),expected_users=24)
    assert len(d)==len(f) and d.user_id.notna().all()
    rec['prepared_table_sha256']='0'*64;receipt.write_text(json.dumps(rec))
    with pytest.raises(ValueError,match='AUTHENTICATION'):m.load_source(native,receipt,roster,expected_n=len(f),expected_users=24)


def test_all_eight_endpoint_registry_and_constant_target():
    m=get_module();assert tuple(m.FEATURES)==('angle','lightness','chroma','hue','aspect','circularity','solidity','width_profile')
    p,_=m.make_pairs(fixture_frame(),'angle',1,20261101)
    for side in ('a','b'):p[m.FEATURES['angle'][0]+'_'+side]=1.
    r=m.evaluate(p,'angle',1,seed=2,trees=8)
    assert r['pooled']['rf']['skill'] is None
    assert r['pooled']['ridge']['skill'] is None


def test_duplicate_observations_rejected_before_partition():
    m=get_module();f=fixture_frame()
    with pytest.raises(ValueError,match='DUPLICATE'):
        m.make_pairs(pd.concat([f,f.iloc[[0]]]),'angle',1,20261101)


def test_missing_climate_is_not_imputed():
    m=get_module();f=fixture_frame();f.loc[:,m.PREDICTORS[0]]=np.nan
    p,a=m.make_pairs(f,'angle',1,20261101)
    assert p.empty and a['candidate_records']==0


def test_no_photos_from_same_observer_can_create_support():
    m=get_module();f=fixture_frame();f.loc[:,'user_id']='one_person'
    p,a=m.make_pairs(f,'angle',1,20261101)
    assert p.empty and a['selected_records']==96


def test_all_eight_constructs_are_fit_including_morphology():
    m=get_module();f=fixture_frame()
    for trait in m.FEATURES:
        for degree in (1,5):
            p,a=m.make_pairs(f,trait,degree,20261101)
            r=m.evaluate(p,trait,degree,seed=20261101,trees=8)
            assert r['status']=='estimated'
            assert r['spatial_overlap_count']==0
            m.privacy_check(r)
            assert r['pooled']['ridge']['skill']>0


def test_spatial_fold_without_repeated_taxon_is_not_claimed_predicted():
    m=get_module();p,_=m.make_pairs(fixture_frame(),'angle',1,20261101)
    tr=p[p.taxon=='t0'];te=p[p.taxon=='t1']
    r=m.fit_fold(tr,te,m.FEATURES['angle'],'a','b',False,trees=8,seed=1)
    assert r=={}


def test_privacy_guard_rejects_nested_source_rows():
    m=get_module()
    with pytest.raises(ValueError,match='PRIVATE_OUTPUT'):
        m.privacy_check({'results':[{'user_id':'hidden'}]})


def test_paired_rf_increment_uses_same_baseline_and_taxon_draws():
    m=get_module()
    loss=pd.DataFrame({'_taxon':list('abcdef'),'baseline':[1.]*6,
                       'ridge':[.64]*6,'rf':[.25]*6})
    out=m._score(loss,seed=5)
    assert 'rf_minus_ridge' in out, 'paired RF-minus-Ridge increment missing'
    delta=out['rf_minus_ridge']
    assert delta['skill_difference']==pytest.approx(.39)
    assert delta['taxon_equal_MSE_reduction']==pytest.approx(.39)
    assert delta['taxon_bootstrap_95_conditional_on_fits']==pytest.approx([.39,.39])
    assert out['evaluation_support']['n_evaluated_taxa']==6
    assert out['evaluation_support']['at_least_five_evaluated_taxa']


def test_identical_models_have_zero_paired_increment_and_zero_loss_is_undefined():
    m=get_module()
    loss=pd.DataFrame({'_taxon':list('abcde'),'baseline':[1.,2.,3.,4.,5.],
                      'ridge':[.1,.2,.3,.4,.5],'rf':[.1,.2,.3,.4,.5]})
    out=m._score(loss,seed=2)
    assert 'rf_minus_ridge' in out, 'paired increment missing'
    assert out['rf_minus_ridge']['skill_difference']==pytest.approx(0.)
    assert out['rf_minus_ridge']['taxon_bootstrap_95_conditional_on_fits']==pytest.approx([0.,0.])
    loss.loc[:,['baseline','ridge','rf']]=0.
    out=m._score(loss,seed=2)
    assert out['rf_minus_ridge']['skill_difference'] is None
    assert out['rf_minus_ridge']['taxon_bootstrap_95_conditional_on_fits'] is None


def test_small_scored_taxon_cohort_is_explicit_not_confused_with_paired_cohort():
    m=get_module()
    loss=pd.DataFrame({'_taxon':['only_one']*12,'baseline':[1.]*12,
                      'ridge':[.8]*12,'rf':[.6]*12})
    out=m._score(loss,seed=3)
    assert 'evaluation_support' in out, 'evaluated-taxon support audit missing'
    assert out['evaluation_support']['n_evaluated_taxa']==1
    assert not out['evaluation_support']['at_least_five_evaluated_taxa']

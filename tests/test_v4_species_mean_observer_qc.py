"""Synthetic tests only. No source photographs, locations or observer IDs."""
from __future__ import annotations
import importlib.util
import json
import numpy as np
import pandas as pd
import pytest

MODULE='analysis.v4.species_mean_observer_qc'

def test_species_mean_module_exists():
    assert importlib.util.find_spec(MODULE) is not None, 'species-mean QC has not been implemented'

def mod():
    import importlib
    return importlib.import_module(MODULE)

def frame(taxa=12, users=60, photos=2, signal=4., noise=.1):
    rng=np.random.default_rng(123)
    rows=[]
    # Shared observers across taxa deliberately test globally exclusive assignment.
    for t in range(taxa):
        for u in range(users):
            for p in range(photos):
                x=signal*t+rng.normal(0,noise)
                theta=np.pi-.05+.018*t+rng.normal(0,.002)
                row={'obs_id':f'obs_{t}_{u}_{p}','accepted_key':f'PRIVATE_SPECIES_{t}',
                     'user_id':f'PRIVATE_USER_{u}'}
                for c in ('orientation_image_vertical_angle','corolla_lab_lightness','corolla_lab_chroma',
                          'capitulum_outline_aspect_ratio','capitulum_outline_circularity',
                          'capitulum_outline_solidity','capitulum_width_profile_cv'):
                    row[c]=x;row[c+'__eligible']=True
                for c,v in [('corolla_hue_sin',np.sin(theta)),('corolla_hue_cos',np.cos(theta))]:
                    row[c]=v;row[c+'__eligible']=True
                rows.append(row)
    return pd.DataFrame(rows)

def test_global_observers_and_high_species_mean_agreement():
    m=mod();d=frame();r=m.audit_construct(d,'angle',replicates=30)
    assert r['n_common_taxa_all_splits_and_methods']==12
    assert len(r['splits'])==6
    assert {s['method'] for s in r['splits']}=={'photo_mean','observer_equal_mean'}
    for s in r['splits']:
        assert s['metrics']['pearson_r']>.999
        assert s['metrics']['lin_ccc']>.999
        assert s['global_observer_overlap']==0
        assert s['support']['n_taxa']==12

def test_perfect_correlation_does_not_hide_absolute_bias():
    m=mod();a=np.arange(10.)[:,None];b=a+20
    r=m.agreement(a,b,replicates=40,seed=3)
    assert r['metrics']['pearson_r']==pytest.approx(1)
    assert r['metrics']['lin_ccc']<.05
    assert r['metrics']['bias_b_minus_a']==pytest.approx(20)
    assert r['metrics']['rmse']==pytest.approx(20)
    assert r['bootstrap_95']['rmse']['lower']==pytest.approx(20)

def test_scale_shift_detected_by_concordance():
    m=mod();a=np.arange(1.,11.)[:,None]
    r=m.agreement(a,2*a,replicates=20,seed=1)['metrics']
    assert r['pearson_r']==pytest.approx(1)
    assert r['lin_ccc']<.8

def test_constant_species_values_not_called_reliable():
    r=mod().agreement(np.ones((9,1)),np.ones((9,1)),replicates=20,seed=2)
    assert r['metrics']['pearson_r'] is None
    assert r['metrics']['lin_ccc'] is None
    assert r['metrics']['rmse']==0
    assert r['bootstrap_95']['pearson_r']['status']=='insufficient_finite_replicates'

def test_hue_wraparound_error_is_two_degrees_not_358():
    a=np.deg2rad(np.full(10,179.));b=np.deg2rad(np.full(10,-179.))
    av=np.c_[np.sin(a),np.cos(a)];bv=np.c_[np.sin(b),np.cos(b)]
    r=mod().agreement(av,bv,circular=True,replicates=20,seed=3)['metrics']
    assert r['angular_rmse_degrees']==pytest.approx(2)
    assert 'lin_ccc' not in r and 'spearman_rho' not in r
    assert r['mean_cosine_agreement']>.999

def test_row_order_does_not_change_results():
    m=mod();d=frame(taxa=8)
    assert m.audit_construct(d,'angle',replicates=20)==m.audit_construct(d.sample(frac=1,random_state=2),'angle',replicates=20)

def test_missing_quality_flag_is_not_assumed_true():
    d=frame(taxa=6).drop(columns='orientation_image_vertical_angle__eligible')
    with pytest.raises(ValueError,match='MISSING_REQUIRED_COLUMNS'):mod().audit_construct(d,'angle',replicates=10)

def test_duplicate_observation_is_rejected():
    d=frame(taxa=6);d=pd.concat([d,d.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match='OBSERVATION_IDENTITY'):mod().audit_construct(d,'angle',replicates=10)

def test_missing_observer_is_rejected():
    d=frame(taxa=6);d.loc[0,'user_id']=None
    with pytest.raises(ValueError,match='OBSERVER_IDENTITY'):mod().audit_construct(d,'angle',replicates=10)

def test_prolific_observer_weight_is_a_sensitivity_not_silent_replacement():
    m=mod();d=frame(taxa=6,noise=0)
    target='PRIVATE_USER_0';mask=d.user_id==target
    extra=pd.concat([d.loc[mask].assign(obs_id=lambda v:v.obs_id+'_extra'+str(i)) for i in range(15)],ignore_index=True)
    extra['orientation_image_vertical_angle']+=80
    d=pd.concat([d,extra],ignore_index=True)
    pairs=m.species_pairs(d,'angle',m.SEEDS[0])
    assert set(pairs)=={'photo_mean','observer_equal_mean'}
    p=pairs['photo_mean'];q=pairs['observer_equal_mean']
    assert p.index.equals(q.index)
    assert not np.allclose(p[['value_a','value_b']],q[['value_a','value_b']])

def test_bad_qc_is_excluded_without_environment_or_geography_gates():
    m=mod();d=frame(taxa=6)
    d.loc[d.obs_id.str.endswith('_0'),'orientation_image_vertical_angle__eligible']=False
    r=m.audit_construct(d,'angle',replicates=20)
    assert r['n_trait_qc_observations']==len(d)//2
    assert r['n_common_taxa_all_splits_and_methods']==6

def test_supported_taxa_are_intersected_across_seeds_and_weightings():
    m=mod();d=frame(taxa=8)
    # One taxon has only one user in one side of seed1; it must not enter ANY seed.
    one=[u for u in d.user_id.unique() if m.side(u,m.SEEDS[0])==1][0]
    d=d.loc[(d.accepted_key!='PRIVATE_SPECIES_0') | (d.user_id.map(lambda x:m.side(x,m.SEEDS[0]))==0) | (d.user_id==one)]
    r=m.audit_construct(d,'angle',replicates=20)
    assert r['n_common_taxa_all_splits_and_methods']==7
    assert all(s['support']['n_taxa']==7 for s in r['splits'])

def test_insufficient_taxa_reported_not_fitted():
    r=mod().audit_construct(frame(taxa=4),'angle',replicates=10)
    assert r['status']=='insufficient_supported_taxa'
    assert r['splits']==[]

def test_all_eight_traits_are_preserved_and_public_result_has_no_ids():
    m=mod();d=frame(taxa=6)
    for c in m.FEATURES:
        r=m.audit_construct(d,c,replicates=10)
        assert len(r['splits'])==6
        assert 'PRIVATE_' not in json.dumps(r)
        m.privacy_check(r)
    assert len(m.FEATURES)==8

def test_source_hashes_counts_and_roster_join(tmp_path):
    m=mod();d=frame(taxa=6);f=tmp_path/'native.csv';u=tmp_path/'roster.csv';rp=tmp_path/'receipt.json'
    d.drop(columns='user_id').to_csv(f,index=False);d[['obs_id','user_id']].to_csv(u,index=False)
    rp.write_text(json.dumps({'prepared_table_sha256':m.sha256(f),'native_observations':len(d),
                             'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    opts=dict(expected_n=len(d),expected_users=60,source_sha=m.sha256(f),roster_sha=m.sha256(u))
    loaded,_=m.load_source(f,rp,u,**opts)
    assert len(loaded)==len(d)
    with pytest.raises(ValueError,match='SOURCE_HASH'):m.load_source(f,rp,u,**{**opts,'source_sha':'0'*64})
    with pytest.raises(ValueError,match='ROSTER_HASH'):m.load_source(f,rp,u,**{**opts,'roster_sha':'0'*64})
    bad=d[['obs_id','user_id']].copy();bad.loc[0,'obs_id']='foreign';bad.to_csv(u,index=False)
    with pytest.raises(ValueError,match='SOURCE_OR_ROSTER_IDENTITY'):m.load_source(f,rp,u,**{**opts,'roster_sha':m.sha256(u)})

def test_production_source_expectations_are_pinned_to_existing_receipts():
    m=mod()
    assert m.EXPECTED_N==313350 and m.EXPECTED_USERS==92839
    assert m.SOURCE_SHA=='342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed'
    assert m.ROSTER_SHA=='a89a586d7aaff7ab4b6c936aa4bfc30c597c152a07bdb49ec7abddb7fee04fef'

def test_invalid_construct_and_public_rows_rejected():
    m=mod()
    with pytest.raises(ValueError,match='UNKNOWN_CONSTRUCT'):m.audit_construct(frame(taxa=5),'invented',replicates=5)
    with pytest.raises(ValueError,match='PRIVATE_OUTPUT_FIELD'):m.privacy_check({'nested':[{'taxon':'x'}]})

def test_regenerated_roster_requires_explicit_upstream_receipt(tmp_path):
    m=mod();d=frame(taxa=6);f=tmp_path/'native.csv';u=tmp_path/'roster.csv.gz';rp=tmp_path/'receipt.json';ur=tmp_path/'roster_receipt.json'
    d.drop(columns='user_id').to_csv(f,index=False)
    d[['obs_id','user_id']].to_csv(u,index=False,compression={'method':'gzip','mtime':123})
    rp.write_text(json.dumps({'prepared_table_sha256':m.sha256(f),'native_observations':len(d),
        'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    ur.write_text(json.dumps({'schema':'azami_authenticated_observer_roster_receipt_v1',
        'native_decoded_sha256':m.sha256(f),'roster_table_sha256':m.sha256(u),
        'source_encrypted_sha256':'dfa98fa162a8029cb6477adf4c26ffc5819cc916a0513d32da14d78d845ab4b0',
        'source_plaintext_zip_sha256':'2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae',
        'n_observations':len(d),'n_observers':60,'authentication_status':'original_source_verified_upstream'}))
    opts=dict(expected_n=len(d),expected_users=60,source_sha=m.sha256(f))
    with pytest.raises(ValueError,match='ROSTER_HASH'):m.load_source(f,rp,u,**opts)
    got,info=m.load_source(f,rp,u,roster_receipt=ur,**opts)
    assert len(got)==len(d)
    assert info['roster_source_validation']=='upstream_attestation_checked_not_independently_reauthenticated'
    bad=json.loads(ur.read_text());bad['source_plaintext_zip_sha256']='0'*64;ur.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match='ROSTER_RECEIPT'):m.load_source(f,rp,u,roster_receipt=ur,**opts)

def test_mean_definitions_match_independent_manual_grouping():
    m=mod();d=frame(taxa=6,users=60,photos=3)
    c='orientation_image_vertical_angle';p=m.species_pairs(d,'angle',m.SEEDS[0])
    tax='PRIVATE_SPECIES_0'
    for half,suffix in [(0,'a'),(1,'b')]:
        sub=d.loc[(d.accepted_key==tax) & d.user_id.map(lambda v:m.side(v,m.SEEDS[0])==half)]
        assert p['photo_mean'].loc[tax,'value_'+suffix]==pytest.approx(sub[c].sum()/len(sub))
        by_user=sub.groupby('user_id')[c].mean()
        assert p['observer_equal_mean'].loc[tax,'value_'+suffix]==pytest.approx(by_user.sum()/len(by_user))

def test_minimum_observer_support_is_not_replaced_by_many_photos():
    m=mod();d=frame(taxa=6,users=2,photos=100)
    r=m.audit_construct(d,'angle',replicates=10)
    assert r['status']=='insufficient_supported_taxa'
    assert r['n_common_taxa_all_splits_and_methods']==0

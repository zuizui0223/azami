"""Observer-disjoint original-photo S1 tests; only synthetic fake user IDs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_photographic_observer_disjoint_reliability import (
    _observer_side, _observer_cell_pairs, audit_observer_half, run,sha256,
)
from analysis.v4.rf_photographic_split_half_reliability import _spatial_photo_cohort


def synthetic(*,taxa=25,cell=10,nphoto=12,shared_observer=False,biological_effect=8.,noise=.2):
    rng=np.random.default_rng(88)
    values=[]
    for t in range(taxa):
        for k in range(cell):
            for p in range(nphoto):
                user=f'author_{t}_{k}_{p}' if not shared_observer else f'author_{t}_{k}_{p//(nphoto//2)}'
                q=biological_effect*(k-4.5)+rng.normal(0,noise)
                theta=3.11+(k-4.5)*.015+rng.normal(0,noise*.01)
                entry={'obs_id':f'originalphoto_{t}_{k}_{p}', 'user_id':user,
                       'accepted_key':f'taxon_{t}',
                       'latitude':-43+6*k+p*.00001,'longitude':-110+8*(t%3)+p*.00001,
                       'orientation_image_vertical_angle':q,
                       'corolla_lab_lightness':q,'corolla_lab_chroma':q,
                       'corolla_hue_sin':np.sin(theta),'corolla_hue_cos':np.cos(theta),
                       'capitulum_outline_aspect_ratio':q,'capitulum_outline_circularity':q,
                       'capitulum_outline_solidity':q,'capitulum_width_profile_cv':q}
                for name in [c for cs in FEATURES.values() for c in cs]:
                    entry[name+'__eligible']=True
                values.append(entry)
    return pd.DataFrame(values)


def observer_map(df):
    return df[['obs_id','user_id']].copy()


def test_global_observer_side_is_stable_and_exclusive_across_cells():
    ids=pd.Series(['observer_a','observer_b','observer_a','observer_b'])
    side=_observer_side(ids,0)
    assert side[0]==side[2] and side[1]==side[3]
    assert np.all(np.isin(side,[0,1]))
    assert _observer_side(ids,0).tolist()==side.tolist()


def test_estimated_source_signal_with_many_disjoint_observers():
    d=synthetic()
    summary=audit_observer_half(d,observer_map(d),'angle',degrees=1)
    assert summary['status']=='estimated_observer_disjoint_spatial_photo_QC_not_RF'
    assert summary['n_successful_observer_disjoint_splits']==8
    assert summary['cross_observer_half_skill_vs_train_observer_half_taxon_median']['q50']>.8
    assert summary['cross_observer_half_taxon_centered_cell_correlation']['q50']>.9
    assert summary['n_splits_positive_skill']==8
    assert 'author_0' not in json.dumps(summary) and 'taxon_0' not in json.dumps(summary)


def test_observer_imbalance_removes_cells_and_cannot_silent_relax():
    d=synthetic(shared_observer=True)
    summary=audit_observer_half(d,observer_map(d),'angle',degrees=1)
    assert summary['status'] in ('insufficient_geographical_observer_independent_replication',
                                'estimated_observer_disjoint_spatial_photo_QC_not_RF')
    assert summary['n_cells_with_two_observer_halves_by_split']['q50']<250
    if summary['status']=='insufficient_geographical_observer_independent_replication':
        assert summary['no_skill_claim_from_insufficient_support']


def test_different_user_id_mappings_refuse_missing_or_duplicate_join():
    d=synthetic(taxa=5)
    photo,_=_spatial_photo_cohort(d,'angle',1)
    mapping=observer_map(d)
    missing=mapping.iloc[1:].copy()
    with pytest.raises(ValueError,match='JOIN_LOST'):
        _observer_cell_pairs(photo,missing,'angle',rep=0)
    duplicate=pd.concat([mapping,mapping.iloc[:1]],ignore_index=True)
    with pytest.raises(ValueError,match='DUPLICATE'):
        _observer_cell_pairs(photo,duplicate,'angle',rep=0)


def test_hue_joint_circular_works_for_observer_splits():
    d=synthetic()
    out=audit_observer_half(d,observer_map(d),'hue',degrees=5)
    assert out['status']=='estimated_observer_disjoint_spatial_photo_QC_not_RF'
    assert out['cross_observer_half_taxon_centered_cell_correlation']['q50']>.8


def test_no_signal_does_not_improve_observer_disjoint_test():
    d=synthetic(biological_effect=0.,noise=4.)
    out=audit_observer_half(d,observer_map(d),'chroma',degrees=1)
    assert out['status']=='estimated_observer_disjoint_spatial_photo_QC_not_RF'
    assert out['cross_observer_half_skill_vs_train_observer_half_taxon_median']['q50']<0


def test_authenticated_source_gate_and_private_aggregate(monkeypatch,tmp_path):
    import analysis.v4.rf_photographic_observer_disjoint_reliability as mod
    df=synthetic()
    raw=tmp_path/'raw.csv.gz';df.drop(columns=['user_id']).to_csv(raw,index=False,compression='gzip')
    receipt=tmp_path/'receipt.json'
    receipt.write_text(json.dumps({'prepared_table_sha256':sha256(raw),
                                   'native_observations':len(df),
                                   'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    sealed=tmp_path/'sealed.aesgcm';sealed.write_bytes(b'not a real secret')
    monkeypatch.setenv('AZAMI_METADATA_KEY','a'*64)
    monkeypatch.setattr(mod,'EXPECTED_NATIVE_OBSERVERS',df.user_id.nunique())
    monkeypatch.setattr(mod,'original_observer_roster',lambda sealed,key,ids:observer_map(df))
    target=tmp_path/'aggregate.json'
    out=run(raw,receipt,sealed,target,expected_n=len(df),
            constructs=('angle','chroma'),grains=(1,5))
    assert out['n_original_distinct_observers']==df.user_id.nunique()
    assert len(out['results'])==4
    assert out['original_frozen_GEB_S1_S2_H3c_and_pr167_RF_untouched']
    for forbidden in ('originalphoto_','taxon_0','author_0','"obs_id"','"accepted_key"',
                      '"user_id"','"latitude"','"longitude"','"observer_id"'):
        assert forbidden not in target.read_text()
    receipt.write_text(json.dumps({'prepared_table_sha256':'0'*64,
                                   'native_observations':len(df),
                                   'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    with pytest.raises(ValueError,match='NOT_AUTHENTICATED'):
        run(raw,receipt,sealed,target,expected_n=len(df),constructs=('angle',))


def test_invalid_design_scope():
    d=synthetic(taxa=5)
    with pytest.raises(ValueError,match='UNREGISTERED'):
        audit_observer_half(d,observer_map(d),'hue',degrees=2)
    with pytest.raises(ValueError,match='UNREGISTERED'):
        audit_observer_half(d,observer_map(d),'angle',degrees=1,reps=1)
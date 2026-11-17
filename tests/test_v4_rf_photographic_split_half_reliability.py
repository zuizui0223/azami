"""Synthetic-only original photo QC. No private S1 source, no raster network."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_photographic_split_half_reliability import (
    _spatial_photo_cohort, _paired_half_cells, audit_split_half, run, sha256,
)


def sample_cohort(*, taxa=25, cells=10, photos=6, geo_strength=4., noise=.2):
    rng=np.random.default_rng(22)
    data=[]
    for t in range(taxa):
        for k in range(cells):
            for p in range(photos):
                val=geo_strength*(k-(cells-1)/2)+rng.normal(0,noise)
                theta=np.pi-.11+k*.035 + rng.normal(0,noise*.01)
                row={"obs_id":f"privatephoto_{t}_{k}_{p}","accepted_key":f"Taxon_{t}",
                     "latitude":-43+6*k+p*.0001,"longitude":-110+5*(t%3)+p*.0001,
                     "orientation_image_vertical_angle":val,
                     "corolla_lab_lightness":val,
                     "corolla_lab_chroma":val,
                     "corolla_hue_sin":float(np.sin(theta)),
                     "corolla_hue_cos":float(np.cos(theta)),
                     "capitulum_outline_aspect_ratio":val,
                     "capitulum_outline_circularity":val,
                     "capitulum_outline_solidity":val,
                     "capitulum_width_profile_cv":val}
                for n in (x for v in FEATURES.values() for x in v):
                    row[n+"__eligible"]=True
                data.append(row)
    return pd.DataFrame.from_records(data)


def test_hash_split_reproducible_disjoint_balanced_three_photo_per_half():
    df=sample_cohort()
    photos,first =_spatial_photo_cohort(df,'angle',1)
    assert first['n_qualified_taxa']==25
    assert first['n_qualified_cells']==250
    assert first['n_original_photos_in_supported_cells']==1500
    a,st=_paired_half_cells(photos,'angle',rep=0)
    b,st2=_paired_half_cells(photos,'angle',rep=0)
    assert a.equals(b) and st==st2
    assert len(a)==250
    assert st['no_same_observation_in_both_halves']


def test_geographic_signal_reproduces_across_disjoint_photo_halves():
    df=sample_cohort()
    res=audit_split_half(df,'angle',degrees=1)
    assert res['status']=='estimated_source_split_half_QC_not_RF'
    assert res['n_successful_predeclared_splits']==8
    assert res['n_splits_with_positive_spatial_cell_skill']==8
    assert res['within_taxon_cross_half_cell_agreement_correlation']['q50']>.95
    assert res['two_way_cross_half_skill_vs_train_taxon_median']['q50']>.95
    assert res['gain_over_shuffled_cell_correspondence']['q50']>0.5
    assert 'Taxon_0' not in json.dumps(res) and 'privatephoto_' not in json.dumps(res)


def test_non_spatial_noise_has_no_claimed_geographic_reliability():
    df=sample_cohort(geo_strength=0.,noise=3.)
    res=audit_split_half(df,'angle',degrees=1)
    assert res['status']=='estimated_source_split_half_QC_not_RF'
    assert res['two_way_cross_half_skill_vs_train_taxon_median']['q50']<0
    assert res['n_splits_with_positive_spatial_cell_skill']<8


def test_joint_circular_hue_across_pi_wrap_boundary():
    df=sample_cohort(geo_strength=0.,noise=.15)
    res=audit_split_half(df,'hue',degrees=5)
    assert res['status']=='estimated_source_split_half_QC_not_RF'
    assert res['support']['n_qualified_taxa']==25
    assert res['within_taxon_cross_half_cell_agreement_correlation']['q50']>.9
    assert res['two_way_cross_half_skill_vs_train_taxon_median']['q50']>0


def test_under_four_photos_per_cell_fail_closed_without_reusing_photo():
    df=sample_cohort(photos=3)
    r=audit_split_half(df,'angle',degrees=1)
    assert r['status']=='insufficient_repeated_taxa_for_split_half'
    assert r['support']['n_qualified_cells']==0


def test_eligibility_fails_closed_and_detects_duplicate_ids():
    df=sample_cohort()
    df.loc[df.obs_id.str.endswith('_0'),'orientation_image_vertical_angle__eligible']=False
    r=audit_split_half(df,'angle',degrees=1)
    assert r['status']=='estimated_source_split_half_QC_not_RF'
    assert r['support']['photos_per_supported_cell']['q50']==5
    wrong=df.copy()
    wrong.loc[0,'obs_id']=wrong.loc[1,'obs_id']
    with pytest.raises(ValueError,match='DUPLICATE'):
        audit_split_half(wrong,'angle',degrees=1)


def test_source_hash_grain_and_aggregate_no_private_keys(tmp_path):
    df=sample_cohort()
    f=tmp_path/'original.csv.gz';df.to_csv(f,index=False,compression='gzip')
    receipt=tmp_path/'receipt.json'
    receipt.write_text(json.dumps({'prepared_table_sha256':sha256(f),
                                  'native_observations':len(df),
                                  'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    out=tmp_path/'aggregate.json'
    d=run(f,receipt,out,expected_n=len(df),constructs=('angle','hue'),grains=(1,5))
    assert len(d['results'])==4
    assert all(x['status']=='estimated_source_split_half_QC_not_RF' for x in d['results'])
    assert d['original_GEB_S1_S2_Holm_H3c_and_RF_untouched']
    content=out.read_text()
    for value in ('privatephoto','Taxon_0','"obs_id"','"accepted_key"','"latitude"','"longitude"','"user_id"','"photographer_id"'):
        assert value not in content
    receipt.write_text(json.dumps({'prepared_table_sha256':'0'*64,'native_observations':len(df),
                                  'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    with pytest.raises(ValueError,match='NOT_AUTHENTICATED'):
        run(f,receipt,out,expected_n=len(df),constructs=('angle',),grains=(1,))


def test_invalid_scope_and_repeated_seed_control():
    d=sample_cohort(taxa=8)
    with pytest.raises(ValueError,match='FIXED_EIGHT'):
        audit_split_half(d,'angle',degrees=1,reps=2)
    with pytest.raises(ValueError,match='FIXED_EIGHT'):
        audit_split_half(d,'angle',degrees=2)
    with pytest.raises(ValueError,match='FIXED_EIGHT'):
        audit_split_half(d,'fictional',degrees=1)
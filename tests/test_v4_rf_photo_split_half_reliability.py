"""Synthetic checks only; original private S1 must never be downloaded here."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.rf_photo_split_half_reliability import (
    _hash_obs_ids, audit_construct, run, sha256,
)


def cohort(*, taxa=12, cell_effect=3., noise=.2, photos=8, hue=False):
    rng=np.random.default_rng(20261101)
    rows=[]
    for t in range(taxa):
        for cell in range(10):
            for photo in range(photos):
                x=cell_effect*(cell-4.5)+rng.normal(0.,noise)
                d={'obs_id':f'private_{t}_{cell}_{photo}',
                   'accepted_key':f'secret_taxon_{t}',
                   'latitude': -62.+6.*cell+.01*photo,
                   'longitude':-130.+3.*t+.01*photo,
                   'orientation_image_vertical_angle':x,
                   'orientation_image_vertical_angle__eligible':True}
                if hue:
                    theta=np.pi-.12+0.06*(cell-4.5)*cell_effect+rng.normal(0.,noise/10.)
                    d.update(corolla_hue_sin=np.sin(theta),corolla_hue_cos=np.cos(theta),
                             corolla_hue_sin__eligible=True,corolla_hue_cos__eligible=True)
                rows.append(d)
    return pd.DataFrame(rows)


def test_hashed_balanced_halves_recover_known_geographical_pattern():
    df=_hash_obs_ids(cohort(),(20261101,20261103,20261105))
    res=audit_construct(df,'angle',degrees=1)
    assert res['n_original_photo_trait_QC_eligible']==960
    for s in res['splits']:
        assert s['status']=='estimated_exploratory_reliability'
        assert s['n_taxa']==12 and s['n_qualified_cells']==120
        assert s['n_contributing_photos']==960
        assert s['taxon_equal_spatial_alignment']['q50']>.95
        assert s['median_alignment_taxon_bootstrap_95'][0]>.8
        assert s['taxon_equal_spearman']['q50']>.9
        assert s['within_taxon_cell_labels_permuted_null_median_alignment']['q95']<.5
        assert s['photos_per_half_cell']['q50']==4.


def test_no_geographical_pattern_does_not_exhibit_strong_reproducibility():
    df=_hash_obs_ids(cohort(cell_effect=0.,noise=20.,photos=8,taxa=25), (20261101,))
    out=audit_construct(df,'angle',degrees=5,seeds=(20261101,))['splits'][0]
    assert out['status']=='estimated_exploratory_reliability'
    assert out['taxon_equal_spatial_alignment']['q50'] < .5
    assert out['half_difference_over_among_cell_RMSE']['q50'] > .7


def test_hue_wraparound_reproducibility_uses_joint_vectors():
    df=_hash_obs_ids(cohort(hue=True,cell_effect=1.,noise=.1), (20261101,))
    r=audit_construct(df,'hue',degrees=1,seeds=(20261101,))['splits'][0]
    assert r['status']=='estimated_exploratory_reliability'
    assert r['taxon_equal_spatial_alignment']['q50']>.8
    assert r['taxon_equal_spearman'] is None


def test_minimum_per_half_excludes_under_supported_cells():
    df=_hash_obs_ids(cohort(photos=5), (20261101,))
    r=audit_construct(df,'angle',degrees=1,seeds=(20261101,))['splits'][0]
    assert r['status']=='insufficient_independently_supported_cells'
    assert r['n_taxa']==0


def test_source_receipt_check_unique_id_gate_and_no_private_ids_in_output(tmp_path):
    df=cohort(hue=True)
    file=tmp_path/'private.csv.gz';df.to_csv(file,index=False,compression='gzip')
    receipt=tmp_path/'receipt.json'; receipt.write_text(json.dumps({
        'prepared_table_sha256':sha256(file),'native_observations':len(df),
        'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    outpath=tmp_path/'aggregate.json'
    out=run(file,receipt,outpath,expected_n=len(df),constructs=('angle','hue'),
            grains=(1,),seeds=(20261101,))
    assert len(out['results'])==2
    assert out['source_taxon_cell_photo_location_data_not_exported']
    txt=outpath.read_text()
    for private in ['private_0_0_0','secret_taxon_0','"obs_id"',
                    '"accepted_key"', '"latitude"', '"longitude"', '"taxon"']:
        assert private not in txt
    receipt.write_text(json.dumps({'prepared_table_sha256':'0'*64,
        'native_observations':len(df),
        'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    with pytest.raises(ValueError,match='SHA_NOT_AUTHENTICATED'):
        run(file,receipt,outpath,expected_n=len(df),constructs=('angle',),
            grains=(1,),seeds=(20261101,))
    receipt.write_text(json.dumps({'prepared_table_sha256':sha256(file),
        'native_observations':len(df),
        'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    df=pd.concat([df,df.iloc[[0]]],ignore_index=True)
    df.to_csv(file,index=False,compression='gzip')
    receipt.write_text(json.dumps({'prepared_table_sha256':sha256(file),
        'native_observations':len(df),
        'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    with pytest.raises(ValueError,match='DUPLICATE_OR_MISSING_ORIGINAL_OBS_ID'):
        run(file,receipt,outpath,expected_n=len(df),constructs=('angle',),
            grains=(1,),seeds=(20261101,))


def test_invalid_endpoint_and_grain_fail_closed():
    df=_hash_obs_ids(cohort(), (20261101,))
    with pytest.raises(ValueError,match='INVALID_REGISTERED'):
        audit_construct(df,'invented',degrees=1,seeds=(20261101,))
    with pytest.raises(ValueError,match='INVALID_REGISTERED'):
        audit_construct(df,'angle',degrees=2,seeds=(20261101,))

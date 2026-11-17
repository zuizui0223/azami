"""Synthetic-only authenticated global-observer-disjoint split checks."""
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.rf_photo_observer_disjoint_reliability import (
    hash_side, audit, run, sha256,
)


def cohort(*, taxa=12, observers=24, cell_effect=1., with_hue=False):
    rng=np.random.default_rng(20261101)
    records=[]
    for t in range(taxa):
        camera=rng.normal(0,.3,size=observers)
        for cell in range(10):
            for u in range(observers):
                val=cell_effect*(cell-4.5)+camera[u]+rng.normal(0,.4)
                row={'obs_id':f'private_{t}_{cell}_{u}',
                     'user_id':f'secret_user_{t}_{u}',
                     'accepted_key':f'secret_taxon_{t}',
                     'latitude':-60.+6.*cell+.001*u,
                     'longitude':-120.+2.*t+.001*u,
                     'orientation_image_vertical_angle':val,
                     'orientation_image_vertical_angle__eligible':True}
                if with_hue:
                    theta=np.pi-.1+.035*cell_effect*(cell-4.5)+rng.normal(0,.02)
                    row.update(corolla_hue_sin=np.sin(theta),corolla_hue_cos=np.cos(theta),
                               corolla_hue_sin__eligible=True,corolla_hue_cos__eligible=True)
                records.append(row)
    return pd.DataFrame(records)


def test_global_user_partition_never_crosses_halves():
    ids=['a','a','b','b','c']
    for seed in (20261101,20261103,20261105):
        assert hash_side(ids[0],seed)==hash_side(ids[1],seed)
        assert hash_side(ids[2],seed)==hash_side(ids[3],seed)


def test_true_geographic_signal_repeats_between_unshared_observers():
    df=cohort(with_hue=True)
    for construct in ('angle','hue'):
        r=audit(df,construct,degrees=1)
        for s in r['splits']:
            assert s['status']=='estimated_exploratory_reliability'
            assert s['n_taxa']>=10
            assert s['taxon_equal_spatial_alignment']['q50']>.85
            assert s['median_alignment_taxon_bootstrap_95'][0]>.6


def test_null_geographic_signal_cannot_fake_large_alignment():
    df=cohort(taxa=24,cell_effect=0.)
    r=audit(df,'angle',degrees=1,seeds=(20261101,))['splits'][0]
    assert r['status']=='estimated_exploratory_reliability'
    assert r['taxon_equal_spatial_alignment']['q50']<.5


def test_single_observer_per_cell_is_excluded_not_leaked():
    df=cohort(observers=1)
    r=audit(df,'angle',degrees=1,seeds=(20261101,))['splits'][0]
    assert r['status']=='insufficient_independently_supported_cells'


def test_authenticated_mapping_and_no_private_export(tmp_path):
    df=cohort(with_hue=True)
    src=tmp_path/'source.csv.gz'; df.drop(columns=['user_id']).to_csv(src,index=False,compression='gzip')
    roster=tmp_path/'roster.csv.gz'
    df[['obs_id','user_id']].to_csv(roster,index=False,compression='gzip')
    receipt=tmp_path/'receipt.json'; receipt.write_text(json.dumps({
        'prepared_table_sha256':sha256(src),
        'native_observations':len(df),
        'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1'}))
    out=tmp_path/'out.json'
    report=run(src,receipt,roster,out,expected_n=len(df),
               expected_observers=df.user_id.nunique(),
               constructs=('angle','hue'),grains=(1,),seeds=(20261101,))
    assert len(report['results'])==2
    txt=out.read_text()
    for secret in ('private_0','secret_user','secret_taxon','"obs_id"',
                   '"user_id"','"latitude"','"longitude"'):
        assert secret not in txt
    modified=df[['obs_id','user_id']].copy();modified.loc[0,'user_id']=None
    modified.to_csv(roster,index=False,compression='gzip')
    with pytest.raises(ValueError,match='OBSERVER_ROSTER_COUNT_OR_IDENTITY_DRIFT'):
        run(src,receipt,roster,out,expected_n=len(df),
            expected_observers=df.user_id.nunique(),constructs=('angle',),
            grains=(1,),seeds=(20261101,))

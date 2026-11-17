"""Synthetic-original-only photo dispersion checks; no private S1 in this suite."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import pytest

from analysis.v4.rf_photographic_within_cell_dispersion import (
    audit_construct, run, sha256,
)


def fake_cohort(*, taxa=25, cell_effect=3., within=1., hue=False):
    rng=np.random.default_rng(111)
    records=[]
    for taxon in range(taxa):
        for cell in range(9):
            for photo in range(5):
                x=cell_effect*(cell-4)+rng.normal(0,within)
                r={'obs_id':f'p_{taxon}_{cell}_{photo}',
                   'accepted_key':f't_{taxon}',
                   'latitude':-46.+7.2*cell+.00001*photo,
                   'longitude':-131.+1.1*taxon+.00001*photo,
                   'orientation_image_vertical_angle':x,
                   'orientation_image_vertical_angle__eligible':True}
                if hue:
                    theta=np.pi-.17+(cell-4)*.045+rng.normal(0,.025)
                    r.update(corolla_hue_sin=np.sin(theta),corolla_hue_cos=np.cos(theta),
                             corolla_hue_sin__eligible=True,
                             corolla_hue_cos__eligible=True)
                records.append(r)
    return pd.DataFrame.from_records(records)


def test_between_cell_spread_exceeds_within_when_spatial_effect_is_large():
    df=fake_cohort()
    for grain in (1,5):
        result=audit_construct(df,'angle',degrees=grain)
        assert result['status']=='estimated_QC_only'
        assert result['n_qualified_taxa']==25
        assert result['n_qualified_cells']==225
        assert result['n_photo_records_in_supported_cells']==1125
        assert result['among_over_within_RMSE_per_taxon']['q50']>3.
        assert result['n_taxa_among_spread_below_within_spread']==0
        assert 't_0' not in json.dumps(result)
        assert 'p_0' not in json.dumps(result)


def test_photo_dispersion_dominates_when_local_cell_medians_are_flat():
    d=audit_construct(fake_cohort(cell_effect=0.,within=5.),'angle',degrees=1)
    assert d['n_qualified_taxa']==25
    assert d['among_over_within_RMSE_per_taxon']['q50']<1.
    assert d['n_taxa_among_spread_below_within_spread']>=20


def test_circular_hue_uses_normalised_chord_distance_not_wrapped_angle():
    r=audit_construct(fake_cohort(hue=True), 'hue', degrees=5)
    assert r['status']=='estimated_QC_only'
    assert r['distance_unit']=='circular_chord_on_unit_sin_cos'
    assert r['n_qualified_taxa']==25
    assert 0<r['taxon_equal_within_cell_RMSE']['q50']<0.12
    assert r['among_over_within_RMSE_per_taxon']['q50']>1.


def test_ineligible_phenotypes_do_not_enter_photo_summary():
    df=fake_cohort()
    mask=df['obs_id'].str.endswith('_0')
    df.loc[mask,'orientation_image_vertical_angle__eligible']=False
    r=audit_construct(df,'angle',degrees=1,min_photos_taxon=35)
    assert r['n_eligible_photo_QC_before_cell_gate']==900
    assert r['n_photo_records_in_supported_cells']==900
    assert r['photos_per_supported_taxon_cell']['q50']==4.


def test_insufficient_geographical_taxon_repeat_is_fail_closed():
    r=audit_construct(fake_cohort(taxa=3),'angle',degrees=1,min_cells=12)
    assert r['status']=='insufficient_repeated_taxon_cells'
    assert r['n_qualified_taxa']==0


def test_authenticated_raw_archive_gate_and_aggregate_privacy(tmp_path):
    df=fake_cohort()
    file=tmp_path/'source.csv.gz';df.to_csv(file,index=False,compression='gzip')
    rec=tmp_path/'receipt.json'
    rec.write_text(json.dumps({'prepared_table_sha256':sha256(file),
                               'native_observations':len(df),
                               'environment_decoding_status':
                                   'physical_units_gsp_invalid_encoding_masked_v1'}))
    output=tmp_path/'aggregate.json'
    out=run(file,rec,output,expected_n=len(df),constructs=('angle',),grains=(1,5))
    assert out['source_original_s1_observations']==len(df)
    assert len(out['results'])==2
    assert out['not_a_new_Holm_p_value_family_or_frozen_GEB_result']
    public=output.read_text()
    for secret in ('p_0_0_0','t_0','"obs_id"','"latitude"','"longitude"','"accepted_key"'):
        assert secret not in public
    rec.write_text(json.dumps({'prepared_table_sha256':'0'*64,
                               'native_observations':len(df),
                               'environment_decoding_status':
                                    'physical_units_gsp_invalid_encoding_masked_v1'}))
    with pytest.raises(ValueError,match='SHA_NOT_AUTHENTICATED'):
        run(file,rec,output,expected_n=len(df),constructs=('angle',),grains=(1,))


def test_refuse_unknown_endpoint_or_grid_grain():
    df=fake_cohort()
    with pytest.raises(ValueError,match='INVALID_REGISTERED'):
        audit_construct(df,'none',degrees=1)
    with pytest.raises(ValueError,match='INVALID_REGISTERED'):
        audit_construct(df,'angle',degrees=2)
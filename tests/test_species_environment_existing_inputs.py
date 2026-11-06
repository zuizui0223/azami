"""Synthetic files only; no original protected observations or remote calls."""
from importlib.util import find_spec
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from analysis.v4 import species_environment_scan as scan

MOD = 'analysis.v4.species_environment_existing_inputs'
NATIVE = ('chelsa_bio01','chelsa_bio04','chelsa_bio12','chelsa_bio15',
          'chelsa_rsds_mean','chelsa_vpd_mean','chelsa_sfcwind_mean')
SEVEN = ('chelsa_bio05','chelsa_bio06','chelsa_bio13','chelsa_bio14',
         'chelsa_bio17','chelsa_bio18','topo_elevation')
NINE = tuple('chelsa_bio'+s for s in ('02','03','07','08','09','10','11','16','19'))
UNITS = {p: ('percent_isothermality' if p.endswith('03') else
             'kg_m-2_month-1' if p.endswith(('16','19')) else 'degree_C') for p in NINE}

def module():
    assert find_spec(MOD) is not None, 'existing-source adapter not implemented'
    from analysis.v4 import species_environment_existing_inputs as m
    return m


def files(tmp_path):
    rng=np.random.default_rng(60); n=150
    d=pd.DataFrame({'obs_id':[f'private_{i}' for i in range(n)],
                    'accepted_key':[f'taxon_{i//6:02}' for i in range(n)]})
    for i,p in enumerate(scan.PREDICTORS): d[p]=i+np.repeat(np.arange(25),6)+rng.normal(0,.2,n)
    d['orientation_image_vertical_angle']=30+.5*d.chelsa_bio01+rng.normal(0,1,n)
    d['orientation_image_vertical_angle__eligible']=True
    paths={k:tmp_path/(k+('.csv.gz' if k in ('native','seven','nine') else '.json'))
           for k in ('native','seven','nine','decoded_receipt','source_receipt','nine_manifest')}
    d.drop(columns=list(SEVEN+NINE)).to_csv(paths['native'],index=False)
    d[['obs_id',*SEVEN]].sample(frac=1,random_state=1).to_csv(paths['seven'],index=False)
    d[['obs_id',*NINE]].sample(frac=1,random_state=2).to_csv(paths['nine'],index=False)
    dec={'prepared_table_sha256':scan.sha256(paths['native']), 'native_observations':n,
         'environment_decoding_status':'physical_units_gsp_invalid_encoding_masked_v1',
         'environment_decoding':{p:{'source_sha256':'a'*64,'scale':.1,'offset':0.} for p in NATIVE}}
    layers={p:{'dataset':'CHELSA-bioclim-v2.1','temporal_support':'1981-2010',
                'source_bundle_sha256':'b'*64,'unit':UNITS[p],
                'finite_sampled_observations':n} for p in NINE}
    man={'schema':'azami_rf_extra32_source_manifest_v1',
         'table_sha256':scan.sha256(paths['nine']), 'layers':layers}
    paths['decoded_receipt'].write_text(json.dumps(dec));paths['nine_manifest'].write_text(json.dumps(man))
    src={'schema':'azami_bio21_source_reconstruction_receipt_20261031_v1',
         'status':'verified_source_pixels_before_fitting','source_observations':n,
         'native_s1_table_sha256':scan.sha256(paths['native']),
         'original_seven_private_table_sha256':scan.sha256(paths['seven']),
         'extra_nine_private_table_sha256':scan.sha256(paths['nine']),
         'nine_manifest_sha256':scan.sha256(paths['nine_manifest']),
         'source_variables':{p:{'source_sha256':'b'*64,'finite_count':n,
                    'source_datetime_tag':None if p=='topo_elevation' else '1981-2010',
                    'source_version_tag':None if p=='topo_elevation' else '2.1'} for p in SEVEN+NINE}}
    paths['source_receipt'].write_text(json.dumps(src))
    return paths


def reseal(paths):
    """Fixture helper, never a production authorisation/receipt generator."""
    dec=json.loads(paths['decoded_receipt'].read_text())
    dec['prepared_table_sha256']=scan.sha256(paths['native'])
    paths['decoded_receipt'].write_text(json.dumps(dec))
    man=json.loads(paths['nine_manifest'].read_text());man['table_sha256']=scan.sha256(paths['nine'])
    paths['nine_manifest'].write_text(json.dumps(man))
    src=json.loads(paths['source_receipt'].read_text())
    for k,f in [('native_s1_table_sha256','native'),('original_seven_private_table_sha256','seven'),
                ('extra_nine_private_table_sha256','nine'),('nine_manifest_sha256','nine_manifest')]:
        src[k]=scan.sha256(paths[f])
    for group,cols in [('seven',SEVEN),('nine',NINE)]:
        df=pd.read_csv(paths[group])
        for p in cols:
            src['source_variables'][p]['finite_count']=int(pd.to_numeric(df[p],errors='coerce').notna().sum())
    paths['source_receipt'].write_text(json.dumps(src))


def load(paths):
    return module().load_existing_inputs(**paths, source_receipt_sha256=scan.sha256(paths['source_receipt']),expected_n=150)


def test_existing_inputs_use_all23_and_preserve_original_order(tmp_path):
    paths=files(tmp_path);d,provenance=load(paths)
    assert set(scan.PREDICTORS)<=set(d) and len(d)==150
    assert d.obs_id.tolist()==[f'private_{i}' for i in range(150)]
    assert provenance['original_raster_bytes_rechecked'] is False
    assert provenance['source_row_count']==150


def test_pinned_source_receipt_mismatch_is_rejected(tmp_path):
    p=files(tmp_path)
    with pytest.raises(ValueError,match='SOURCE_RECEIPT_SHA'):
        module().load_existing_inputs(**p,source_receipt_sha256='0'*64,expected_n=150)


def test_changed_private_table_bytes_rejected_before_reading(tmp_path):
    p=files(tmp_path);p['seven'].write_bytes(p['seven'].read_bytes()+b'changed')
    with pytest.raises(ValueError,match='INPUT_TABLE_SHA'):load(p)


@pytest.mark.parametrize('kind',['duplicate','missing','replaced'])
def test_bad_external_identity_cannot_silently_drop_rows(tmp_path,kind):
    p=files(tmp_path);d=pd.read_csv(p['nine'])
    if kind=='duplicate': d.loc[0,'obs_id']=d.loc[1,'obs_id']
    elif kind=='missing': d=d.iloc[:-1]
    else:d.loc[0,'obs_id']='other_record'
    d.to_csv(p['nine'],index=False);reseal(p)
    with pytest.raises(ValueError,match='SOURCE_IDENTITY'):load(p)


def test_wrong_manifest_time_support_rejected(tmp_path):
    p=files(tmp_path);d=json.loads(p['nine_manifest'].read_text())
    d['layers'][NINE[0]]['temporal_support']='1991-2020'
    p['nine_manifest'].write_text(json.dumps(d));reseal(p)
    with pytest.raises(ValueError,match='LAYER_CONTRACT'):load(p)


def test_wrong_manifest_physical_unit_rejected(tmp_path):
    p=files(tmp_path);d=json.loads(p['nine_manifest'].read_text())
    d['layers']['chelsa_bio03']['unit']='degree_C'
    p['nine_manifest'].write_text(json.dumps(d));reseal(p)
    with pytest.raises(ValueError,match='LAYER_CONTRACT'):load(p)


def test_raw_not_decoded_input_rejected(tmp_path):
    p=files(tmp_path);d=json.loads(p['decoded_receipt'].read_text())
    d['environment_decoding_status']='raw'
    p['decoded_receipt'].write_text(json.dumps(d))
    with pytest.raises(ValueError,match='NATIVE_DECODE'):load(p)


def test_native_layer_source_declaration_required(tmp_path):
    p=files(tmp_path);d=json.loads(p['decoded_receipt'].read_text())
    del d['environment_decoding']['chelsa_vpd_mean']
    p['decoded_receipt'].write_text(json.dumps(d))
    with pytest.raises(ValueError,match='NATIVE_LAYER'):load(p)


def test_nonnumeric_environment_does_not_become_silent_missingness(tmp_path):
    p=files(tmp_path);d=pd.read_csv(p['seven']);d['chelsa_bio18']=d.chelsa_bio18.astype(object)
    d.loc[0,'chelsa_bio18']='not-a-number';d.to_csv(p['seven'],index=False);reseal(p)
    with pytest.raises(ValueError,match='ENVIRONMENT_VALUES'):load(p)


def test_missing_environment_values_preserve_pair_specific_support(tmp_path):
    p=files(tmp_path);d=pd.read_csv(p['nine']);d.loc[:7,'chelsa_bio19']=np.nan
    d.to_csv(p['nine'],index=False);reseal(p)
    loaded,info=load(p)
    assert len(loaded)==150 and loaded.chelsa_bio19.isna().sum()==8
    assert len(scan.aggregate_pair(loaded,'angle','chelsa_bio01'))==25
    assert info['finite_environment_counts']['chelsa_bio19']==142


def test_aggregate_run_keeps_all552_rows_and_no_private_records(tmp_path):
    p=files(tmp_path);out=tmp_path/'aggregate.json'
    r=module().run_existing_inputs(**p,source_receipt_sha256=scan.sha256(p['source_receipt']),
                                  expected_n=150,output=out)
    assert len(r['results'])==552 and r['source_row_count']==150
    assert sum(x['panel']=='marginal' and x['status']=='estimated' for x in r['results'])==23
    assert r['joint_correction_status']=='not_estimated'
    for private in ('private_0','taxon_00','"obs_id"','"accepted_key"','"latitude"','"longitude"'):
        assert private not in out.read_text()


def test_output_cannot_overwrite_an_input(tmp_path):
    p=files(tmp_path);before=p['native'].read_bytes()
    with pytest.raises(ValueError,match='OUTPUT_OVERWRITES_INPUT'):
        module().run_existing_inputs(**p,source_receipt_sha256=scan.sha256(p['source_receipt']),
                                    expected_n=150,output=p['native'])
    assert p['native'].read_bytes()==before

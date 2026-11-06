"""Local-only adapter from existing S1/BIO21 files to the 23-axis atlas.

No network, credentials, decryption, data recovery or new source extraction.
Only already-authorized local files are read. A caller-pinned public receipt
links private file bytes and manifest; it does not reauthenticate raw rasters.
Only the non-identifying scan summary is written, never a joined source table.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd
from analysis.v4 import species_environment_scan as atlas

NATIVE7=('chelsa_bio01','chelsa_bio04','chelsa_bio12','chelsa_bio15',
         'chelsa_rsds_mean','chelsa_vpd_mean','chelsa_sfcwind_mean')
SEVEN=('chelsa_bio05','chelsa_bio06','chelsa_bio13','chelsa_bio14',
       'chelsa_bio17','chelsa_bio18','topo_elevation')
NINE=tuple('chelsa_bio'+s for s in ('02','03','07','08','09','10','11','16','19'))
UNITS={p:('percent_isothermality' if p.endswith('03') else
          'kg_m-2_month-1' if p.endswith(('16','19')) else 'degree_C') for p in NINE}
HEX=re.compile(r'[a-f0-9]{64}\Z')
assert set(NATIVE7+SEVEN+NINE)==set(atlas.PREDICTORS)
assert len(NATIVE7+SEVEN+NINE)==23


def _json(path:Path)->dict:
    value=json.loads(Path(path).read_text(encoding='utf-8'))
    if not isinstance(value,dict):raise ValueError('INPUT_RECEIPT_NOT_OBJECT')
    return value


def _finite_numeric(frame:pd.DataFrame,columns:tuple[str,...])->None:
    if not set(columns)<=set(frame):raise ValueError('ENVIRONMENT_COLUMNS_MISSING')
    for col in columns:
        v=pd.to_numeric(frame[col],errors='coerce')
        if (frame[col].notna() & v.isna()).any() or np.isinf(v.to_numpy(float)).any():
            raise ValueError('ENVIRONMENT_VALUES_INVALID')
        frame[col]=v


def _external(path:Path,columns:tuple[str,...],ids:pd.Series,expected_n:int)->pd.DataFrame:
    d=pd.read_csv(path,dtype={'obs_id':str},low_memory=False)
    if (set(d)!={'obs_id',*columns} or len(d)!=expected_n or
            d.obs_id.isna().any() or not d.obs_id.is_unique or
            set(d.obs_id)!=set(ids)):
        raise ValueError('EXTERNAL_SOURCE_IDENTITY_INVALID')
    _finite_numeric(d,columns)
    return d


def load_existing_inputs(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                         source_receipt:Path,nine_manifest:Path,*,
                         source_receipt_sha256:str,expected_n:int=atlas.EXPECTED_N
                         )->tuple[pd.DataFrame,dict]:
    """PRIVATE in-memory joined frame; metadata return contains aggregates only."""
    paths={k:Path(v) for k,v in dict(native=native,decoded_receipt=decoded_receipt,
            seven=seven,nine=nine,source_receipt=source_receipt,nine_manifest=nine_manifest).items()}
    if (not HEX.fullmatch(str(source_receipt_sha256)) or
            atlas.sha256(paths['source_receipt'])!=source_receipt_sha256):
        raise ValueError('SOURCE_RECEIPT_SHA_MISMATCH')
    src=_json(paths['source_receipt']);dec=_json(paths['decoded_receipt'])
    man=_json(paths['nine_manifest'])
    if (src.get('schema')!='azami_bio21_source_reconstruction_receipt_20261031_v1' or
            src.get('status')!='verified_source_pixels_before_fitting' or
            src.get('source_observations')!=expected_n):
        raise ValueError('SOURCE_RECEIPT_CONTRACT_INVALID')
    digests={k:atlas.sha256(paths[k]) for k in ('native','seven','nine','nine_manifest','decoded_receipt')}
    checks={'native':'native_s1_table_sha256','seven':'original_seven_private_table_sha256',
            'nine':'extra_nine_private_table_sha256','nine_manifest':'nine_manifest_sha256'}
    if any(src.get(field)!=digests[name] for name,field in checks.items()):
        raise ValueError('INPUT_TABLE_SHA_MISMATCH')
    if (dec.get('prepared_table_sha256')!=digests['native'] or
            dec.get('native_observations')!=expected_n or
            dec.get('environment_decoding_status')!='physical_units_gsp_invalid_encoding_masked_v1'):
        raise ValueError('NATIVE_DECODE_CONTRACT_INVALID')
    for p in NATIVE7:
        meta=dec.get('environment_decoding',{}).get(p,{})
        try:scale,offset=float(meta['scale']),float(meta['offset'])
        except (KeyError,TypeError,ValueError):raise ValueError('NATIVE_LAYER_METADATA_INVALID') from None
        if not HEX.fullmatch(str(meta.get('source_sha256',''))) or not np.isfinite([scale,offset]).all() or scale<=0:
            raise ValueError('NATIVE_LAYER_METADATA_INVALID')
    if (man.get('schema')!='azami_rf_extra32_source_manifest_v1' or
            man.get('table_sha256')!=digests['nine'] or set(man.get('layers',{}))!=set(NINE) or
            set(src.get('source_variables',{}))!=set(SEVEN+NINE)):
        raise ValueError('LAYER_CONTRACT_INVALID')
    for p in SEVEN+NINE:
        spec=src['source_variables'][p]
        if not HEX.fullmatch(str(spec.get('source_sha256',''))):raise ValueError('LAYER_CONTRACT_INVALID')
        if p!='topo_elevation' and (spec.get('source_datetime_tag') not in (None,'','1981-2010') or
                                   spec.get('source_version_tag') not in (None,'','2.1')):
            raise ValueError('LAYER_CONTRACT_INVALID')
        if p in NINE:
            layer=man['layers'][p]
            if (layer.get('dataset')!='CHELSA-bioclim-v2.1' or
                    layer.get('temporal_support')!='1981-2010' or layer.get('unit')!=UNITS[p] or
                    layer.get('source_bundle_sha256')!=spec['source_sha256']):
                raise ValueError('LAYER_CONTRACT_INVALID')
    endpoints=[v for group in atlas.FEATURES.values() for v in group]
    columns={'obs_id','accepted_key',*endpoints,*(v+'__eligible' for v in endpoints),*NATIVE7}
    header=set(pd.read_csv(paths['native'],nrows=0).columns)
    if header & set(SEVEN+NINE):raise ValueError('NATIVE_ALREADY_HAS_EXTERNAL_COLUMNS')
    d=pd.read_csv(paths['native'],usecols=lambda c:c in columns,
                  dtype={'obs_id':str,'accepted_key':str},low_memory=False)
    atlas.validate_observations(d)
    if len(d)!=expected_n:raise ValueError('NATIVE_SOURCE_IDENTITY_COUNT_INVALID')
    _finite_numeric(d,NATIVE7)
    a=_external(paths['seven'],SEVEN,d.obs_id,expected_n)
    b=_external(paths['nine'],NINE,d.obs_id,expected_n)
    ids=d.obs_id.copy()
    for block in (a,b):d=d.merge(block,on='obs_id',how='left',validate='one_to_one',sort=False)
    if len(d)!=expected_n or not d.obs_id.equals(ids):raise ValueError('JOIN_SOURCE_IDENTITY_CHANGED')
    counts={p:int(d[p].notna().sum()) for p in atlas.PREDICTORS}
    if any(src['source_variables'][p].get('finite_count')!=counts[p] for p in SEVEN+NINE):
        raise ValueError('SOURCE_COVERAGE_RECEIPT_MISMATCH')
    return d,{'schema':'azami_existing_local_inputs_for_species_scan_v1',
        'source_row_count':int(expected_n),'source_receipt_sha256':source_receipt_sha256,
        'input_digests':digests,'finite_environment_counts':counts,
        'global_complete_case_filter_applied':False,
        'original_raster_bytes_rechecked':False,
        'scope':'byte/identity checks against caller-pinned upstream extraction receipt; no source recovery',
        'joined_source_table_written':False}


def run_existing_inputs(native:Path,decoded_receipt:Path,seven:Path,nine:Path,
                        source_receipt:Path,nine_manifest:Path,output:Path,*,
                        source_receipt_sha256:str,expected_n:int=atlas.EXPECTED_N,
                        covariance_path:Path|None=None)->dict:
    inputs=[native,decoded_receipt,seven,nine,source_receipt,nine_manifest]
    if covariance_path is not None:inputs.append(covariance_path)
    if Path(output).resolve() in {Path(p).resolve() for p in inputs}:
        raise ValueError('OUTPUT_OVERWRITES_INPUT')
    d,meta=load_existing_inputs(native,decoded_receipt,seven,nine,source_receipt,nine_manifest,
              source_receipt_sha256=source_receipt_sha256,expected_n=expected_n)
    cov=None if covariance_path is None else atlas.covariance_from_payload(_json(covariance_path))
    result=atlas.scan(d,covariance=cov)
    result.update(source_row_count=int(expected_n),existing_input_provenance=meta,
                  covariance_file_sha256=None if covariance_path is None else atlas.sha256(covariance_path))
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('native','decoded-receipt','seven','nine','source-receipt','nine-manifest','output'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--source-receipt-sha256',required=True)
    p.add_argument('--covariance',type=Path)
    args=vars(p.parse_args());args['covariance_path']=args.pop('covariance')
    r=run_existing_inputs(**args)
    print(json.dumps({'status':r['analysis_status'],'planned_pairs':r['n_planned_trait_exposure_pairs'],
                      'output_rows':len(r['results']),'joint_correction_status':r['joint_correction_status']}))

if __name__=='__main__':main()

"""Native-first 665,115-observation source and CHELSA replay (exploratory v4).

Never publish observation identifiers, image URLs or exact coordinates. Metadata,
WCVP authority and raster caches are checksum-pinned. 5-degree cells are storage
groups only; original-pixel values and nodata are retained without imputation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
from analysis.v4.native_scope import classify_roster
from analysis.v4.native_climate_cache import sample_tiles

ORIGINAL_SHA='2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae'
AUTHORITY_SHA='a5bf23efffbf8e8f9459d5c3391a0136adcedd8f9dbff0117caa969ab545c2bd'
NATIVE_N=313350
SOURCE_N=665115
PHOTO_N=1122854
CLIMATE={
'bio01':('chelsa_bio01','0e3b83ebfb7de1865662cde2c2c2395fa7737460996e66354f87391f8f0867ab'),
'bio04':('chelsa_bio04','dac232683b2a72b9c4f62eab0f341be7927e5e3287a29dbf3a162b893c2307e3'),
'bio12':('chelsa_bio12','71e7edeb0d4c0042e973ed97672a1593bef98a608c917792006e22722e288ada'),
'bio15':('chelsa_bio15','9f255a7c756fc83a5dd81e6c13d8fee63c53698d93bffb65ae13288e4eebb445'),
'rsds_mean':('chelsa_rsds_mean','2193bd04b979f6684f7bce7c670074af63072d270a20c4dc2fd0bd5934862df4'),
'vpd_mean':('chelsa_vpd_mean','2b92152e2ed33436ea6c242497d2b0e61c34711e958a7f4eb3620d8bf3ad1265'),
'sfcwind_mean':('chelsa_sfcwind_mean','f383fa07776984ff193864fa8d44a3b0baa079b64a798e014b9e364b9b1f67b0'),
'gsp':('chelsa_gsp','8259581896660edd7c07667b97640f71f174cab484eeadfe1292588ab8a96e8d'),
'npp':('chelsa_npp','cd6ab1fd63568ced9b6b18fe4e0d36f65ccdf00b2403e7fbcba9ab9578e76c5a')}
FIELDS=['obs_id','photo_id','taxon_name','taxon_rank','observed_on','captive',
        'latitude','longitude','coordinate_usable_for_environment','obscured',
        'quality_grade','positional_accuracy','geoprivacy']


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda:source.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def add_explicit_non_species_placeholders(master,authority):
    authority=dict(authority)
    resolved=list(authority['resolutions'])
    queried={r['input_name'] for r in resolved}
    missing=master.loc[~master.taxon_name.isin(queried),['taxon_name','taxon_rank']]
    if missing.taxon_rank.eq('species').any():
        raise ValueError('Missing exact species authority query')
    for name in sorted(missing.taxon_name.unique()):
        resolved.append({'input_name':name,'resolution_status':'unresolved_taxon',
                         'accepted_key':None,'accepted_name':'',
                         'n_exact_records':0,'n_accepted_targets':0})
    authority['resolutions']=resolved
    return authority,len(queried),int(missing.taxon_name.nunique())


def classify_all(metadata_zip,authority_json,private_dir):
    metadata_zip=Path(metadata_zip);authority_json=Path(authority_json)
    private_dir=Path(private_dir);private_dir.mkdir(parents=True,exist_ok=True)
    if digest(metadata_zip)!=ORIGINAL_SHA:raise ValueError('Original upstream ZIP changed')
    if digest(authority_json)!=AUTHORITY_SHA:raise ValueError('350 species authority JSON changed')
    authority=json.loads(authority_json.read_text())
    if len(authority['resolutions'])!=350:raise ValueError('Expected all 350 species queried')
    with zipfile.ZipFile(metadata_zip) as z,z.open('photo_metadata_merged.csv') as stream:
        photos=pd.read_csv(stream,usecols=FIELDS,dtype={'obs_id':str,'photo_id':str},low_memory=False)
    if len(photos)!=PHOTO_N or photos.obs_id.nunique()!=SOURCE_N or not photos.photo_id.is_unique:
        raise ValueError('Unexpected source photo/observation identity')
    for col in ['taxon_name','taxon_rank','latitude','longitude','captive','obscured']:
        if photos.groupby('obs_id')[col].nunique(dropna=False).gt(1).any():
            raise ValueError('Conflicting metadata within observation: '+col)
    master=photos.drop_duplicates('obs_id').drop(columns=['photo_id']).copy()
    master['source_photo_count']=master.obs_id.map(photos.groupby('obs_id').size())
    authority,queried,unresolved=add_explicit_non_species_placeholders(master,authority)
    classified=classify_roster(master,authority)
    native=classified.loc[classified.native_primary].copy()
    if len(native)!=NATIVE_N:
        raise ValueError(f'Native count shifted: {len(native)} vs {NATIVE_N}')
    if native.accepted_key.nunique()!=306:raise ValueError('Accepted taxon census shifted')
    native[['obs_id','taxon_name','accepted_key','latitude','longitude','observed_on',
            'source_photo_count']].to_csv(private_dir/'native350_roster.csv.gz',index=False)
    classified.to_csv(private_dir/'all_upstream_native350_status.csv.gz',index=False)
    aggregate={'source_observations':SOURCE_N,'source_photos':PHOTO_N,
       'queried_species_names':queried,'non_species_names_explicitly_unresolved':unresolved,
       'native_observations':len(native),'native_photos':int(native.source_photo_count.sum()),
       'native_source_names':int(native.taxon_name.nunique()),
       'native_accepted_taxa':int(native.accepted_key.nunique()),
       'range_status_counts':{str(k):int(v) for k,v in classified.native_range_status.value_counts().items()},
       'source_metadata_sha256':ORIGINAL_SHA,'native_authority_sha256':AUTHORITY_SHA,
       'native_roster_sha256':digest(private_dir/'native350_roster.csv.gz'),
       'private_observation_rows_uploaded_to_public_CI':False,
       'trait_or_detector_values_used_for_native_selection':False}
    (private_dir/'source_receipt.json').write_text(json.dumps(aggregate,indent=2)+'\n')
    return aggregate


def sample_all(native_roster,cache_root,private_dir):
    native_roster=Path(native_roster);cache_root=Path(cache_root)
    private_dir=Path(private_dir);private_dir.mkdir(parents=True,exist_ok=True)
    roster=pd.read_csv(native_roster,dtype={'obs_id':str})
    if len(roster)!=NATIVE_N or not roster.obs_id.is_unique:
        raise ValueError('Expected whole 313350 native roster, no repeats')
    out=roster[['obs_id','taxon_name','accepted_key','latitude','longitude']].copy()
    tile_sets=[];variable_counts={}
    for var,(key,_) in CLIMATE.items():
        directory=cache_root/var
        meta=json.loads((directory/'grid.json').read_text())
        if meta['resampling']!='none' or len(meta['tiles'])!=316:
            raise ValueError('Changed original resolution or old tile contract')
        for filename,expected in meta['cache_sha256'].items():
            if digest(directory/filename)!=expected:
                raise ValueError('Tile digest changed: '+var+'/'+filename)
        value=sample_tiles(directory,out.latitude.to_numpy(),out.longitude.to_numpy())
        out[key]=value
        variable_counts[var]={'complete':int(np.isfinite(value).sum()),
                              'missing':int(np.isnan(value).sum())}
        tile_sets.append(set(meta['tiles']))
    if not all(s==tile_sets[0] for s in tile_sets):raise ValueError('Raster variable tile sets differ')
    cells=np.floor((out.latitude+90)/5).astype(int).astype(str)+'_'+np.floor((out.longitude+180)/5).astype(int).astype(str)
    outside=int((~cells.isin(tile_sets[0])).sum())
    complete=np.isfinite(out[[p[0] for p in CLIMATE.values()]].to_numpy()).all(1)
    missing=int((~complete).sum())
    if (outside,missing)!=(106,110):
        raise ValueError(f'Unexpected 5-degree tile/nodata gap: {outside}/{missing}')
    path=private_dir/'native350_chelsa_9.csv.gz'
    out.to_csv(path,index=False)
    receipt={'native_source_observations':len(out),'original_pixel_cache_cells':len(tile_sets[0]),
       'outside_cached_cells_observations':outside,
       'outside_cached_cells_unique':int((~cells.isin(tile_sets[0])).groupby(cells).any().sum()),
       'complete_all_nine':int(complete.sum()),'missing_any':missing,
       'complete_fraction':float(complete.mean()),'variables':variable_counts,
       'five_degree_tiles_are_storage_only':True,'nodata_imputed':False,
       'models_recomputed_in_this_step':False,'private_coordinates_uploaded_to_public_CI':False,
       'environment_csv_sha256':digest(path)}
    (private_dir/'climate_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def measurement_frame(roster, environment, measured, archived, contract):
    """Left-join all native observations; each endpoint retains its own QC."""
    for label, table in [('roster', roster), ('environment', environment), ('measurements', measured)]:
        if table.obs_id.isna().any() or not table.obs_id.is_unique:
            raise ValueError('Duplicate or missing observation ID: '+label)
    if set(roster.obs_id) != set(environment.obs_id):
        raise ValueError('Environment must cover the complete native roster')
    frame = roster.merge(environment[['obs_id']+[v[0] for v in CLIMATE.values()]],
                         on='obs_id', how='left', validate='one_to_one')
    frame['archived_measurement_row'] = frame.obs_id.isin(measured.obs_id)
    measured = measured.set_index('obs_id')
    archived = archived.loc[archived.obs_id.isin(frame.obs_id)].copy()
    if archived.duplicated(['obs_id','endpoint_id']).any():
        raise ValueError('Duplicate archived endpoint identity')
    for source in (measured.reset_index(), archived):
        identity = source[['obs_id','taxon_name']].drop_duplicates().merge(
            roster[['obs_id','taxon_name']], on='obs_id', suffixes=('_source',''))
        if identity.taxon_name_source.ne(identity.taxon_name).any():
            raise ValueError('Measurement and native source taxon labels differ')
    counts = {}
    for row in contract.to_dict('records'):
        endpoint, column = row['endpoint_id'], row['observation_variable']
        if column in measured:
            values = frame.obs_id.map(measured[column])
            usable = frame.obs_id.map(measured[row['measurement_variable']+'_n_usable_heads'])
            eligible = np.isfinite(values) & usable.gt(0)
            origin = 'full_archived_observation_table'
        else:
            records = archived.loc[archived.endpoint_id.eq(endpoint)].set_index('obs_id')
            values = frame.obs_id.map(records.value)
            flags = records.analysis_eligible.astype(str).str.lower()
            if not flags.isin(['true','false']).all():
                raise ValueError('Unrecognised endpoint eligibility')
            eligible = frame.obs_id.map(flags.eq('true')).eq(True) & np.isfinite(values)
            origin = 'zenodo_endpoint_table'
        frame[endpoint] = values
        frame[endpoint+'__eligible'] = eligible
        counts[endpoint] = {'finite':int(np.isfinite(values).sum()),
                            'eligible':int(eligible.sum()), 'source':origin}
    main_environment = ['chelsa_bio01','chelsa_rsds_mean','chelsa_vpd_mean','chelsa_bio12']
    climate_ok = np.isfinite(frame[main_environment]).all(axis=1)
    # Exposure means are covariates, never averaged trait responses.
    for name in main_environment:
        frame[name+'__between'] = frame.groupby('accepted_key')[name].transform('mean')
        frame[name+'__within'] = frame[name]-frame[name+'__between']
    for endpoint in counts:
        counts[endpoint]['eligible_with_four_climate'] = int(
            (frame[endpoint+'__eligible'] & climate_ok).sum())
    def predictor_summary(mask):
        x = frame.loc[mask & climate_ok,main_environment].to_numpy(float)
        result = {'n':len(x),'variables':main_environment}
        if len(x)<2 or (x.std(axis=0)<=1e-12).any():
            return dict(result,status='insufficient_variation')
        correlation = np.corrcoef(x,rowvar=False)
        full_rank = np.linalg.matrix_rank(correlation)==len(main_environment)
        return dict(result,status='estimated' if full_rank else 'singular',
                    correlation=correlation.tolist(),
                    vif=np.diag(np.linalg.inv(correlation)).tolist() if full_rank else None)
    predictor_diagnostics = {'all_native':predictor_summary(pd.Series(True,index=frame.index))}
    for group, names in {
        'colour_any_coordinate':['corolla_lab_lightness','corolla_lab_chroma','corolla_hue_sin','corolla_hue_cos'],
        'image_orientation':['orientation_image_vertical_angle'],
        'legacy_bract_projection_not_opening':['bract_projection_p95','bract_projection_peak_density'],
    }.items():
        columns = [n+'__eligible' for n in names if n in counts]
        if columns:
            predictor_diagnostics[group] = predictor_summary(frame[columns].any(axis=1))
    report = {'native_observations':len(frame),
              'with_archived_measurement_row':int(frame.archived_measurement_row.sum()),
              'without_archived_measurement_row':int((~frame.archived_measurement_row).sum()),
              'complete_four_climate':int(climate_ok.sum()),
              'endpoints':counts, 'trait_responses_averaged':False,
              'four_predictor_diagnostics':predictor_diagnostics,
              'all_endpoint_complete_case_required':False,
              'phyllary_opening_angle_measured':False,
              'unit':'inherited observation-level values; no further trait aggregation',
              'status':'inputs_prepared_not_models_or_new_measurements'}
    return frame, report


def join_measurements(private_dir, measurements, traits_zip):
    from reproducibility.run_current_analysis import INPUTS
    private_dir = Path(private_dir)
    expected = 'b49a093c849928c456d6747ec51772b5a6b78debb2ae5eb5dbe2a91b08f66dd7'
    if digest(measurements) != expected or digest(traits_zip) != INPUTS['continuous'][1]:
        raise ValueError('Changed measurement input bytes')
    roster = pd.read_csv(private_dir/'native350_roster.csv.gz',dtype={'obs_id':str})
    environment = pd.read_csv(private_dir/'native350_chelsa_9.csv.gz',dtype={'obs_id':str})
    measured = pd.read_csv(measurements,dtype={'obs_id':str},low_memory=False)
    with zipfile.ZipFile(traits_zip) as z, z.open('universe/continuous_trait_universe_observation_long.csv') as f:
        archived = pd.read_csv(f, usecols=['obs_id','taxon_name','endpoint_id','value','analysis_eligible'],
                               dtype={'obs_id':str},low_memory=False)
    root = Path(__file__).resolve().parents[2]
    contract = pd.read_csv(root/'legacy/ch1_global/v2/ontology/ch1_continuous_trait_contract.csv')
    frame, report = measurement_frame(roster,environment,measured,archived,contract)
    target = private_dir/'native_observation_measurements.csv.gz'
    frame.to_csv(target,index=False,compression={'method':'gzip','mtime':0})
    report.update(measurements_sha256=expected, zenodo_traits_sha256=INPUTS['continuous'][1],
                  prepared_table_sha256=digest(target))
    (private_dir/'measurement_receipt.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('classify');a.add_argument('--metadata-zip',type=Path,required=True)
    a.add_argument('--authority-json',type=Path,required=True);a.add_argument('--private-dir',type=Path,required=True)
    a=sub.add_parser('sample');a.add_argument('--native-roster',type=Path,required=True)
    a.add_argument('--cache-root',type=Path,required=True);a.add_argument('--private-dir',type=Path,required=True)
    a=sub.add_parser('join');a.add_argument('--private-dir',type=Path,required=True)
    a.add_argument('--measurements',type=Path,required=True);a.add_argument('--traits-zip',type=Path,required=True)
    args=p.parse_args()
    if args.command=='classify':
        r=classify_all(args.metadata_zip,args.authority_json,args.private_dir)
    elif args.command=='sample':r=sample_all(args.native_roster,args.cache_root,args.private_dir)
    else:r=join_measurements(args.private_dir,args.measurements,args.traits_zip)
    print(json.dumps({k:v for k,v in r.items() if k not in ['range_status_counts','variables']},indent=2))

if __name__=='__main__':main()

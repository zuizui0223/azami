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


def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('classify');a.add_argument('--metadata-zip',type=Path,required=True)
    a.add_argument('--authority-json',type=Path,required=True);a.add_argument('--private-dir',type=Path,required=True)
    a=sub.add_parser('sample');a.add_argument('--native-roster',type=Path,required=True)
    a.add_argument('--cache-root',type=Path,required=True);a.add_argument('--private-dir',type=Path,required=True)
    args=p.parse_args()
    if args.command=='classify':
        r=classify_all(args.metadata_zip,args.authority_json,args.private_dir)
    else:r=sample_all(args.native_roster,args.cache_root,args.private_dir)
    print(json.dumps({k:v for k,v in r.items() if k not in ['range_status_counts','variables']},indent=2))

if __name__=='__main__':main()

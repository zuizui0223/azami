"""Build native-first full-source master tables before any phenotype selection.

Native means WCVP's TDWG3 region-level status, NOT proof of local ancestry.
Source snapshot 8066010557; all ranks retained in audit, only source species
rank eligible. Missing/ambiguous native evidence never becomes introduced.
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
from analysis.v4.native_climate_cache import sample_tiles, VARIABLES

METADATA_SHA = '2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae'
COLS = ['obs_id','taxon_name','taxon_rank','captive','latitude','longitude',
        'coordinate_usable_for_environment','obscured','positional_accuracy',
        'quality_grade','geoprivacy','observed_on','user_id']
IDENTITY = ['taxon_name','taxon_rank','captive','latitude','longitude',
            'coordinate_usable_for_environment','obscured','geoprivacy']


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for a in iter(lambda:f.read(1<<20),b''): h.update(a)
    return h.hexdigest()


def unique_roster(frame):
    if frame.obs_id.isna().any(): raise ValueError('Missing observation ID')
    cols=[c for c in IDENTITY if c in frame]
    if frame.groupby('obs_id')[cols].nunique(dropna=False).gt(1).any().any():
        raise ValueError('Conflicting metadata for a source observation')
    return frame.drop_duplicates('obs_id').reset_index(drop=True)


def classify_all(roster, source):
    if not roster.obs_id.is_unique: raise ValueError('Observation IDs must be unique')
    species=roster.loc[roster.taxon_rank.eq('species')].copy()
    classified=classify_roster(species,source)
    # A contradictory privacy field must not bypass the explicit-coordinate gate.
    private=classified.geoprivacy.fillna('').astype(str).str.lower().isin(['private','obscured']) if 'geoprivacy' in classified else pd.Series(False,index=classified.index)
    classified.loc[private,'native_primary']=False
    classified.loc[private,'eligibility_reason']='private_or_obscured_geoprivacy'
    other=roster.loc[~roster.taxon_rank.eq('species')].copy()
    other['native_primary']=False
    other['native_range_status']='not_evaluated_rank'
    other['resolution_status']='not_evaluated_rank'
    other['eligibility_reason']='rank_not_species'
    return pd.concat([classified,other],ignore_index=True).sort_values('obs_id').reset_index(drop=True)


def physical_values(raw, metadata):
    return np.asarray(raw,float)*float(metadata['scale'])+float(metadata['offset'])


def native_measurement_join(roster, traits):
    if not roster.obs_id.is_unique or not traits.obs_id.is_unique:
        raise ValueError('Both join tables must have unique observation IDs')
    primary=roster.loc[roster.native_primary].copy()
    return primary.merge(traits,on='obs_id',how='inner',validate='one_to_one',suffixes=('','_legacy'))


def build(metadata_zip, source_json, outdir, climate=None):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    if digest(metadata_zip)!=METADATA_SHA: raise ValueError('Wrong master source archive')
    source=json.loads(Path(source_json).read_text())
    pieces=[];photo_rows=0
    with zipfile.ZipFile(metadata_zip) as z, z.open('photo_metadata_merged.csv') as f:
        for a in pd.read_csv(f,usecols=COLS,dtype={'obs_id':str,'user_id':str},chunksize=100000,low_memory=False):
            photo_rows+=len(a);pieces.append(unique_roster(a))
    roster=unique_roster(pd.concat(pieces,ignore_index=True))
    assert photo_rows==1122854 and len(roster)==665115
    classified=classify_all(roster,source)
    print('Native classification complete',len(classified),flush=True)
    native=classified.loc[classified.native_primary].copy()
    native['analysis_taxon']=native.accepted_key.map(lambda x:f'wcvp:{int(x)}')
    # Never export private coordinates in the audit ledger.
    unsafe=classified.coordinate_usable_for_environment.astype(str).str.lower().ne('true') | classified.obscured.astype(str).str.lower().ne('false') | classified.geoprivacy.fillna('').isin(['private','obscured'])
    classified.loc[unsafe,['latitude','longitude']]=np.nan
    classified.to_csv(outdir/'all_observation_native_ledger.csv.gz',index=False,compression={'method':'gzip','compresslevel':1})
    native.to_csv(outdir/'native_observations.csv.gz',index=False,compression={'method':'gzip','compresslevel':1})
    native.to_pickle(outdir/'native_observations.pkl')
    # Keep every photo in eligible observations, with licensing separate.
    ids=set(native.obs_id);queue=[];allcols=None
    with zipfile.ZipFile(metadata_zip) as z,z.open('photo_metadata_merged.csv') as f:
        for a in pd.read_csv(f,dtype={'obs_id':str,'photo_id':str},chunksize=100000,low_memory=False):
            take=a.loc[a.obs_id.isin(ids),['obs_id','photo_id','photo_index','photo_license_code','photo_attribution','raw_image_url','large_image_url']].copy()
            take['acquisition_license_eligible']=take.photo_license_code.isin(['cc0','cc-by','cc-by-sa','cc-by-nc','cc-by-nc-sa'])
            queue.append(take)
    photos=pd.concat(queue,ignore_index=True)
    print('Native photo rows',len(photos),flush=True)
    photos.to_csv(outdir/'native_photo_manifest.csv.gz',index=False,compression={'method':'gzip','compresslevel':1})
    env=native[['obs_id','analysis_taxon','taxon_name','latitude','longitude']].copy()
    climate_report={}
    if climate:
        for var in VARIABLES:
            directory=Path(climate)/var
            m=json.loads((directory/'grid.json').read_text())
            for fn,sha in m['cache_sha256'].items():
                if digest(directory/fn)!=sha: raise ValueError(f'Climate checksum mismatch: {var}/{fn}')
            raw=sample_tiles(directory,native.latitude,native.longitude)
            print('Climate',var,'finite',int(np.isfinite(raw).sum()),flush=True)
            env['chelsa_'+var]=physical_values(raw,m)
            a=np.floor((native.latitude.to_numpy()+90)/5).astype(int);b=np.floor((native.longitude.to_numpy()+180)/5).astype(int)
            cached=np.array([f'{x}_{y}' in m['tiles'] for x,y in zip(a,b)])
            climate_report[var]={'n_finite':int(np.isfinite(raw).sum()),'outside_cached_cells':int((~cached).sum()),'nodata_inside_cached_cells':int((cached & ~np.isfinite(raw)).sum()),'source_sha256':m['source_sha256'],'grid_sha256':digest(directory/'grid.json'),'scale':m['scale'],'offset':m['offset'],'unit':m['tags'].get('variable_unit'),'resolution_degrees':m['transform'][0]}
        env['all9_available']=env[['chelsa_'+v for v in VARIABLES]].notna().all(axis=1)
        env.to_csv(outdir/'native_environment.csv.gz',index=False,compression={'method':'gzip','compresslevel':1})
        env.to_pickle(outdir/'native_environment.pkl')
        env.loc[~env.all9_available].to_csv(outdir/'native_environment_gaps.csv',index=False)
    counts=classified.groupby(['taxon_name','taxon_rank','eligibility_reason']).size().unstack(fill_value=0).reset_index()
    counts.to_csv(outdir/'native_taxon_flow.csv',index=False)
    report={'scope':'all historical upstream metadata, native-first; no phenotype input',
            'source_photo_rows':photo_rows,'source_observations':len(roster),'species_rank_observations':int(roster.taxon_rank.eq('species').sum()),
            'source_species_names':int(roster.loc[roster.taxon_rank.eq('species'),'taxon_name'].nunique()),
            'native_observations':len(native),'native_photo_rows':len(photos),'native_source_taxa':native.taxon_name.nunique(),'native_accepted_taxa':native.accepted_key.nunique(),
            'range_status_counts':classified.native_range_status.value_counts().to_dict(),'eligibility_counts':classified.eligibility_reason.value_counts().to_dict(),
            'native_top_taxa':native.taxon_name.value_counts().head(10).to_dict(),'photos_by_license':photos.photo_license_code.fillna('missing').value_counts().to_dict(),
            'coordinate_accuracy_known':int(pd.to_numeric(native.positional_accuracy,errors='coerce').notna().sum()),
            'metadata_zip_sha256':METADATA_SHA,'native_source_sha256':digest(source_json),'climate':climate_report,
            'environment_all9':int(env.all9_available.sum()) if climate else None,
            'new_pixel_measurements':0,'native_definition':'WCVP TDWG3 native region, source noncaptive and usable non-obscured coordinates; not proof of local indigenous population',
            'subspecies_or_genus_coercion':False,'licence_used_for_native_filter':False,'phenotype_used_for_native_filter':False}
    (outdir/'native_full_receipt.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--metadata-zip',type=Path,required=True);p.add_argument('--native-sources',type=Path,required=True);p.add_argument('--out-dir',type=Path,required=True);p.add_argument('--climate-root',type=Path)
    a=p.parse_args();print(json.dumps(build(a.metadata_zip,a.native_sources,a.out_dir,a.climate_root),indent=2))

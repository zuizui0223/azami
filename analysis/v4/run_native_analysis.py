"""Replay native-first analyses from ALL upstream observation metadata.

Raw photographs are not remeasured. Supplied measurements are the archived
image-derived observation values; source evidence and climate values are frozen.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v4.native_scope import classify_roster
from analysis.v4.native_statistics import integration
from analysis.v4.native_environment import analyse,VARS

SOURCE_SHA='d3bf06b5237a371094ed4c52322d270fb653b689056886bbd62145105ec3d046'
TRAITS=['orientation_image_vertical_angle','corolla_lab_lightness','corolla_lab_chroma','corolla_hue_sin','corolla_hue_cos']
RENAMES={'orientation_angle_degrees_median':TRAITS[0],**{c+'_median':c for c in TRAITS[1:]}}
ROSTER_COLUMNS=['obs_id','taxon_name','quality_grade','latitude','longitude','coordinate_usable_for_environment',
                'positional_accuracy','obscured','observed_on','captive','user_id']


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def run(roster_path,measurement_path,evidence_path,environment_path,outdir,replicates=1000):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    if sha(evidence_path)!=SOURCE_SHA:raise ValueError('Native source snapshot does not match frozen evidence')
    source=pd.read_csv(roster_path,dtype={'obs_id':str},usecols=ROSTER_COLUMNS)
    if len(source)!=460036 or source.obs_id.nunique()!=460036 or source.taxon_name.nunique()!=288:
        raise ValueError('Expected all 460036 upstream observations and 288 source labels')
    classified=classify_roster(source,json.loads(Path(evidence_path).read_text()))
    classified.to_csv(outdir/'all_observations_native_status.csv.gz',index=False)
    n=classified[classified.native_primary].copy();native=set(n.obs_id)
    n.to_csv(outdir/'native_observation_roster.csv.gz',index=False)
    old=pd.read_csv(measurement_path,dtype={'obs_id':str},low_memory=False).rename(columns=RENAMES)
    if not old.obs_id.is_unique:raise ValueError('Duplicate source measurement observation')
    f=old[old.obs_id.isin(native)].merge(n[['obs_id','taxon_name','accepted_key','accepted_name','user_id']],
                                       on=['obs_id','taxon_name'],validate='one_to_one')
    if len(f)!=old.obs_id.isin(native).sum():raise ValueError('Source label mismatch in measurement join')
    f['source_taxon_name']=f.taxon_name;f['taxon_name']='wcvp:'+f.accepted_key.astype(str)
    lon=np.floor(((f.longitude+180)%360)/2).astype(int)
    lat=np.minimum(89,np.floor((np.sin(np.radians(f.latitude))+1)*45)).astype(int)
    f['cell']=lon.astype(str)+'_'+lat.astype(str)
    f.to_csv(outdir/'native_measurements.csv.gz',index=False)
    e=pd.read_csv(environment_path,dtype={'obs_id':str},low_memory=False)
    if not e.obs_id.is_unique or set(e.obs_id)!=native:raise ValueError('Climate table must include every native target, including nodata rows')
    check=n[['obs_id','latitude','longitude']].merge(e[['obs_id','latitude','longitude']],on='obs_id',suffixes=('_source','_env'),validate='one_to_one')
    for c in ['latitude','longitude']:
        if not np.allclose(check[c+'_source'],check[c+'_env'],atol=1e-10,rtol=0):raise ValueError('Environment coordinates differ')
    finite=np.isfinite(e[VARS].to_numpy(float)).all(1)
    if finite.mean()<.99:raise ValueError('Do not silently fall back to old sparse environment subset')
    report={'analysis_scope':'retrospective exploratory native-only upstream numerical analysis',
            'all_observations':len(source),'source_taxa':source.taxon_name.nunique(),
            'native_observations':len(n),'native_source_taxa':n.taxon_name.nunique(),
            'native_accepted_taxa':n.accepted_key.nunique(),'range_counts':classified.native_range_status.value_counts().to_dict(),
            'eligibility_counts':classified.eligibility_reason.value_counts().to_dict(),
            'native_measurement_observations':len(f),'climate_complete':int(finite.sum()),'climate_nodata':int((~finite).sum()),
            'files_sha256':{str(Path(p).name):sha(p) for p in [roster_path,measurement_path,evidence_path,environment_path]},
            'bootstrap':{},'images_remeasured':0,'phylogeny':'reserved for a separate sensitivity'}
    for metric in ['components','joint_hue']:
        for minimum in [5,10,20,50]:
            report['bootstrap'][f'{metric}_min{minimum}']=integration(f,TRAITS,replicates=replicates,min_n=minimum,colour_metric=metric)
    data=f.merge(e[['obs_id',*VARS]],on='obs_id',validate='one_to_one')
    report['environment']=analyse(data)
    (outdir/'native_full_analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for flag in ['roster','measurements','evidence','environment','out-dir']:p.add_argument('--'+flag,type=Path,required=True)
    p.add_argument('--replicates',type=int,default=1000);a=p.parse_args()
    r=run(a.roster,a.measurements,a.evidence,a.environment,a.out_dir,a.replicates)
    print(json.dumps({k:r[k] for k in ['all_observations','native_observations','native_measurement_observations','climate_complete']},indent=2))

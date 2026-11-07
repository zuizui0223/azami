"""Original S1 observer-disjoint geographic split-half QC (aggregate only).

Sealed user_id mapping must be recovered and authenticated separately inside
protected Actions, then matched one-to-one to the 313350 frozen S1 obs IDs.
Each global observer is assigned exclusively to one half per fixed seed. At most
one original photo per observer x taxon x grid cell enters either half.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_photo_split_half_reliability import (
    EXPECTED_N, SEEDS, GRAINS, KEYS, MIN_HALF, MIN_CELLS, MIN_PHOTOS,
    _unit, _evaluate_pairs, sha256,
)

EXPECTED_OBSERVERS = 92839


def hash_side(user: str, seed: int) -> int:
    return hashlib.sha256(f'{seed}\0{user}'.encode('utf-8')).digest()[-1] & 1


def pairs_by_observer(good: pd.DataFrame, cols: list[str], *, seed: int,
                      hue: bool) -> pd.DataFrame:
    """One photo per observer/taxon-cell; user never appears in both sides."""
    d = good.copy()
    users = d.user_id.drop_duplicates()
    assignment = {u: hash_side(u, seed) for u in users}
    d['half'] = d.user_id.map(assignment).to_numpy(np.int8)
    # Caller selected one hashed photo per user x taxon-cell before assigning.
    grp = d.groupby(KEYS + ['half'], sort=False, observed=True)
    med = grp[cols].median()
    med['n'] = grp.size()
    halves = set(d.half.unique())
    if halves != {0, 1}:
        return pd.DataFrame()
    a = med.xs(0, level='half').add_suffix('_a')
    b = med.xs(1, level='half').add_suffix('_b')
    paired = a.join(b, how='inner').reset_index()
    paired = paired.loc[(paired.n_a >= MIN_HALF) & (paired.n_b >= MIN_HALF)].copy()
    if hue and len(paired):
        norm_a = np.linalg.norm(paired[[c+'_a' for c in cols]].to_numpy(float), axis=1)
        norm_b = np.linalg.norm(paired[[c+'_b' for c in cols]].to_numpy(float), axis=1)
        paired = paired.loc[(norm_a >= .20) & (norm_b >= .20)].copy()
    if paired.empty:
        return paired
    tax = paired.groupby('taxon', observed=True).agg(
        n_cells=('n_a','size'), n_a=('n_a','sum'), n_b=('n_b','sum'))
    keep = tax.index[(tax.n_cells >= MIN_CELLS) &
                     (tax.n_a + tax.n_b >= MIN_PHOTOS)]
    return paired.loc[paired.taxon.isin(keep)].copy()


def audit(frame: pd.DataFrame, construct: str, *, degrees: int,
          seeds: tuple[int, ...] = SEEDS) -> dict:
    if construct not in FEATURES or degrees not in GRAINS:
        raise ValueError('INVALID_REGISTERED_OBSERVER_SPLIT_DESIGN')
    cols=list(FEATURES[construct]); hue=construct=='hue'
    needed={'obs_id','accepted_key','user_id','latitude','longitude',*cols,
            *(c+'__eligible' for c in cols)}
    if not needed.issubset(frame):
        raise ValueError('MISSING_OBSERVER_SPLIT_SOURCE_COLUMNS')
    lat=pd.to_numeric(frame.latitude,errors='coerce').to_numpy(float)
    lon=pd.to_numeric(frame.longitude,errors='coerce').to_numpy(float)
    y=frame[cols].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    valid=(np.isfinite(y).all(axis=1)&np.isfinite(lat)&np.isfinite(lon)&
           (np.abs(lat)<=90)&(np.abs(lon)<=180)&
           frame.accepted_key.notna().to_numpy()&frame.user_id.notna().to_numpy())
    for c in cols:valid &= frame[c+'__eligible'].eq(True).to_numpy(bool)
    if hue:valid &= np.linalg.norm(y,axis=1)>=.20
    good=frame.loc[valid,['obs_id','accepted_key','user_id']].copy()
    good.rename(columns={'accepted_key':'taxon'},inplace=True)
    good['taxon']=good.taxon.astype(str)
    good['user_id']=good.user_id.astype(str)
    good['lat_cell']=np.floor((lat[valid]+90)/degrees).astype(int)
    good['lon_cell']=np.floor((lon[valid]+180)/degrees).astype(int)
    v=_unit(y[valid]) if hue else y[valid]
    for i,c in enumerate(cols):good[c]=v[:,i]
    # Stable trait-blind choice of one original observation per observer/cell.
    good['priority']=np.fromiter(
        (int.from_bytes(hashlib.sha256(s.encode('utf-8')).digest()[:8],'big')
         for s in good.obs_id),dtype=np.uint64,count=len(good))
    good=good.sort_values(['priority','obs_id'],kind='mergesort').drop_duplicates(
        KEYS+['user_id'],keep='first')
    selected_photos=int(len(good))
    reports=[]
    for seed in seeds:
        pairs=pairs_by_observer(good,cols,seed=seed,hue=hue)
        reports.append({'seed':int(seed),**_evaluate_pairs(pairs,cols,hue=hue,seed=seed)})
    return {'construct':construct,'degrees':int(degrees),
            'n_original_photo_QC_eligible':int(valid.sum()),
            'n_after_one_photo_per_observer_taxon_cell':selected_photos,
            'splits':reports}


def run(native:Path,receipt_path:Path,roster:Path,output:Path,*,
        expected_n:int=EXPECTED_N,expected_observers:int=EXPECTED_OBSERVERS,
        constructs:tuple[str,...]|None=None,grains:tuple[int,...]=GRAINS,
        seeds:tuple[int,...]=SEEDS)->dict:
    receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
    source_sha=sha256(native)
    if source_sha!=receipt.get('prepared_table_sha256') or \
       receipt.get('native_observations')!=expected_n or \
       receipt.get('environment_decoding_status')!='physical_units_gsp_invalid_encoding_masked_v1':
        raise ValueError('SOURCE_S1_DECODE_AND_SHA_NOT_AUTHENTICATED')
    chosen=tuple(FEATURES) if constructs is None else tuple(constructs)
    if (not chosen or len(chosen)!=len(set(chosen)) or not set(chosen)<=set(FEATURES) or
        not grains or len(grains)!=len(set(grains)) or not set(grains)<=set(GRAINS) or
        not seeds or len(seeds)!=len(set(seeds))):
        raise ValueError('INVALID_OBSERVER_SPLIT_SELECTION')
    columns=['obs_id','accepted_key','latitude','longitude']
    for c in chosen:
        for v in FEATURES[c]:columns.extend([v,v+'__eligible'])
    raw=pd.read_csv(native,usecols=list(dict.fromkeys(columns)),
                    dtype={'obs_id':str},low_memory=False)
    obs=pd.read_csv(roster,usecols=['obs_id','user_id'],dtype=str,low_memory=False)
    if (len(raw)!=expected_n or raw.obs_id.isna().any() or not raw.obs_id.is_unique or
        len(obs)!=expected_n or obs.obs_id.isna().any() or not obs.obs_id.is_unique or
        obs.user_id.isna().any() or obs.user_id.nunique()!=expected_observers):
        raise ValueError('ORIGINAL_OBSERVER_ROSTER_COUNT_OR_IDENTITY_DRIFT')
    merged=raw.merge(obs,on='obs_id',validate='one_to_one',how='left',sort=False)
    if len(merged)!=expected_n or merged.user_id.isna().any() or not merged.obs_id.equals(raw.obs_id):
        raise ValueError('ORIGINAL_SOURCE_OBSERVER_JOIN_FAILURE')
    results=[audit(merged,c,degrees=degree,seeds=seeds)
             for c in chosen for degree in grains]
    out={
        'schema':'azami_original_S1_global_observer_disjoint_split_half_v1',
        'status':'authenticated_exploratory_photographer_QC_not_RF_or_GEB',
        'source_original_S1_observations':int(expected_n),
        'distinct_authenticated_source_observers':int(expected_observers),
        'source_decoded_sha256':source_sha,
        'private_observer_roster_sha256':sha256(roster),
        'design':{'constructs':list(chosen),'degrees':list(grains),'seeds':list(seeds),
                  'partition_key':'global_observer_identity_SHA256_not_photo_ID',
                  'one_observation_per_observer_taxon_cell':True,
                  'min_independent_observers_per_half_cell':MIN_HALF,
                  'min_cells_per_taxon':MIN_CELLS,'min_photos_per_taxon':MIN_PHOTOS},
        'interpretation':'photo phenotype geography recurring between independent observer groups; camera/season/geography confounding and true biology remain unresolved',
        'not_an_independent_biological_individual_or_camera_calibration':True,
        'different_support_from_photo_split_and_RF_environment_complete_cases':True,
        'frozen_GEB_S1_S2_H3c_unchanged':True,
        'private_observer_ids_and_photo_rows_not_exported':True,
        'results':results}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('native','receipt','roster','output'):
        p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args()
    result=run(a.native,a.receipt,a.roster,a.output)
    print({'status':result['status'],'result_conditions':len(result['results'])},flush=True)

if __name__=='__main__':main()

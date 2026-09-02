#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

CLIMATE=['chelsa_bio01','chelsa_bio04','chelsa_bio12','chelsa_bio15']
TOPO=['topo_elevation','topo_slope','topo_roughness']
SOIL_PROPS=['bdod','cec','cfvo','clay','sand','silt','nitrogen','phh2o','soc','ocd']
DEPTHS=[('0_5cm',5.0),('5_15cm',10.0),('15_30cm',15.0)]

def parse_args():
    p=argparse.ArgumentParser(); p.add_argument('--observation',required=True); p.add_argument('--layers-root',required=True); p.add_argument('--out-dir',required=True); p.add_argument('--min-raw-layers',type=int,default=6); return p.parse_args()

def main():
    a=parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    base=pd.read_csv(a.observation,dtype={'obs_id':str},low_memory=False)
    if base.obs_id.duplicated().any(): raise SystemExit('Duplicate obs_id in source')
    root=Path(a.layers_root); coverage=[]; missing=[]
    expected=CLIMATE+TOPO+[f'soil_{p}_{d}' for p in SOIL_PROPS for d,_ in DEPTHS]
    for pid in expected:
        hits=list(root.rglob(f'{pid}.csv'))
        if len(hits)==0:
            missing.append(pid); coverage.append({'predictor':pid,'n_nonmissing':0,'complete_fraction':0.0,'role':'missing_raw_layer'}); continue
        if len(hits)!=1: raise SystemExit(f'{pid}: expected at most one layer CSV, found {len(hits)}')
        layer=pd.read_csv(hits[0],dtype={'obs_id':str})
        if layer.obs_id.duplicated().any() or set(layer.obs_id)!=set(base.obs_id): raise SystemExit(f'{pid}: obs_id mismatch')
        vals=pd.to_numeric(layer.set_index('obs_id').loc[base.obs_id,'value'],errors='coerce').to_numpy(float)
        base[pid]=vals
        coverage.append({'predictor':pid,'n_nonmissing':int(np.isfinite(vals).sum()),'complete_fraction':float(np.isfinite(vals).mean()),'role':'raw_layer'})
    merged=[p for p in expected if p in base.columns]
    if len(merged)<a.min_raw_layers: raise SystemExit(f'Only {len(merged)} raw layers available; minimum is {a.min_raw_layers}')
    for prop in SOIL_PROPS:
        numer=np.zeros(len(base)); denom=np.zeros(len(base)); used=[]
        for depth,w in DEPTHS:
            col=f'soil_{prop}_{depth}'
            if col not in base.columns: continue
            v=pd.to_numeric(base[col],errors='coerce').to_numpy(float); ok=np.isfinite(v); numer[ok]+=v[ok]*w; denom[ok]+=w; used.append(col)
        if not used: continue
        comp=np.full(len(base),np.nan); ok=denom>0; comp[ok]=numer[ok]/denom[ok]
        col=f'soil_{prop}_0_30cm'; base[col]=comp
        coverage.append({'predictor':col,'n_nonmissing':int(np.isfinite(comp).sum()),'complete_fraction':float(np.isfinite(comp).mean()),'role':'partial_depth_weighted_0_30cm','source_layers':'|'.join(used)})
    base.to_csv(out/'strict_spatial_thinned_with_climate_topography_soil.csv',index=False,encoding='utf-8-sig')
    cov=pd.DataFrame(coverage); cov.to_csv(out/'environment_predictor_coverage.csv',index=False,encoding='utf-8-sig')
    report={'observations':len(base),'species':int(base.taxon_name.nunique()),'expected_raw_layers':len(expected),'raw_layers_merged':len(merged),'missing_raw_layers':missing,'available_composites':[c for c in base.columns if c.startswith('soil_') and c.endswith('_0_30cm')]}
    (out/'parallel_environment_merge_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__': main()

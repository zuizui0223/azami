#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--observation',required=True)
    p.add_argument('--predictor-id',required=True)
    p.add_argument('--url',required=True)
    p.add_argument('--out-dir',required=True)
    p.add_argument('--sample-batch-size',type=int,default=1024)
    return p.parse_args()


def main():
    a=parse_args(); out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(a.observation,usecols=['obs_id','latitude','longitude'],dtype={'obs_id':str},low_memory=False)
    if df['obs_id'].duplicated().any(): raise SystemExit('Duplicate obs_id')
    lat=pd.to_numeric(df.latitude,errors='coerce').to_numpy(float)
    lon=pd.to_numeric(df.longitude,errors='coerce').to_numpy(float)
    if not (np.isfinite(lat)&np.isfinite(lon)).all(): raise SystemExit('Invalid coordinates')
    env=rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',CPL_VSIL_CURL_ALLOWED_EXTENSIONS='.tif,.tiff,.vrt',GDAL_HTTP_MULTIRANGE='YES',GDAL_HTTP_TIMEOUT='300',GDAL_HTTP_MAX_RETRY='8',GDAL_HTTP_RETRY_DELAY='5',VSI_CACHE='TRUE',VSI_CACHE_SIZE='200000000')
    env.__enter__(); src=None
    try:
        src=rasterio.open(a.url if a.url.startswith('/vsicurl/') else '/vsicurl/'+a.url)
        tr=Transformer.from_crs('EPSG:4326',src.crs,always_xy=True)
        xs,ys=tr.transform(lon.tolist(),lat.tolist())
        coords=list(zip(xs,ys)); values=np.full(len(coords),np.nan,float)
        for start in range(0,len(coords),a.sample_batch_size):
            stop=min(start+a.sample_batch_size,len(coords))
            for j,sample in enumerate(src.sample(coords[start:stop],indexes=1,masked=True)):
                v=sample[0]
                if np.ma.is_masked(v): continue
                try: v=float(v)
                except Exception: continue
                if math.isfinite(v): values[start+j]=v
        coverage=float(np.isfinite(values).mean())
        pd.DataFrame({'obs_id':df.obs_id.astype(str),'predictor_id':a.predictor_id,'value':values}).to_csv(out/f'{a.predictor_id}.csv',index=False,encoding='utf-8-sig')
        meta={'predictor_id':a.predictor_id,'url':a.url,'coverage':coverage,'n_observations':len(df),'crs':str(src.crs),'dtype':str(src.dtypes[0]),'nodata':None if src.nodata is None else float(src.nodata),'scales':list(src.scales),'offsets':list(src.offsets)}
        (out/f'{a.predictor_id}.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps(meta,ensure_ascii=False))
    finally:
        if src is not None: src.close()
        env.__exit__(None,None,None)

if __name__=='__main__': main()

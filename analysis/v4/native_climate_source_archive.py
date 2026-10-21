"""Lossless CHELSA cache for the native-only roster's 5-degree occupied cells.

Five degrees controls file grouping ONLY. Pixels retain the original ~1 km grid,
values, nodata and scale/offset. No coarse environmental averaging or imputation.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import numpy as np

VARIABLES=('bio01','bio04','bio12','bio15','rsds_mean','vpd_mean','sfcwind_mean','gsp','npp')


def choose_url_paths(variable):
    if variable not in VARIABLES:raise ValueError('Unregistered climate variable')
    label='sfcWind_mean' if variable=='sfcwind_mean' else variable
    filename=f'CHELSA_{label}_1981-2010_V.2.1.tif'
    if variable=='rsds_mean':filename='CHELSA_rsds_1981-2010_mean_V.2.1.tif'
    old='https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/'
    new='https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/'
    if variable.startswith('bio'):
        return [new+f'bioclim/{variable}/1981-2010/'+filename,old+filename]
    return [old+filename,new+'climatologies/1981-2010/bio/'+filename,
            'https://os.unil.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/'+filename]


def tile_window(src,ilat,ilon):
    from rasterio.windows import Window,from_bounds
    south=-90+5*ilat;west=-180+5*ilon
    w=from_bounds(west,south,west+5,south+5,src.transform)
    left=int(np.floor(w.col_off+1e-7));top=int(np.floor(w.row_off+1e-7))
    right=int(np.ceil(w.col_off+w.width-1e-7));bottom=int(np.ceil(w.row_off+w.height-1e-7))
    return Window(left,top,right-left,bottom-top).intersection(Window(0,0,src.width,src.height))


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def extract_tiles(source,cells,outdir):
    import rasterio
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    with rasterio.open(source) as src:
        if src.crs is None or src.crs.to_epsg()!=4326:raise ValueError('Expected geographic CHELSA grid')
        meta={'transform':list(src.transform)[:6],'crs':str(src.crs),'width':src.width,'height':src.height,
              'nodata':src.nodata,'dtype':src.dtypes[0],'scale':src.scales[0],'offset':src.offsets[0],
              'tags':src.tags(),'source_sha256':digest(source),'tiles':{},'resampling':'none'}
        for start in range(0,len(cells),32):
            arrays={};name=f'tiles-{start//32:03d}.npz'
            for a,b in cells[start:start+32]:
                key=f'{a}_{b}';w=tile_window(src,int(a),int(b));data=src.read(1,window=w)
                arrays[key]=data
                meta['tiles'][key]={'file':name,'row_off':int(w.row_off),'col_off':int(w.col_off),
                                    'height':int(w.height),'width':int(w.width)}
            np.savez_compressed(outdir/name,**arrays)
    meta['cache_sha256']={p.name:digest(p) for p in outdir.glob('tiles-*.npz')}
    (outdir/'grid.json').write_text(json.dumps(meta,indent=2))
    return meta


def sample_tiles(directory,lat,lon):
    from affine import Affine
    from rasterio.transform import rowcol
    d=Path(directory);m=json.loads((d/'grid.json').read_text())
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    valid=np.isfinite(lat)&np.isfinite(lon)&(abs(lat)<=90)&(abs(lon)<=180)
    out=np.full(len(lat),np.nan);ids=np.flatnonzero(valid)
    rr,cc=rowcol(Affine(*m['transform']),lon[valid],lat[valid]);rr=np.asarray(rr);cc=np.asarray(cc)
    a=np.floor((lat[valid]+90)/5).astype(int);b=np.floor((lon[valid]+180)/5).astype(int)
    keys=np.array([f'{x}_{y}' for x,y in zip(a,b)])
    loaded={}
    for key in np.unique(keys):
        if key not in m['tiles']:continue
        t=m['tiles'][key];fn=t['file']
        if fn not in loaded:loaded[fn]=np.load(d/fn)
        arr=loaded[fn][key];sel=np.flatnonzero(keys==key);r=rr[sel]-t['row_off'];c=cc[sel]-t['col_off']
        good=(r>=0)&(r<arr.shape[0])&(c>=0)&(c<arr.shape[1]);sel=sel[good];v=arr[r[good],c[good]].astype(float)
        if m['nodata'] is not None:v[v==m['nodata']]=np.nan
        out[ids[sel]]=v
    for z in loaded.values():z.close()
    return out


def fetch(variable,cells_path,outdir):
    import requests
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    temporary=outdir/'source.tif';errors=[];receipt=None
    for url in choose_url_paths(variable):
        try:
            with requests.get(url,stream=True,timeout=(20,90),headers={'User-Agent':'azami-native-upstream/20261020'}) as r:
                r.raise_for_status();size=0;started=time.monotonic()
                with temporary.open('wb') as f:
                    for part in r.iter_content(1<<20):
                        f.write(part);size+=len(part)
                        if size>5_000_000_000:raise ValueError('Unexpectedly large source')
                        if time.monotonic()-started>600:raise TimeoutError('Source transfer exceeded budget')
                receipt={'url':url,'etag':r.headers.get('ETag'),'last_modified':r.headers.get('Last-Modified'),'bytes':size}
            import rasterio
            with rasterio.open(temporary) as src:
                if src.count!=1:raise ValueError('Expected a single climate band')
            break
        except Exception as exc:
            errors.append({'url':url,'error':str(exc)});temporary.unlink(missing_ok=True);receipt=None
    if receipt is None:
        (outdir/'acquisition_failure.json').write_text(json.dumps(errors,indent=2))
        raise RuntimeError('No approved CHELSA route succeeded')
    contract=json.loads(Path(cells_path).read_text());meta=extract_tiles(temporary,contract['cells'],outdir)
    receipt.update(variable=variable,completed_utc=datetime.now(timezone.utc).isoformat(),
                   cells=len(contract['cells']),source_sha256=meta['source_sha256'],failed_routes=errors,
                   resolution_preserved=True,native_target_observations=contract['native_target_observations'])
    (outdir/'source_receipt.json').write_text(json.dumps(receipt,indent=2));temporary.unlink()
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--variable',choices=VARIABLES,required=True)
    p.add_argument('--cells',type=Path,required=True);p.add_argument('--out-dir',type=Path,required=True)
    args=p.parse_args();fetch(args.variable,args.cells,args.out_dir)

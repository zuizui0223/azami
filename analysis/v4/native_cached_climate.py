"""Read ONLY the five already cached original-pixel CHELSA v2.1 rasters.

No network, downloading, 5-degree averaging, interpolating missing introduced
regions or fitting responses here. Climate archive SHA-256 and per-tile bytes
are locked to the previous original CHELSA provenance. 5-degree tiles are
only storage groups; pixel lookup uses the original affine transform.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np

VARIABLES=("chelsa_bio01","chelsa_bio04","chelsa_bio12",
           "chelsa_bio15","chelsa_rsds_mean")


def digest_file(path:Path)->str:
    sha=hashlib.sha256()
    with path.open("rb") as f:
        for blob in iter(lambda:f.read(1<<22),b""):
            sha.update(blob)
    return sha.hexdigest()


def verified_cache(directory:Path,contract:dict)->dict:
    root=Path(directory)
    with (root/"grid.json").open() as f:
        grid=json.load(f)
    if grid.get("source_sha256")!=contract.get("source_sha256"):
        raise ValueError("Original CHELSA source changed")
    if grid.get("resampling")!="none":
        raise ValueError("CHELSA original pixel values were resampled")
    if len(grid["tiles"])!=316:
        raise ValueError("Original global cache cell census changed")
    for key in ("scale","offset","nodata","dtype","width","height"):
        if str(grid[key])!=str(contract[key]):
            # 0 versus 0.0 are equivalent floating values and harmless.
            if key not in ("scale","offset","nodata") or float(grid[key])!=float(contract[key]):
                raise ValueError("Original CHELSA raster contract changed: "+key)
    if not np.allclose(grid["transform"],contract["transform"],atol=1e-11,rtol=0):
        raise ValueError("Original CHELSA affine transform changed")
    for member,expected in grid["cache_sha256"].items():
        if not member.startswith("tiles-") or not member.endswith(".npz"):
            raise ValueError("Unexpected original raster tile member")
        if digest_file(root/member)!=expected:
            raise ValueError("Cached CHELSA tile checksum changed: "+member)
    return grid


def sample_cached_pixels(directory:Path,grid:dict,
                         latitude,longitude)->np.ndarray:
    lat=np.asarray(latitude,dtype=float)
    lon=np.asarray(longitude,dtype=float)
    if lat.shape!=lon.shape:
        raise ValueError("Latitude/longitude coordinate mismatch")
    valid=np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)&(np.abs(lon)<=180)
    out=np.full(len(lat),np.nan,dtype=float)
    a,b,c,d,e,f=[float(x) for x in grid["transform"]]
    # Invert the full affine grid transform, not a 5-degree site proxy.
    det=a*e-b*d
    if not np.isfinite(det) or abs(det)<1e-16:
        raise ValueError("Singular source raster transform")
    ids=np.flatnonzero(valid)
    la=lat[valid];lo=lon[valid]
    column=np.floor((e*(lo-c)-b*(la-f))/det).astype(int)
    row=np.floor((-d*(lo-c)+a*(la-f))/det).astype(int)
    tile_lat=np.floor((la+90)/5).astype(int)
    tile_lon=np.floor((lo+180)/5).astype(int)
    label=np.array([str(i)+"_"+str(j) for i,j in zip(tile_lat,tile_lon)])
    loaded={}
    try:
        for key in np.unique(label):
            if key not in grid["tiles"]:
                continue
            tile=grid["tiles"][key]
            name=tile["file"]
            if name not in loaded:
                loaded[name]=np.load(directory/name)
            indices=np.flatnonzero(label==key)
            image=loaded[name][key]
            yy=row[indices]-int(tile["row_off"])
            xx=column[indices]-int(tile["col_off"])
            ok=(yy>=0)&(xx>=0)&(yy<image.shape[0])&(xx<image.shape[1])
            values=image[yy[ok],xx[ok]].astype(float)
            if grid["nodata"] is not None:
                values[values==float(grid["nodata"])]=np.nan
            out[ids[indices[ok]]]=values
    finally:
        for cached in loaded.values():
            cached.close()
    return out


def physical_climate_for_focal(frame,root:Path,contract_path:Path)->tuple[dict,dict]:
    with Path(contract_path).open() as f:
        source=json.load(f)
    if not all(variable in source for variable in VARIABLES):
        raise ValueError("Five previously selected climatic variables missing from fixed metadata")
    if len(frame)<1000:
        raise ValueError("Original focal matched sample is unexpectedly small")
    if frame.obs_id.duplicated().any():
        raise ValueError("Original focal source observations duplicated")
    if frame.latitude.isna().any() or frame.longitude.isna().any():
        raise ValueError("Wild focal source record without usable coordinate")
    out={}
    audit={}
    for variable in VARIABLES:
        meta=source[variable]
        directory=Path(root)/variable
        grid=verified_cache(directory,meta)
        unscaled=sample_cached_pixels(directory,grid,
                                      frame.latitude.to_numpy(),
                                      frame.longitude.to_numpy())
        final=unscaled*float(grid["scale"])+float(grid["offset"])
        audit[variable]={
            "source_raster_sha256":str(grid["source_sha256"]),
            "grid_cached_cells":len(grid["tiles"]),
            "pixels_resampled":False,
            "n_finite":int(np.isfinite(final).sum()),
            "n_missing":int(np.isnan(final).sum()),
            "original_scale":grid["scale"],
            "original_offset":grid["offset"],
            "missing_tiles_imputed":False,
        }
        out[variable]=final
    return out,audit

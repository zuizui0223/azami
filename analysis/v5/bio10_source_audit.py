"""Validate official CHELSA v2.1 BIO10 source before using it in V5.

This is deliberately a *source-only* audit: no individual observations, heads,
raster pixels or outcome data are used. It emits a portable file hash and
georeferencing receipt, to be pinned before the comparative biological fit.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import requests
import rasterio
from analysis.v4.s1_s2_extreme_environment import CHELSA_TEMPLATE, CHELSA_LEGACY

OUT=Path("work/submission/public")
CACHE=Path("work/submission/private")
URLS=(CHELSA_TEMPLATE.format(bio="bio10"),CHELSA_LEGACY.format(bio="bio10"))

def acquire():
    CACHE.mkdir(parents=True,exist_ok=True)
    OUT.mkdir(parents=True,exist_ok=True)
    target=CACHE/"CHELSA_bio10_1981-2010_V.2.1.tif"
    errors=[]
    for url in URLS:
        target.unlink(missing_ok=True)
        sha=hashlib.sha256();n=0
        try:
            with requests.get(url,stream=True,timeout=(30,180),
                              headers={"User-Agent":"azami-bio10-source-audit/20261110"}) as res:
                res.raise_for_status()
                with target.open("wb") as h:
                    for chunk in res.iter_content(1<<22):
                        if not chunk:continue
                        n+=len(chunk)
                        if n>2_500_000_000:raise ValueError("BIO10_SOURCE_OVERSIZE")
                        sha.update(chunk)
                        h.write(chunk)
            if n<50_000_000:raise ValueError("BIO10_SOURCE_TOO_SMALL")
            with rasterio.open(target) as ds:
                metadata=dict(
                    width=ds.width,height=ds.height,
                    epsg=ds.crs.to_epsg(),count=ds.count,
                    dtype=ds.dtypes[0],
                    nodata=None if ds.nodata is None else float(ds.nodata),
                    scale=float(ds.scales[0]),offset=float(ds.offsets[0]),
                    transform=[float(q) for q in tuple(ds.transform)[:6]]
                )
            if (metadata["epsg"]!=4326 or metadata["width"]!=43200
                or metadata["height"]!=20880 or metadata["count"]!=1
                or metadata["dtype"]!="uint16"
                or metadata["nodata"]!=65535
                or abs(metadata["scale"]-.1)>1e-10
                or abs(metadata["offset"]+273.15)>1e-10):
                raise ValueError("BIO10_OFFICIAL_RASTER_GEOMETRY_OR_UNITS_UNEXPECTED")
            receipt=dict(
                variable="chelsa_bio10",
                interpretation="mean temperature of warmest quarter",
                source="CHELSA V2.1 1981-2010 bio10 original TIFF",
                url=url,sha256=sha.hexdigest(),size_bytes=n,
                **metadata,
                original_pixels_only=True,
                source_verified_no_biological_samples=True)
            (OUT/"v5_bio10_source_receipt.json").write_text(
                json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
            print(json.dumps(receipt,ensure_ascii=False),flush=True)
            return receipt
        except (requests.RequestException,ValueError,RuntimeError,OSError) as exc:
            errors.append({"url":url,"reason":str(exc)[:200]})
    raise RuntimeError("NO_AUTHENTIC_CHELSA_BIO10_SOURCE "+json.dumps(errors))

if __name__=="__main__":
    try:acquire()
    finally:
        (CACHE/"CHELSA_bio10_1981-2010_V.2.1.tif").unlink(missing_ok=True)

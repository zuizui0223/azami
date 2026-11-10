"""Shared source audit and original CHELSA BIO05/BIO18 pixel sampler for V5 I+II."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from analysis.v4.s1_s2_extreme_environment import urls
from analysis.v5.native_pca_env_fullcoverage import sample_original_pixels

BASE=Path("work/submission")
PRIVATE=BASE/"private"
OUT=BASE/"public"
RECEIPT=Path("manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/sampling_coverage/environment_source_receipt.json")
EXTRA=("chelsa_bio18","chelsa_bio05")


def sample_warm(wild, land, native):
    """Sample the *same original TIFF pixels* for all three WCVP/S1 units."""
    sources=json.loads(RECEIPT.read_text())["sources"]
    known={s["variable"]:s for s in sources}
    if set(EXTRA)-set(known):
        raise ValueError("NO_PINNED_BIO05_BIO18_SOURCE_RECEIPT")
    n_w,n_l,n_s=len(wild),len(land),len(native)
    assert n_w==297622 and n_s==313350 and n_l>20000
    sites=pd.concat(
        [wild[["latitude","longitude"]],
         land[["latitude","longitude"]],
         native[["latitude","longitude"]]],
        ignore_index=True)
    out={}
    for var in EXTRA:
        spec=known[var]
        if not spec["sampled_from_original_raster_pixels"] or spec["crs"]!="EPSG:4326":
            raise ValueError("WARM_SOURCE_PROVENANCE_CHANGED")
        original_urls=urls(var)
        target=PRIVATE/(var+"_v5_pinned_original.tif")
        attempts=[]
        ok=False
        for url in original_urls:
            import requests
            target.unlink(missing_ok=True)
            sha=hashlib.sha256()
            size=0
            try:
                with requests.get(url,stream=True,timeout=(25,120),
                    headers={"User-Agent":"azami-V5-pinned-annual-warm/20261110"}) as res:
                    res.raise_for_status()
                    with target.open("wb") as handle:
                        for chunk in res.iter_content(1<<22):
                            if not chunk:continue
                            size+=len(chunk)
                            if size>3_500_000_000:
                                raise RuntimeError("SOURCE_RASTER_OVERSIZE")
                            sha.update(chunk);handle.write(chunk)
                if sha.hexdigest()!=spec["download_sha256"] or size!=spec["download_bytes"]:
                    raise ValueError("WARM_BIO_SOURCE_SHA_OR_SIZE_MISMATCH")
                with rasterio.open(target) as src:
                    if (src.crs.to_epsg()!=4326 or src.count!=1 or
                        src.width!=spec["width"] or src.height!=spec["height"] or
                        src.dtypes[0]!=spec["source_dtype"] or
                        src.nodata!=spec["source_nodata"] or
                        abs(src.scales[0]-spec["scale"])>1e-10 or
                        abs(src.offsets[0]-spec["offset"])>1e-10):
                        raise ValueError("PINNED_WARM_RASTER_LAYOUT_CHANGED")
                ok=True;break
            except (requests.RequestException, ValueError, RuntimeError, OSError) as exc:
                attempts.append({"reason":type(exc).__name__,"message":str(exc)[:200]})
        if not ok:
            target.unlink(missing_ok=True)
            raise ValueError("ORIGINAL_WARM_RASTER_UNAVAILABLE_OR_MISMATCH_"+var+json.dumps(attempts))
        try:
            z=sample_original_pixels(target,sites)
            z=z*float(spec["scale"])+float(spec["offset"])
        finally:
            target.unlink(missing_ok=True)
        if len(z)!=(n_w+n_l+n_s):
            raise ValueError("ORIGINAL_WARM_PIXEL_SAMPLER_LENGTH_DRIFT")
        # V4 verified all source 469,958 unique locations; full original
        # CHELSA is allowed to lack coverage, but cohort parity requires
        # the exact intersection of available original observations.
        for frame, part in ((wild,z[:n_w]),(land,z[n_w:n_w+n_l]),(native,z[n_w+n_l:])):
            frame[var]=part
        out[var]={
            "source_sha256":sha.hexdigest(),
            "download_bytes":size,
            "wild_finite":int(np.isfinite(z[:n_w]).sum()),
            "background_finite":int(np.isfinite(z[n_w:n_w+n_l]).sum()),
            "native_s1_finite":int(np.isfinite(z[n_w+n_l:]).sum()),
            "original_raster_pixels_only":True,
            "raster_deleted_after_sampling":True,
        }
    return wild,land,native,out

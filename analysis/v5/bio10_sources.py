"""Original CHELSA BIO10 mean-of-warmest-quarter climate for V5 I and II.

V5 original BIO05 and BIO18 originals are checksum pinned separately.
This new BIO10 source is available at the same official CHELSA 1981–2010
bioclim archive; a source-only audit records its hash and georeferencing.
Before considering biological estimates final, pin the emitted SHA256 and
rerun on the exact same source. See PINNED_BIO10_SHA256 below.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from analysis.v5.bio10_source_audit import acquire, CACHE
from analysis.v5.native_pca_env_fullcoverage import sample_original_pixels

BIO10="chelsa_bio10"
PINNED_BIO10_SHA256=None  # Must be set after independently checked provenance-only run.


def sample_bio10(wild, land, native):
    receipt=acquire()
    if receipt["variable"] != BIO10 or not receipt["original_pixels_only"]:
        raise ValueError("BIO10_SOURCE_AUDIT_INVALID")
    if PINNED_BIO10_SHA256 is not None and receipt["sha256"] != PINNED_BIO10_SHA256:
        raise ValueError("BIO10_SOURCE_HASH_CHANGED_FROM_PIN")
    source=CACHE/"CHELSA_bio10_1981-2010_V.2.1.tif"
    nw,nl,ns=len(wild),len(land),len(native)
    if nw!=297622 or ns!=313350 or nl<20000:
        raise ValueError("V5_BIO10_ORIGINAL_SOURCE_GRAIN_CHANGED")
    combined=pd.concat([
        wild[["latitude","longitude"]],
        land[["latitude","longitude"]],
        native[["latitude","longitude"]],
    ],ignore_index=True)
    try:
        raw=sample_original_pixels(source,combined)
        physical=raw*receipt["scale"]+receipt["offset"]
        if len(physical)!=nw+nl+ns:
            raise ValueError("BIO10_ORIGINAL_RASTER_SAMPLING_LENGTH_DRIFT")
        spans=(physical[:nw],physical[nw:nw+nl],physical[nw+nl:])
        for frame,values in zip((wild,land,native),spans):
            frame[BIO10]=values
        fraction=[float(np.isfinite(v).mean()) for v in spans]
        if fraction[0]<.98 or fraction[1]<.95 or fraction[2]<.98:
            raise ValueError("BIO10_V5_PHOTO_OR_BACKGROUND_COVERAGE_TOO_LOW")
        return wild,land,native,{
            "sha256":receipt["sha256"],
            "pinned_sha256":PINNED_BIO10_SHA256,
            "sha256_checked_against_pin":PINNED_BIO10_SHA256 is not None,
            "url":receipt["url"],
            "scale":receipt["scale"],"offset":receipt["offset"],
            "source_size_bytes":receipt["size_bytes"],
            "n_wild_finite":int(np.isfinite(spans[0]).sum()),
            "n_land_finite":int(np.isfinite(spans[1]).sum()),
            "n_native_finite":int(np.isfinite(spans[2]).sum()),
            "original_pixels_only":True,
            "no_interpolation":True,
        }
    finally:
        source.unlink(missing_ok=True)

"""Rehydrate V4's already-authenticated five original-pixel CHELSA S2 grids.

Reuses exactly the five original raster caches, archived SHA256 and per-tile
integrity checks, for 297,622 source-authenticated wild native/introduced
observations. This is not an additional climate hypothesis or a new download
from uncontrolled raster providers. No individual records leave this runner.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import subprocess
import zipfile
import numpy as np
from analysis.v4.native_cached_climate import physical_climate_for_focal

CACHES = {
    "chelsa_bio01": (12243016948, "0e3b83ebfb7de1865662cde2c2c2395fa7737460996e66354f87391f8f0867ab"),
    "chelsa_bio04": (12241949860, "dac232683b2a72b9c4f62eab0f341be7927e5e3287a29dbf3a162b893c2307e3"),
    "chelsa_bio12": (12242703088, "71e7edeb0d4c0042e973ed97672a1593bef98a608c917792006e22722e288ada"),
    "chelsa_bio15": (12242508971, "9f255a7c756fc83a5dd81e6c13d8fee63c53698d93bffb65ae13288e4eebb445"),
    "chelsa_rsds_mean": (12242917323, "2193bd04b979f6684f7bce7c670074af63072d270a20c4dc2fd0bd5934862df4"),
}


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def recover_original_five(focal,cache_root,metadata_path):
    if len(focal)!=297622:
        raise ValueError("S2_BASE_EXACT_SOURCE_ROSTER_LOST")
    if not os.environ.get("GH_TOKEN"):
        raise ValueError("S2_ORIGINAL_CHELSA_PRIVATE_ARTIFACT_TOKEN_UNAVAILABLE")
    root=Path(cache_root)
    if "private" not in root.parts:
        raise ValueError("S2_CHELSA_CACHE_MUST_REMAIN_PRIVATE")
    root.mkdir(parents=True,exist_ok=True)
    receipts={}
    for variable,(artifact,expected) in CACHES.items():
        dest=root/variable
        dest.mkdir(parents=True,exist_ok=True)
        zip_path=root/(variable+".zip")
        with zip_path.open("wb") as output:
            subprocess.run(["gh","api",
                            f"repos/zuizui0223/azami/actions/artifacts/{artifact}/zip"],
                           stdout=output,check=True)
        if sha256(zip_path)!=expected:
            zip_path.unlink(missing_ok=True)
            raise ValueError("S2_FIVE_CHELSA_SOURCE_ARCHIVE_HASH_MISMATCH_"+variable)
        with zipfile.ZipFile(zip_path) as archive:
            members=archive.namelist()
            if any(Path(n).is_absolute() or ".." in Path(n).parts for n in members):
                raise ValueError("UNSAFE_CHELSA_CACHE_MEMBER")
            archive.extractall(dest)
        zip_path.unlink()
        receipts[variable]={"archive_id":artifact,"archive_sha256":expected,
                            "files_restored":len(members)}
    values,audit=physical_climate_for_focal(focal,root,metadata_path)
    for variable,climate in values.items():
        a=np.asarray(climate,float)
        if len(a)!=len(focal):
            raise ValueError("S2_BASE_CLIMATE_SOURCE_ROW_MISMATCH")
        if float(np.isfinite(a).mean())<.95:
            raise ValueError("S2_BASE_CLIMATE_PIXEL_COVERAGE_BELOW_95PC_"+variable)
        focal[variable]=a
    return focal,{"status":"EXACT_FIVE_ORIGINAL_CHELSA_RASTERS_REUSED",
                  "variable_archives":receipts,
                  "original_pixel_audit":audit,
                  "original_individual_ids_or_coordinates_exported":False}

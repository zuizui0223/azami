"""Restore only original photographer/quality metadata in a PRIVATE Actions job.

No raw original image/source ZIP is written, published or uploaded. The existing
sealed metadata is read from the archived authenticated Release asset with the
current repository secret; no new iNaturalist/WCVP/climate/photo acquisition.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import os
import zipfile
from pathlib import Path

import pandas as pd

SEALED_SHA256="dfa98fa162a8029cb6477adf4c26ffc5819cc916a0513d32da14d78d845ab4b0"
SOURCE_ZIP_SHA256="2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae"
AAD=b"azami-original-metadata-v1"
EXPECTED_SOURCE_ROWS=1122854
EXPECTED_SOURCE_OBSERVATIONS=665115
EXPECTED_NATIVE_N=313350
EXPECTED_NATIVE_OBSERVERS=92839


def original_observer_roster(sealed:bytes, key:bytes, native_ids:pd.Series) -> pd.DataFrame:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    if hashlib.sha256(sealed).hexdigest()!=SEALED_SHA256:
        raise ValueError("Original encrypted metadata source identity changed")
    if len(key)!=32:
        raise ValueError("Original metadata key must have 32 bytes")
    try:
        clear=AESGCM(key).decrypt(sealed[:12],sealed[12:],AAD)
    except Exception:
        raise ValueError("Original source metadata AES-GCM authentication failed") from None
    if hashlib.sha256(clear).hexdigest()!=SOURCE_ZIP_SHA256:
        raise ValueError("Original metadata plaintext checksum mismatch")
    with zipfile.ZipFile(io.BytesIO(clear)) as z:
        if "photo_metadata_merged.csv" not in z.namelist():
            raise ValueError("Original metadata column source missing")
        with z.open("photo_metadata_merged.csv") as handle:
            photos=pd.read_csv(handle,usecols=["obs_id","user_id","quality_grade"],
                               dtype={"obs_id":str,"user_id":str},low_memory=False)
    if len(photos)!=EXPECTED_SOURCE_ROWS or photos.obs_id.nunique()!=EXPECTED_SOURCE_OBSERVATIONS:
        raise ValueError("Original metadata changed population")
    if photos.groupby("obs_id").user_id.nunique(dropna=False).gt(1).any():
        raise ValueError("Same source observation attributed to multiple observers")
    if photos.groupby("obs_id").quality_grade.nunique(dropna=False).gt(1).any():
        raise ValueError("Same source observation attributed to inconsistent quality grades")
    observers=photos.drop_duplicates("obs_id")[["obs_id","user_id","quality_grade"]]
    if native_ids.isna().any() or not native_ids.is_unique or len(native_ids)!=EXPECTED_NATIVE_N:
        raise ValueError("Native roster is not the frozen 313350 observations")
    mapped=pd.DataFrame({"obs_id":native_ids.astype(str)}).merge(
        observers,on="obs_id",validate="one_to_one",how="left")
    if mapped.user_id.isna().any():
        raise ValueError("Missing original observer in native subset")
    if mapped.user_id.nunique()!=EXPECTED_NATIVE_OBSERVERS:
        raise ValueError("Original native observer population changed")
    return mapped


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sealed",type=Path,required=True)
    p.add_argument("--native",type=Path,required=True)
    p.add_argument("--out-private",type=Path,required=True)
    a=p.parse_args()
    keyhex=os.environ.get("AZAMI_METADATA_KEY","")
    if len(keyhex)!=64:
        raise ValueError("AZAMI_METADATA_KEY unavailable or malformed")
    ids=pd.read_csv(a.native, usecols=["obs_id"],dtype={"obs_id":str}).obs_id
    result=original_observer_roster(a.sealed.read_bytes(),bytes.fromhex(keyhex),ids)
    a.out_private.parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(a.out_private,index=False,compression="gzip")
    print({"status":"private_observer_source_verified",
           "native_observations":len(result),
           "distinct_original_observers":result.user_id.nunique(),
           "privacy":"NO HUMAN IDs EXPORTED IN PUBLIC ARTIFACT"},flush=True)


if __name__=="__main__":
    main()

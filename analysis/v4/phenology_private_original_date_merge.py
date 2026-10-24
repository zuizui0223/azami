"""Private reconstruction of original native S1 flowering-date inputs.

The eligible flowering-month proxy is NOT true anthesis or an observation of
pollen presentation. A trusted original-source date must come from exact native
source obs_id or the authenticated release metadata; never infer from photo
upload date, photo outcome or contemporary iNaturalist observations.

Reads AES256-GCM original metadata using a GitHub Actions secret only. Raw
obs IDs, photographers, dates and lat/lon are never written to public artifacts.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import json
import os
from pathlib import Path
import zipfile

import pandas as pd

from analysis.v4.native_observers import (
    SEALED_SHA256, SOURCE_ZIP_SHA256, AAD, EXPECTED_SOURCE_ROWS,
    EXPECTED_SOURCE_OBSERVATIONS, EXPECTED_NATIVE_N,
    EXPECTED_NATIVE_OBSERVERS,
)

POSSIBLE_DATE_FIELDS=("observed_on","observed_on_string")
META_CSV="photo_metadata_merged.csv"


def original_dates(native_path:Path,sealed_path:Path,private_output:Path):
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    encoded=sealed_path.read_bytes()
    if hashlib.sha256(encoded).hexdigest()!=SEALED_SHA256:
        raise ValueError("Original encrypted iNaturalist event evidence drift")
    token=os.environ.get("AZAMI_METADATA_KEY","")
    if len(token)!=64:
        raise ValueError("Private original source decryption key unavailable")
    try:
        plaintext=AESGCM(bytes.fromhex(token)).decrypt(encoded[:12],encoded[12:],AAD)
    except Exception:
        raise ValueError("Original archived observation metadata failed GCM") from None
    if hashlib.sha256(plaintext).hexdigest()!=SOURCE_ZIP_SHA256:
        raise ValueError("Decrypted original release archive SHA does not match source")
    # Inspect only trusted manifest names, not raw individuals.
    with zipfile.ZipFile(io.BytesIO(plaintext)) as archive:
        if META_CSV not in archive.namelist():
            raise ValueError("Original source observation metadata CSV absent")
        with archive.open(META_CSV) as obj:
            header=list(pd.read_csv(obj,nrows=0).columns)
        datecols=[c for c in POSSIBLE_DATE_FIELDS if c in header]
        read_cols=["obs_id","user_id","quality_grade",*datecols]
        if not all(x in header for x in read_cols):
            raise ValueError("Original metadata lacks source IDs or original photographer")
        with archive.open(META_CSV) as obj:
            df=pd.read_csv(obj,usecols=read_cols,
                           dtype={"obs_id":str,"user_id":str},
                           low_memory=False)
    if (len(df)!=EXPECTED_SOURCE_ROWS or
        df.obs_id.nunique()!=EXPECTED_SOURCE_OBSERVATIONS):
        raise ValueError("Original release observation universe not matched")
    if df.groupby("obs_id").user_id.nunique(dropna=False).gt(1).any():
        raise ValueError("One original observation assigned multiple photographers")
    # Metadata is photo-level and must aggregate to ORIGINAL observation
    # by a deterministic, conflict-preserving equality gate.
    local=df.drop_duplicates("obs_id")
    for datecol in datecols:
        values=df[["obs_id",datecol]].copy()
        values[datecol]=values[datecol].fillna("").astype(str).str.strip()
        conflicts=values.groupby("obs_id")[datecol].nunique(dropna=False)
        if (conflicts>1).any():
            raise ValueError("Conflicting original observed dates among photos of one event")
    orig=pd.read_csv(native_path,compression="infer",
                     usecols=lambda c:c in {"obs_id","observed_on","observed_on_string"},
                     dtype={"obs_id":str},low_memory=False)
    if len(orig)!=EXPECTED_NATIVE_N or not orig.obs_id.is_unique:
        raise ValueError("S1 native observational cohort changed")
    # All original dates are independently checked; never substitute upload
    # or taxon flowering phenology from a presumed species calendar.
    ncols=[c for c in POSSIBLE_DATE_FIELDS if c in orig.columns]
    master=orig[["obs_id",*ncols]].merge(
        local[["obs_id","user_id",*datecols]],on="obs_id",how="left",
        validate="one_to_one",suffixes=("_native","_archive"))
    if master.user_id.isna().any() or master.user_id.nunique()!=EXPECTED_NATIVE_OBSERVERS:
        raise ValueError("Original S1 native observer identity missing")
    original_candidates=[c for c in master.columns
                         if c in POSSIBLE_DATE_FIELDS or c.startswith("observed_on_")]
    # Prioritise observed_on from the original S1 native rows if available.
    # Ingest only original event dates, NEVER photographer timestamps or
    # upload times as proxies. If both sources exist, disagreements fail closed.
    match=[c for c in ("observed_on_native","observed_on","observed_on_string_native",
                       "observed_on_archive","observed_on_string_archive","observed_on_string")
           if c in original_candidates]
    if not match:
        return {"status":"STOP_NO_ORIGINAL_OBSERVATION_DATE_IN_AUTHENTICATED_SOURCES",
                "native_observations":EXPECTED_NATIVE_N,
                "archive_source_observation_count":EXPECTED_SOURCE_OBSERVATIONS,
                "date_field_candidates":[],
                "source_has_no_defensible_flowering_month":True},None
    date_series=None
    matches=[]
    for col in match:
        dt=pd.to_datetime(master[col],errors="coerce",format="mixed")
        if date_series is None:
            date_series=dt
        else:
            differing=(dt.notna()&date_series.notna()&
                       (dt.dt.date!=date_series.dt.date))
            if differing.any():
                raise ValueError("Original immutable source metadata disagree on dated event")
            date_series=date_series.fillna(dt)
        matches.append(col)
    if date_series is None:
        raise ValueError("No original observation date found")
    available=(date_series.dt.year.between(2000,2026)&date_series.dt.month.between(1,12))
    merged=master[["obs_id","user_id"]].copy()
    merged["observed_on"]=date_series.dt.strftime("%Y-%m-%d")
    private_output.parent.mkdir(parents=True,exist_ok=True)
    merged.to_csv(private_output,index=False,compression="gzip")
    report={
        "status":"ORIGINAL_OBSERVATION_DATE_SOURCE_VERIFIED_PRIVATE_ONLY",
        "native_observations":len(merged),
        "source_date_columns_used":matches,
        "date_valid_native_observations":int(available.sum()),
        "native_observations_missing_valid_event_date":int((~available).sum()),
        "month_distribution_original_events":{str(m):int(((date_series.dt.month==m)&available).sum())
                                              for m in range(1,13)},
        "original_observer_count":int(merged.user_id.nunique()),
        "source_encrypted_archive_sha256":SEALED_SHA256,
        "source_original_release_zip_sha256":SOURCE_ZIP_SHA256,
        "original_observation_datetime_not_photo_upload_time":True,
        "privacy_raw_observer_ids_event_dates_and_coords_not_public":True,
        "flowering_phenology_true_anthesis_observed":False,
        "S1_camera_angle_not_in_date_source_reconstruction":True,
    }
    return report,private_output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native",type=Path,required=True)
    p.add_argument("--sealed",type=Path,required=True)
    p.add_argument("--out-private",type=Path,required=True)
    p.add_argument("--out-public",type=Path,required=True)
    x=p.parse_args()
    summary,private=original_dates(x.native,x.sealed,x.out_private)
    x.out_public.parent.mkdir(parents=True,exist_ok=True)
    x.out_public.write_text(json.dumps(summary,ensure_ascii=False,
                                      indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("ORIGINAL_S1_PHENOLOGY_DATE_AVAILABILITY="+json.dumps(summary,ensure_ascii=False))
    if private is None:
        raise SystemExit("Authentic event dates absent: no flowering proxy inferred")


if __name__=="__main__":
    main()

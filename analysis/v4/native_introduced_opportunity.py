"""Audit 665115 original Cirsium observations by WCVP-native/introduced status.

SOURCE: original 2026-07-03 encrypted photo metadata and the exact June-2026
WCVP-350/TDWG3 polygons; original archived image-derived observation measures.
This executes on a private GitHub Actions runner. Never publishes original IDs,
point coordinates, users, photos, access tokens, or decrypted source archives.

'outside native' is NOT a synonym for introduced. Original native-first cohort
is reproduced exactly as an invariant; full 665115 observations are classified
before any trait/detector/photographic screening. Introduced source observations
are identified only when establishmentMeans explicitly says INTRODUCED and
the same wild/georeferenced source-inclusion rules are met.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_range_status import classify_roster

ORIGINAL_SHA = "2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae"
ORIGINAL_SEALED_SHA = "dfa98fa162a8029cb6477adf4c26ffc5819cc916a0513d32da14d78d845ab4b0"
ORIGINAL_AAD = b"azami-original-metadata-v1"
AUTHORITY_SHA = "a5bf23efffbf8e8f9459d5c3391a0136adcedd8f9dbff0117caa969ab545c2bd"
AUTHORITY_ARTIFACT_ZIP_SHA = "ec0e3b67c2fb81d4d1dd8a9215585e4523abcbf8085959854ec80dc826a56ca1"
MEASUREMENT_SEALED_SHA = "111c8fb668265a623c948b635dea68b5182f0a2be04d91ec03bda056b608eb60"
MEASUREMENT_AAD = b"azami-observation-measurements-v1"
MEASUREMENT_SHA = "b49a093c849928c456d6747ec51772b5a6b78debb2ae5eb5dbe2a91b08f66dd7"
SOURCE_N = 665115
SOURCE_PHOTOS = 1122854
NATIVE_N = 313350
NATIVE_PHOTOS = 534186
NATIVE_ACCEPTED = 306
TARGET = ("Cirsium arvense","Cirsium vulgare")
CONTINENTS = {
    1:"Europe", 2:"Africa", 3:"Asia temperate",
    4:"Asia tropical", 5:"Australasia", 6:"Pacific",
    7:"North America", 8:"South America", 9:"Antarctic"
}
MIN_ANGLE = 30
MIN_CELLS = 5
MIN_OBSERVERS = 10


def digest(raw:bytes)->str:
    return hashlib.sha256(raw).hexdigest()


def decrypt_checked(raw:bytes,key_hex:str,sha_sealed:str,sha_clear:str,aad:bytes)->bytes:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
    if digest(raw)!=sha_sealed:
        raise ValueError("Encrypted original source identity changed")
    if len(key_hex)!=64:
        raise ValueError("Source decryption key unavailable to authorised runner")
    try:
        output=AESGCM(bytes.fromhex(key_hex)).decrypt(raw[:12],raw[12:],aad)
    except (ValueError,KeyError):
        raise ValueError("Original source authentication failed") from None
    if digest(output)!=sha_clear:
        raise ValueError("Decrypted original source SHA mismatch")
    return output


def authority_source(blob:bytes)->dict:
    if digest(blob)!=AUTHORITY_ARTIFACT_ZIP_SHA:
        raise ValueError("Changed June-2026 WCVP authority artifact ZIP")
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        raw=z.read("native_sources_350.json")
    if digest(raw)!=AUTHORITY_SHA:
        raise ValueError("Changed authority 350 source or TDWG geometry")
    source=json.loads(raw)
    if len(source["resolutions"])!=350 or "geometry" not in source:
        raise ValueError("WCVP/TDWG authority unexpectedly incomplete")
    return source


def source_photos(blob:bytes)->pd.DataFrame:
    columns=[
        "obs_id","photo_id","taxon_name","taxon_rank","observed_on",
        "captive","obscured","latitude","longitude",
        "coordinate_usable_for_environment","quality_grade",
        "positional_accuracy","geoprivacy","user_id"
    ]
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if "photo_metadata_merged.csv" not in z.namelist():
            raise ValueError("Expected original iNaturalist observation photos absent")
        with z.open("photo_metadata_merged.csv") as fh:
            data=pd.read_csv(fh,usecols=columns,
                dtype={"obs_id":str,"photo_id":str,"user_id":str},low_memory=False)
    if len(data)!=SOURCE_PHOTOS or data.obs_id.nunique()!=SOURCE_N:
        raise ValueError("Original full Cirsium photo/observation census changed")
    if data.photo_id.isna().any() or not data.photo_id.is_unique:
        raise ValueError("Original photo identity inconsistent")
    if (data.groupby("obs_id").user_id.nunique(dropna=False)>1).any():
        raise ValueError("Original observation attributed to multiple observers")
    return data


def taxon_authority_with_placeholders(master:pd.DataFrame,auth:dict)->dict:
    source=dict(auth)
    res=list(source["resolutions"])
    queries={r["input_name"] for r in res}
    not_queried=master.loc[~master.taxon_name.isin(queries),
                          ["taxon_name","taxon_rank"]].drop_duplicates()
    if (not_queried.taxon_rank=="species").any():
        raise ValueError("Unqueried species-rank name in original population")
    for name in sorted(not_queried.taxon_name.unique()):
        res.append({
            "input_name":str(name),"resolution_status":"unresolved_taxon",
            "accepted_key":None,"accepted_name":"",
            "n_exact_records":0,"n_accepted_targets":0
        })
    if len(not_queried)!=117:
        raise ValueError("Original non-species placeholder source number changed")
    source["resolutions"]=res
    return source


def source_continent(source:dict)->dict:
    features=source["geometry"]["features"]
    values={}
    for row in features:
        p=row["properties"]
        level=str(p["LEVEL3_COD"]).upper()
        group=int(p["LEVEL1_COD"])
        if level in values and values[level]!=group:
            raise ValueError("TDWG level 3 assigned two continents")
        if group not in CONTINENTS:
            raise ValueError("Unrecognised TDWG level-1 geographic region")
        values[level]=group
    return values


def source_roster(photos:pd.DataFrame,auth:dict)->pd.DataFrame:
    duplicated=photos.groupby("obs_id",sort=False).size()
    master=photos.drop_duplicates("obs_id").drop(columns=["photo_id"]).copy()
    master["source_photo_count"]=master.obs_id.map(duplicated).astype(int)
    with_placeholders=taxon_authority_with_placeholders(master,auth)
    classified=classify_roster(master,with_placeholders)
    if len(classified)!=SOURCE_N or classified.obs_id.nunique()!=SOURCE_N:
        raise ValueError("Source classification changed original observations")
    native=classified.loc[classified.native_primary]
    if (len(native),int(native.source_photo_count.sum()),
        int(native.accepted_key.nunique()))!=(
        NATIVE_N,NATIVE_PHOTOS,NATIVE_ACCEPTED):
        raise ValueError("Pre-detector native denominator/source identity changed")
    level1=source_continent(auth)
    classified["tdwg1"]=classified.tdwg_level3_code.map(level1)
    if ((classified.tdwg_level3_code!="")&classified.tdwg1.isna()).any():
        raise ValueError("Some original TDWG3 geographical cells have no parent")
    return classified


def append_observation_measures(frame:pd.DataFrame,clear:bytes|None)->pd.DataFrame:
    data=frame.copy()
    data["has_archived_measurement_row"]=False
    data["orientation_eligible"]=False
    data["colour_eligible"]=False
    if clear is None:
        return data
    fields=["obs_id","orientation_angle_degrees_median",
            "n_usable_heads_orientation","corolla_lab_chroma_median",
            "n_usable_heads_colour"]
    measured=pd.read_csv(io.BytesIO(clear),usecols=fields,dtype={"obs_id":str},
                        low_memory=False)
    if not measured.obs_id.is_unique or measured.obs_id.isna().any():
        raise ValueError("Archived measurements duplicate observation identity")
    measured["orientation_eligible"]=(
        pd.to_numeric(measured["n_usable_heads_orientation"],errors="coerce").gt(0)
        &pd.to_numeric(measured["orientation_angle_degrees_median"],
                       errors="coerce").between(0,180)
    )
    measured["colour_eligible"]=(
        pd.to_numeric(measured["n_usable_heads_colour"],errors="coerce").gt(0)
        &np.isfinite(pd.to_numeric(measured["corolla_lab_chroma_median"],
                                   errors="coerce"))
    )
    original_ids=data.obs_id
    keep=measured[["obs_id","orientation_eligible","colour_eligible"]]
    data=data.drop(columns=["orientation_eligible","colour_eligible",
                            "has_archived_measurement_row"])
    data=data.merge(keep,on="obs_id",how="left",indicator=True,validate="one_to_one")
    data["has_archived_measurement_row"]=data["_merge"].eq("both")
    data["orientation_eligible"]=data.orientation_eligible.fillna(False).astype(bool)
    data["colour_eligible"]=data.colour_eligible.fillna(False).astype(bool)
    data=data.drop(columns=["_merge"])
    if len(data)!=SOURCE_N or not data.obs_id.equals(original_ids):
        raise ValueError("Appended trait records changed source order/grain")
    return data


def angular_cells(lat:pd.Series,lon:pd.Series)->pd.Series:
    a=pd.to_numeric(lat,errors="coerce").to_numpy(float)
    b=pd.to_numeric(lon,errors="coerce").to_numpy(float)
    if not (np.isfinite(a).all() and np.isfinite(b).all()
            and (np.abs(a)<=90).all() and (np.abs(b)<=180).all()):
        raise ValueError("Geographic cell calculation expects valid coordinates")
    y=np.clip(np.floor((np.sin(np.deg2rad(a))+1)/.05).astype(int),0,39)
    x=np.floor(((b+180)%360)/5).astype(int)
    return pd.Series(y*72+x,index=lat.index)


def range_labels(frame:pd.DataFrame)->pd.Series:
    labels=pd.Series("unresolved_or_missing_TDWG",index=frame.index,dtype=str)
    label=frame.native_range_status.astype(str)
    labels.loc[label.eq("native")]="WCVP_native_TDWG"
    labels.loc[label.eq("introduced")]="WCVP_introduced_TDWG"
    labels.loc[label.isin(["ambiguous","uncertain","unlisted"])]=(
        "WCVP_uncertain_or_unlisted")
    labels.loc[label.eq("unresolved_taxon")]="taxon_unresolved"
    labels.loc[label.eq("unmapped_tdwg")]="unmapped_region"
    labels.loc[label.eq("invalid_coordinates")]="invalid_coordinates"
    # All non-wild/obscured cases remain in denominator, never accepted
    # into wild introduced candidates irrespective of WCVP distribution.
    reason=frame.eligibility_reason.astype(str)
    labels.loc[reason.eq("captive")]="captive"
    labels.loc[reason.eq("captive_unknown")]="captive_unknown"
    labels.loc[reason.eq("obscured")]="obscured"
    labels.loc[reason.eq("coordinates_unusable")]="coordinates_unusable"
    return labels


def _record_group(g:pd.DataFrame)->dict:
    geo=g.dropna(subset=["latitude","longitude"])
    valid=geo.latitude.between(-90,90)&geo.longitude.between(-180,180)
    geo=geo.loc[valid]
    valid_dates=pd.to_datetime(g.observed_on,errors="coerce",format="mixed")
    cells=angular_cells(geo.latitude,geo.longitude) if len(geo) else pd.Series(dtype=int)
    q=g.quality_grade.fillna("unknown").astype(str).value_counts()
    return {
        "source_observations":len(g),
        "source_photo_links":int(g.source_photo_count.sum()),
        "distinct_original_observers":int(g.user_id.nunique()),
        "distinct_equal_sin_5degree_cells":int(cells.nunique()),
        "distinct_TDWG3_regions":int(g.tdwg_level3_code.nunique()),
        "years_observed":int(valid_dates.dt.year.nunique()),
        "dated_source_observations":int(valid_dates.notna().sum()),
        "research_grade":int(q.get("research",0)),
        "needs_id":int(q.get("needs_id",0)),
        "casual":int(q.get("casual",0)),
        "archived_measured_obs":int(g.has_archived_measurement_row.sum()),
        "orientation_eligible_obs":int(g.orientation_eligible.sum()),
        "colour_eligible_obs":int(g.colour_eligible.sum()),
    }


def aggregate(frame:pd.DataFrame)->dict:
    if len(frame)!=SOURCE_N:
        raise ValueError("Unexpected full source count")
    status=range_labels(frame)
    counts={str(k):int(v) for k,v in status.value_counts(dropna=False).items()}
    wild_native=frame.native_primary.eq(True)
    wild_intro=frame.eligibility_reason.eq("range_introduced")
    if int(wild_native.sum())!=NATIVE_N:
        raise ValueError("Native status no longer reproduces original analysis")
    if (wild_intro&wild_native).any():
        raise ValueError("Native-introduced target leakage")
    if wild_intro.sum()<=0:
        raise ValueError("WCVP has no introduced wild observations; classifier possibly broken")
    all_rank=frame.taxon_rank.fillna("unspecified").value_counts()
    overall={
        "original_observations":SOURCE_N,
        "original_photo_links":SOURCE_PHOTOS,
        "native_primary_wild_observations":int(wild_native.sum()),
        "introduced_primary_wild_observations":int(wild_intro.sum()),
        "excluded_or_undetermined_observations":int((~(wild_native|wild_intro)).sum()),
        "source_rank_counts":{str(k):int(v) for k,v in all_rank.items()},
        "source_status_counts_with_original_wild_validation":counts,
        "status_partition_conserved":int(wild_native.sum()+wild_intro.sum()+
            (~(wild_native|wild_intro)).sum())==SOURCE_N,
        "native_primary_accepted_taxa":NATIVE_ACCEPTED,
        "source_measurement_coverage":{
            "total_preexisting_measured_rows":int(frame.has_archived_measurement_row.sum()),
            "native_measured_rows":int(frame.loc[wild_native,"has_archived_measurement_row"].sum()),
            "introduced_measured_rows":int(frame.loc[wild_intro,"has_archived_measurement_row"].sum()),
            "introduced_orientation_eligible_obs":int(frame.loc[wild_intro,"orientation_eligible"].sum()),
            "introduced_colour_eligible_obs":int(frame.loc[wild_intro,"colour_eligible"].sum()),
        }
    }
    if overall["source_measurement_coverage"]["native_measured_rows"]!=217104:
        raise ValueError("Published native measurement denominator changed")
    focal=[]
    for name in TARGET:
        species=frame.loc[frame.taxon_name.eq(name)].copy()
        if species.empty:
            raise ValueError("Expected source focal species not recorded")
        for typ,mask in (
            ("native",wild_native.loc[species.index]),
            ("introduced",wild_intro.loc[species.index]),
            ("other",~(wild_native.loc[species.index]|wild_intro.loc[species.index]))
        ):
            block=species.loc[mask].copy()
            focal.append({"source_name":name,"range_state":typ,
                          "TDWG_lvl1_group":"ALL",**_record_group(block)})
            if typ not in {"native","introduced"}:
                continue
            for level,group in block.groupby("tdwg1",dropna=False,sort=True):
                if pd.isna(level):
                    raise ValueError("Native/introduced source has missing TDWG1 parent")
                focal.append({
                    "source_name":name,"range_state":typ,
                    "TDWG_lvl1_group":CONTINENTS[int(level)],
                    "TDWG_lvl1_code":int(level),
                    **_record_group(group)
                })
    suff=[]
    for name in TARGET:
        d={r["range_state"]:r for r in focal
            if r["source_name"]==name and r["TDWG_lvl1_group"]=="ALL"}
        can=all(d[v]["orientation_eligible_obs"]>=MIN_ANGLE and
                d[v]["distinct_equal_sin_5degree_cells"]>=MIN_CELLS and
                d[v]["distinct_original_observers"]>=MIN_OBSERVERS for v in ("native","introduced"))
        supp={
            "source_name":name,
            "both_range_orientation_supported":bool(can),
            "min_orientation_eligible_per_range":MIN_ANGLE,
            "min_independent_5_degree_cells_per_range":MIN_CELLS,
            "min_original_observers_per_range":MIN_OBSERVERS,
            "different_focal_source_names_not_automatically_pooled":True
        }
        suff.append(supp)
    return {
        "schema_version":"v4_source665115_wcvp_introduced_rain_transfer_opportunity_v1",
        "status":"SOURCE_CLASSIFIED_NOT_A_NATIVE_TO_INTRODUCED_ENVIRONMENTAL_FIT",
        "original_metadata_sha256":ORIGINAL_SHA,
        "source_authority_json_sha256":AUTHORITY_SHA,
        "historical_full_observation_measurement_sha256":MEASUREMENT_SHA,
        "snapshot":"2026-07-03 originals; 2026-06-10 WCVP authority",
        "overall":overall,
        "focal_source_species_by_TDWG_lvl1":focal,
        "focal_native_intro_orientation_support":suff,
        "source_wild_status_rule":"Exact WCVP species→accepted taxon and unique TDWG3; native blank-or-NATIVE; introduced explicitly INTRODUCED; no captive/obscured/unusable-coordinate accepted",
        "privacy":"Only species-level and level1 continental aggregates, no user/location/record IDs",
        "next_stage":"Run native-vs-introduced within-species slopes ONLY if both ranges have source-observer and original measurement support and overlapping BIO12 climate records",
        "interpretation_limit":"TDWG establishment classifies native geographic regions, not proof that an individual has invaded or that a trait is genetically inherited."
    }


def execute(original:Path,authority:Path,measurements:Path,output:Path,
            metadata_key_hex:str,measurements_key_hex:str)->dict:
    native_clear=decrypt_checked(
        original.read_bytes(),metadata_key_hex,ORIGINAL_SEALED_SHA,ORIGINAL_SHA,ORIGINAL_AAD)
    master=source_photos(native_clear)
    del native_clear
    authority_data=authority_source(authority.read_bytes())
    frame=source_roster(master,authority_data)
    del master
    measured=decrypt_checked(
        measurements.read_bytes(),measurements_key_hex,
        MEASUREMENT_SEALED_SHA,MEASUREMENT_SHA,MEASUREMENT_AAD)
    frame=append_observation_measures(frame,measured)
    del measured
    result=aggregate(frame)
    output.mkdir(parents=True,exist_ok=True)
    (output/"source665115_introduced_opportunity.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps({
        "status":result["status"],
        "source":SOURCE_N,
        "native":result["overall"]["native_primary_wild_observations"],
        "introduced":result["overall"]["introduced_primary_wild_observations"],
        "supports":result["focal_native_intro_orientation_support"],
    },ensure_ascii=False))
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("original","authority","measurements","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    execute(a.original,a.authority,a.measurements,a.output,
        os.environ.get("AZAMI_METADATA_KEY",""),
        os.environ.get("AZAMI_MEASUREMENTS_KEY",""))


if __name__=="__main__":
    main()

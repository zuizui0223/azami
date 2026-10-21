"""Phenology eligibility *before* flowering-season precipitation extraction.

Source: exactly the original SHA-pinned 313350 native observation table and
original users. No angle values, environmental outcomes, ancestral states or
photographic colour *values* enter season definition. Flowering is only a
PROXY: measurable visible floret colour, not independently scored anthesis.
Original observation dates cannot estimate flowering duration without bias.

Private outputs retain authority keys/region/month weight rows for subsequent
12 CHELSA monthly precipitation rasters; only aggregated diagnostics become
public workflow artifacts. No source IDs, original photographers or precise
coordinates are published.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

from analysis.v4.native_gates_ae import authority_mapping
from analysis.v4.run_native_abcde import EXPECTED_SHA
from analysis.v4.native_observers import EXPECTED_NATIVE_N

SEASON_WINDOW=3
REGION_LAT_DEGREES=20
REGION_LON_DEGREES=90
MIN_EVENTS=30
MIN_UNIQUE_OBSERVERS=10
MIN_INDEPENDENT_5DEG_CELLS=3
MIN_YEARS=2
MIN_OBSERVED_MONTHS=2
MIN_ORIENTATION_ROWS=10
FOCAL=("Cirsium_arvense","Cirsium_vulgare")


def date_columns(frame:pd.DataFrame)->pd.DataFrame:
    if "observed_on" not in frame:
        raise ValueError("Original observation dates are absent")
    dates=pd.to_datetime(frame["observed_on"],errors="coerce",format="mixed")
    result=frame.copy()
    result["_month"]=dates.dt.month
    result["_year"]=dates.dt.year
    valid=(result["_year"].between(2000,2026)
           &result["_month"].between(1,12))
    result.loc[~valid,["_year","_month"]]=np.nan
    lat= pd.to_numeric(result.latitude,errors="coerce")
    lon= pd.to_numeric(result.longitude,errors="coerce")
    if not lat.between(-90,90).all() or not lon.between(-180,180).all():
        raise ValueError("Native coordinates outside Earth's degree ranges")
    result["_latband"]=np.clip(
        np.floor((lat+90)/REGION_LAT_DEGREES).astype(int),0,8)
    result["_lonband"]=np.clip(
        np.floor((lon+180)/REGION_LON_DEGREES).astype(int),0,3)
    # photo-pose and efflorescence evidence are not used to choose region.
    result["_site_cell"]=(
        np.clip(np.floor((np.sin(np.deg2rad(lat))+1)/.05).astype(int),0,39)
        *72 + np.floor(((lon+180)%360)/5).astype(int)
    )
    return result


def flowering_windows(frame:pd.DataFrame,mapping:list[dict],
                      *,min_events:int=MIN_EVENTS,
                      min_observers:int=MIN_UNIQUE_OBSERVERS,
                      min_cells:int=MIN_INDEPENDENT_5DEG_CELLS,
                      min_years:int=MIN_YEARS,
                      min_months:int=MIN_OBSERVED_MONTHS,
                      min_orientation:int=MIN_ORIENTATION_ROWS):
    """Return a private taxon×20°lat×90°lon window and public aggregate audit.

    The sole phenology indicator is colour__eligible (image-visible floral
    proxy), *never the original pixel colour or head orientation magnitude*.
    One observer×cell×year×month source event is counted at most once.
    Three-month windows are circular and fixed; winter across December works.
    """
    if "user_id" not in frame or frame.user_id.isna().any():
        raise ValueError("Original observer identities required; do not impute")
    if frame.obs_id.isna().any() or not frame.obs_id.is_unique:
        raise ValueError("Source observations not unique")
    required={
        "accepted_key","obs_id","taxon_name","corolla_lab_chroma__eligible",
        "orientation_image_vertical_angle__eligible","corolla_lab_chroma",
        "orientation_image_vertical_angle","user_id","latitude","longitude"
    }
    if not required<=set(frame):
        raise ValueError("Missing required source/evidence fields")
    if frame.accepted_key.nunique()!=len(mapping):
        raise ValueError("Original accepted-key census changed")
    mixed={str(row["accepted_key"]) for row in mapping
           if len(row["source_names"])>1}
    if len(mixed)!=7:
        raise ValueError("Frozen seven merged taxon-concept keys changed")
    frame=date_columns(frame)
    valid=frame._month.notna()
    col=frame["corolla_lab_chroma__eligible"].eq(True)&np.isfinite(
        pd.to_numeric(frame.corolla_lab_chroma,errors="coerce"))
    orient=frame["orientation_image_vertical_angle__eligible"].eq(True)&np.isfinite(
        pd.to_numeric(frame.orientation_image_vertical_angle,errors="coerce"))
    # Exclude taxon-concept collisions for climatic inference, retain counts
    # for the pre-trait full native denominator and region-season missingness.
    clean=frame.loc[valid & ~frame.accepted_key.astype(str).isin(mixed)].copy()
    clean["_floret_proxy"]=col.loc[clean.index].to_numpy()
    clean["_orientation_eligible"]=orient.loc[clean.index].to_numpy()
    idcols=["accepted_key","_latband","_lonband"]
    records=[]
    tested=0
    for key,group in clean.groupby(idcols,sort=True):
        test=group.loc[group["_floret_proxy"]]
        # A single observer repeatedly photographing a patch within one
        # month/year/site is *one* seasonal evidence event.
        independent=(test.sort_values("obs_id").drop_duplicates(
            ["user_id","_site_cell","_year","_month"]))
        tested+=1
        events=len(independent)
        observers=independent.user_id.nunique()
        cells=independent._site_cell.nunique()
        years=independent._year.nunique()
        months=independent._month.nunique()
        orientations=int(group["_orientation_eligible"].sum())
        qualified=(events>=min_events and observers>=min_observers
                   and cells>=min_cells and years>=min_years
                   and months>=min_months and orientations>=min_orientation)
        if not qualified:
            continue
        counts=np.bincount(independent._month.to_numpy(int),minlength=13)[1:]
        windows=[sum(counts[(i+j)%12] for j in range(SEASON_WINDOW)) for i in range(12)]
        # If equally supported, deterministic early-calendar tie break;
        # orientation and local rain are never inspected.
        peak=int(np.argmax(windows))
        onset=peak+1
        peak_months=[int((peak+j)%12+1) for j in range(SEASON_WINDOW)]
        off_months=[j for j in range(1,13) if j not in peak_months]
        records.append({
            "accepted_key":str(key[0]),
            "latband_20":int(key[1]),"lonband_90":int(key[2]),
            "month_start":onset,
            "season_month_1":peak_months[0],
            "season_month_2":peak_months[1],
            "season_month_3":peak_months[2],
            "source_floret_proxy_observer_cell_month_events":events,
            "source_unique_observers":int(observers),
            "source_distinct_cells":int(cells),
            "source_distinct_years":int(years),
            "source_distinct_observed_months":int(months),
            "source_orientation_eligible":orientations,
            "season_event_fraction":float(windows[peak]/events),
            "phenology_basis":"visible_floret_image_proxy_not_observed_flowering",
        })
    private=pd.DataFrame.from_records(records)
    if not private.empty and private.duplicated(["accepted_key","latband_20","lonband_90"]).any():
        raise ValueError("Duplicate species-region phenological windows")
    focal=[]
    for name in FOCAL:
        subset=frame.loc[frame.taxon_name.astype(str).str.replace(" ","_",regex=False).eq(name)]
        focal.append({
            "source_name":name.replace("_"," "),
            "native_source_observations":len(subset),
            "dated_observations":int(subset["_month"].notna().sum()),
            "floret_proxy_eligible_dated":int((
                col.loc[subset.index] & subset["_month"].notna()).sum()),
            "orientation_eligible_dated":int((
                orient.loc[subset.index]&subset["_month"].notna()).sum()),
            "native_regions_20_by_90":int(
                subset[["_latband","_lonband"]].drop_duplicates().shape[0]),
        })
    output={
        "status":"phenology_windows_estimated_from_unvalidated_floret_image_proxy",
        "source_native_observations":len(frame),
        "native_accepted_taxa":int(frame.accepted_key.nunique()),
        "valid_native_dates":int(valid.sum()),
        "source_floret_proxy_eligible_dated":int((col&valid).sum()),
        "source_orientation_eligible_dated":int((orient&valid).sum()),
        "source_both_orientation_and_colour_eligible_dated":int((col&orient&valid).sum()),
        "source_undated":int((~valid).sum()),
        "source_dates_was_month_after_2021":int((frame._year>2021).sum()),
        "taxon_concept_collision_keys_excluded":len(mixed),
        "taxon_region_combinations_examined":tested,
        "phenology_qualified_taxon_region_combinations":len(private),
        "phenology_qualified_unique_accepted_keys":int(private.accepted_key.nunique()) if len(private) else 0,
        "reference_month_count":SEASON_WINDOW,
        "photo_floret_indicator":"corolla_lab_chroma__eligible AND finite chroma; not independently validated full anthesis",
        "region_definition":"20-degree latitude band x 90-degree longitude quadrant; all months circular",
        "effort_constraint":"one observer×5-degree equal-sin site cell×year×month evidence event",
        "min_support":{"events":min_events,"observers":min_observers,
            "cells":min_cells,"years":min_years,"distinct_months":min_months,
            "orientation_rows_in_same_native_region":min_orientation},
        "season_coverage_fraction_distribution":{
            "q10":float(private.season_event_fraction.quantile(.1)) if len(private) else None,
            "median":float(private.season_event_fraction.median()) if len(private) else None,
            "q90":float(private.season_event_fraction.quantile(.9)) if len(private) else None,
        },
        "potential_two_focal_invader_source_names":focal,
        "per_photo_head_identity_verified":False,
        "seasonal_precipitation_from_monthly_CHELSA_extracted":False,
        "negative_controls_image_background_measured":False,
        "claim_boundary":[
            "This output describes the observed month of eligible floral images, not independently validated species flowering.",
            "Search effort, flowering phenology and colour measurement suitability are not separable in this frame.",
            "Window choices are not based on environmental rainfall or head angle.",
            "Month-specific precipitation normals are missing from original nine-raster CHELSA caches.",
            "Season is region-specific, so transhemispheric months are never naively pooled.",
            "Photographic image vertical is not botanical gravity inclination.",
            "2022-2026 source observations require 1981-2010 precipitation to be interpreted only as long-term climate, never observed visit-day rain.",
        ],
    }
    return output,private


def execute(table:Path,observer:Path,mapping:Path,output:Path,private_output:Path):
    digest=hashlib.sha256(table.read_bytes()).hexdigest()
    if digest!=EXPECTED_SHA:
        raise ValueError("Unknown native observation input; source changed")
    frame=pd.read_csv(table,dtype={"obs_id":str,"accepted_key":str},low_memory=False)
    obs=pd.read_csv(observer,dtype={"obs_id":str,"user_id":str},low_memory=False)
    if len(frame)!=EXPECTED_NATIVE_N or not frame.obs_id.is_unique:
        raise ValueError("Frozen native source census changed")
    if len(obs)!=len(frame) or obs.obs_id.duplicated().any() or set(obs.obs_id)!=set(frame.obs_id):
        raise ValueError("Original observer IDs not one-to-one with native sample")
    frame=frame.merge(obs[["obs_id","user_id"]],how="left",on="obs_id",validate="one_to_one")
    if frame.user_id.isna().any():
        raise ValueError("Source observer missing")
    mapping_rows=authority_mapping(mapping)
    results,season=flowering_windows(frame,mapping_rows)
    results["native_source_sha256"]=digest
    results["registered_taxon_mapping_git_blob_sha1"]="b662f9687f11d4069965531d5214962b07f6df8c"
    output.mkdir(parents=True,exist_ok=True)
    private_output.mkdir(parents=True,exist_ok=True)
    (output/"flowering_phenology_opportunity.json").write_text(
        json.dumps(results,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    season.to_csv(private_output/"flowering_windows_private.csv",index=False)
    print(json.dumps({
        "status":results["status"],
        "native_records":len(frame),
        "dated_floret_eligible":results["source_floret_proxy_eligible_dated"],
        "qualified_taxon_regions":len(season),
        "qualified_taxa":results["phenology_qualified_unique_accepted_keys"],
        "precipitation":"not_calculated_without_real_12_month_CHELSA_rasters"
    }),flush=True)
    return results


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("table","observer","mapping","output","private_output"):
        p.add_argument("--"+key.replace("_","-"),required=True,type=Path)
    a=p.parse_args()
    execute(a.table,a.observer,a.mapping,a.output,a.private_output)

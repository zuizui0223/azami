"""Complete WCVP-available-environment PCA-env for two source-supported Cirsium.

Reanalyse the original 665115-photo-source native/explicit-introduced regions
BEFORE orientation measurement selection, on their full TDWG3 available
environmental backgrounds (not the previous native-occupied 316 tile cache).

For each of the five already-approved CHELSA v2.1 original TIFFs, download
the full *official source* once to the PRIVATE runner, verify the exact
original SHA256 + affine transform + scale/offset, then immediately extract
ONLY values at source wild sightings and 0.5-degree WCVP land background
centres, and DELETE the original TIFF. No imputed or resampled rasters.
Private coordinates and original human/observation IDs NEVER get uploaded.

Two-dimensional Gaussian-smoothed PCA-env is an exploratory approximation;
not bitwise-identical to the published R ecospat implementation and not
an adaptation, head-orientation or causal invasion experiment.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import rasterio
from rasterio.transform import rowcol

from analysis.v4.native_climate_source_archive import choose_url_paths, tile_window
from analysis.v4.native_cached_climate import VARIABLES, digest_file
from analysis.v4.native_introduced_opportunity import (
    ORIGINAL_SEALED_SHA, ORIGINAL_SHA, ORIGINAL_AAD,
    decrypt_checked, source_photos, authority_source, source_roster,
    MEASUREMENT_SEALED_SHA, MEASUREMENT_SHA, MEASUREMENT_AAD
)
from analysis.v4.native_full_rain_head_response import analyse_photo_orientation
from analysis.v4.native_wet_edge import save as save_wet_edge
from analysis.v4.native_pca_env import (
    SPECIES_KEYS, SOURCE_COUNTS, CONTINENTS, background_mesh,
    _background_species_audit, _pca_env_2d, climatic_shift_summary,
    _native_and_introduced_codes, _weights_per_5degree_cell
)

LIMIT_SOURCE_BYTES=4_500_000_000
RANGE_STATES=("native","introduced")
EXPECTED_ORIGINAL_WILD_TWO=297622
CACHE_ARCHIVE_MISSING="Do not substitute the historic 316 native-occupied tile subset"


def prepare_original_wild_sources(sealed:Path,authority_archive:Path,
                                  secret_hex:str,*,keep_unbalanced:bool=False):
    clear=decrypt_checked(sealed.read_bytes(),secret_hex,
        ORIGINAL_SEALED_SHA,ORIGINAL_SHA,ORIGINAL_AAD)
    photos=source_photos(clear)
    del clear
    authority=authority_source(authority_archive.read_bytes())
    classified=source_roster(photos,authority)
    del photos
    relevant=classified.loc[
        classified.taxon_name.isin(SPECIES_KEYS) &
        classified.eligibility_reason.isin(("include","range_introduced"))
    ].copy()
    relevant["range_state"]=np.where(
        relevant.native_primary.eq(True),"native","introduced")
    if len(relevant)!=EXPECTED_ORIGINAL_WILD_TWO:
        raise ValueError("Two-species original native and introduced source count changed")
    for (taxon,status),target in SOURCE_COUNTS.items():
        sub=relevant.loc[relevant.taxon_name.eq(taxon) &
                          relevant.range_state.eq(status)]
        if len(sub)!=target or set(sub.accepted_key.astype(str))!={
            SPECIES_KEYS[taxon]}:
            raise ValueError("Original source class or taxon WCVP key changed")
    from analysis.v4.native_introduced_opportunity import angular_cells
    relevant["_cell"]=angular_cells(relevant.latitude,relevant.longitude)
    all_wild=relevant if keep_unbalanced else None
    relevant=relevant.sort_values("obs_id").drop_duplicates(
        ["taxon_name","range_state","_cell","user_id"],keep="first").copy()
    if relevant.obs_id.duplicated().any() or relevant.user_id.isna().any():
        raise ValueError("Wild source IDs/observers changed during balancing")
    # DO NOT require orientation to be scoreable! This is a niche comparison.
    return (relevant,authority,all_wild) if keep_unbalanced else (relevant,authority)


def official_source_variable_name(name:str)->str:
    if name not in VARIABLES or not name.startswith("chelsa_"):
        raise ValueError("Only five registered frozen CHELSA variables allowed")
    return name.removeprefix("chelsa_")


def acquire_verified_source(variable:str,metadata:dict,work:Path)->tuple[Path,dict]:
    key=official_source_variable_name(variable)
    urls=choose_url_paths(key)
    expected=metadata["source_sha256"]
    if len(expected)!=64:
        raise ValueError("Original source raster checksum unavailable")
    target=work/("private_"+key+"_original.tif")
    last=[]
    for url in urls:
        try:
            target.unlink(missing_ok=True)
            sha=hashlib.sha256()
            bytes_received=0
            with requests.get(url,stream=True,timeout=(25,90),
                    headers={"User-Agent":"azami-source-hash-pinned-full-PCA-env/20261021"}
            ) as request:
                request.raise_for_status()
                with target.open("wb") as dest:
                    for chunk in request.iter_content(chunk_size=1<<22):
                        if not chunk:continue
                        bytes_received+=len(chunk)
                        if bytes_received>LIMIT_SOURCE_BYTES:
                            raise ValueError("CHELSA official TIFF exceeded registered source size cap")
                        sha.update(chunk)
                        dest.write(chunk)
            if sha.hexdigest()!=expected:
                raise ValueError("CHELSA source SHA mismatch: reject all new introduced pixels")
            with rasterio.open(target) as ds:
                if ds.crs.to_epsg()!=4326 or ds.count!=1:
                    raise ValueError("CHELSA source TIFF projection/band changed")
                if (ds.width,ds.height)!=(metadata["width"],metadata["height"]):
                    raise ValueError("CHELSA original grid dimensions changed")
                actual=list(ds.transform)[:6]
                if not np.allclose(actual,metadata["transform"],rtol=0,atol=1e-10):
                    raise ValueError("CHELSA original affine transform changed")
                if ds.dtypes[0]!=metadata["dtype"] or ds.nodata!=metadata["nodata"]:
                    raise ValueError("CHELSA source dtype or NODATA changed")
                if not np.isclose(ds.scales[0],metadata["scale"],rtol=0,atol=1e-10):
                    raise ValueError("CHELSA variable scale changed")
                if not np.isclose(ds.offsets[0],metadata["offset"],rtol=0,atol=1e-10):
                    raise ValueError("CHELSA physical offset changed")
            return target,{
                "source_url":url,"source_SHA256":expected,
                "downloaded_TIFF_bytes":bytes_received,
                "original_pixels_only":True,
                "has_been_SHA_verified_against_existing_native_source":True,
            }
        except (requests.RequestException,OSError,ValueError) as exc:
            last.append({"url":url,"error":str(exc)[:240]})
            target.unlink(missing_ok=True)
    raise RuntimeError(
        "Full official source CHELSA acquisition unavailable for "+variable+
        ". This is a source coverage failure, NOT absence of a niche difference: "+
        json.dumps(last)
    )


def sample_original_pixels(source:Path,site:pd.DataFrame)->np.ndarray:
    """Read no more than one 5-degree ORIGINAL TIFF window per populated cell.

    Both occurrences and accessible environment reference centres use exact
    original (1km-ish) integer pixel and original TIFF nodata masking.
    Five-degree windows are ONLY an I/O grouping; no window averaging.
    """
    lat=site.latitude.to_numpy(float)
    lon=site.longitude.to_numpy(float)
    if not(np.isfinite(lat).all() and np.isfinite(lon).all()):
        raise ValueError("Private source sites invalid")
    valid=(np.abs(lat)<=90)&(np.abs(lon)<=180)
    if not valid.all():
        raise ValueError("Source sites outside world longitude/latitude")
    result=np.full(len(site),np.nan,dtype=np.float64)
    ilat=np.floor((lat+90)/5).astype(int)
    ilon=np.floor((lon+180)/5).astype(int)
    tags=np.array([f"{a}_{b}" for a,b in zip(ilat,ilon)],dtype=object)
    with rasterio.open(source) as ds:
        row,col=rowcol(ds.transform,lon,lat)
        rr=np.asarray(row,dtype=int);cc=np.asarray(col,dtype=int)
        for tag in sorted(set(tags)):
            indexes=np.flatnonzero(tags==tag)
            k,l=(int(v) for v in tag.split("_"))
            try:
                window=tile_window(ds,k,l)
            except Exception:
                # Valid latitudes south of physical CHELSA cover (if any) are
                # original unavailable values, NOT numerical zeros.
                continue
            if window.width<=0 or window.height<=0:
                continue
            arr=ds.read(1,window=window,masked=False)
            yy=rr[indexes]-int(window.row_off)
            xx=cc[indexes]-int(window.col_off)
            ok=(yy>=0)&(xx>=0)&(yy<arr.shape[0])&(xx<arr.shape[1])
            z=arr[yy[ok],xx[ok]].astype(float)
            if ds.nodata is not None:
                z[z==float(ds.nodata)]=np.nan
            result[indexes[ok]]=z
    return result


def apply_original_five_CHELSA(wild:pd.DataFrame,land:pd.DataFrame,
                              *,source_meta:Path,work:Path,
                              photo_angles:pd.DataFrame|None=None):
    metadata=json.loads(source_meta.read_text(encoding="utf-8"))
    work.mkdir(parents=True,exist_ok=True)
    site=pd.concat([
        wild[["latitude","longitude"]],
        land[["latitude","longitude"]],
        photo_angles[["latitude","longitude"]]
          if photo_angles is not None else wild[["latitude","longitude"]].iloc[0:0]
    ],ignore_index=True)
    nw=len(wild)
    nl=len(land)
    source_audit={}
    for name in VARIABLES:
        path,receipt=acquire_verified_source(name,metadata[name],work)
        try:
            source_raw=sample_original_pixels(path,site)
            source_phys=source_raw*float(metadata[name]["scale"])+float(metadata[name]["offset"])
            wild[name]=source_phys[:nw]
            land[name]=source_phys[nw:nw+nl]
            if photo_angles is not None:
                photo_angles[name]=source_phys[nw+nl:]
                receipt["orientation_eligible_original_pixel_finite"]=int(
                    np.isfinite(photo_angles[name]).sum())
            receipt["wild_pixels_finite"]=int(np.isfinite(source_phys[:nw]).sum())
            receipt["background_pixels_finite"]=int(np.isfinite(source_phys[nw:nw+nl]).sum())
            receipt["background_pixels_nodata"]=int(np.isnan(source_phys[nw:nw+nl]).sum())
            source_audit[name]=receipt
            if receipt["wild_pixels_finite"]<1000 or receipt["background_pixels_finite"]<1000:
                raise ValueError("New original TIFF did not yield enough valid source pixels")
        finally:
            path.unlink(missing_ok=True)
    if list(metadata) and any(name not in source_audit for name in VARIABLES):
        raise ValueError("All five original CHELSA source TIFFs were not verified")
    return (wild,land,source_audit,photo_angles) if photo_angles is not None else (wild,land,source_audit)


def prepare_verified_photo_angles(all_wild:pd.DataFrame,
                                 sealed_measurement:Path,
                                 secret_hex:str)->pd.DataFrame:
    """Join original archive angle per OBSERVATION on the native+intro wild
    population before geographic/photo-eligibility balancing.
    Native niche estimation stays strictly phenotype-unselected.
    """
    clear=decrypt_checked(
        sealed_measurement.read_bytes(),secret_hex,
        MEASUREMENT_SEALED_SHA,MEASUREMENT_SHA,MEASUREMENT_AAD
    )
    cols=[
        "obs_id","orientation_angle_degrees_median",
        "orientation_angle_degrees_n_usable_heads"
    ]
    measured=pd.read_csv(io.BytesIO(clear),usecols=cols,
                         dtype={"obs_id":str},low_memory=False)
    del clear
    if (len(measured)!=406582 or measured.obs_id.isna().any()
        or not measured.obs_id.is_unique):
        raise ValueError("Original 406582-row measurement archive identity changed")
    original=all_wild.merge(
        measured,on="obs_id",how="left",validate="one_to_one"
    )
    angles=pd.to_numeric(original.orientation_angle_degrees_median,errors="coerce")
    count=pd.to_numeric(original.orientation_angle_degrees_n_usable_heads,errors="coerce")
    photo=original.loc[count.gt(0)&angles.between(0,180)].copy()
    expected={
        ("Cirsium arvense","native"):49072,
        ("Cirsium arvense","introduced"):43195,
        ("Cirsium vulgare","native"):33637,
        ("Cirsium vulgare","introduced"):47833
    }
    observed={tuple(key):len(data) for key,data in photo.groupby(
        ["taxon_name","range_state"],sort=True
    )}
    if observed!=expected or len(photo)!=173737:
        raise ValueError("Both species' original native/introduced head-angle output changed")
    return photo


def _collect_background_result(wild:pd.DataFrame,land:pd.DataFrame,
                                authority:dict)->dict:
    complete=np.isfinite(land[list(VARIABLES)].to_numpy(float)).all(axis=1)
    background=land.loc[complete].copy()
    if background.empty:
        raise ValueError("No complete original-grid WCVP climatic background")
    results={}
    for name in SPECIES_KEYS:
        env,audit=_background_species_audit(authority,background,name,land)
        native_en=env.loc[env._range.eq("native")]
        introduced_en=env.loc[env._range.eq("introduced")]
        sw=wild.loc[wild.taxon_name.eq(name)].copy()
        ok=np.isfinite(sw[list(VARIABLES)].to_numpy(float)).all(axis=1)
        sn=sw.loc[ok &sw.range_state.eq("native")].copy()
        si=sw.loc[ok &sw.range_state.eq("introduced")].copy()
        subject={
            "original_observer_cell_balanced":{
                "native":int((sw.range_state=="native").sum()),
                "introduced":int((sw.range_state=="introduced").sum())
            },
            "all_five_CHELSA_original_pixel_complete":{
                "native":len(sn),"introduced":len(si)
            },
            "background":audit,
            "niche_climate_median_shifts":climatic_shift_summary(sn,si),
            "overall_pca_env":None,
            "by_introduced_continent":{},
        }
        def analyse_if_supported(a,b,o,p,label):
            if (len(a)<200 or len(b)<200 or len(o)<100 or len(p)<100):
                return {"status":"STOP_INSUFFICIENT_EXACT_BACKGROUND_OR_SIGHTINGS"}
            available_fraction=label
            if available_fraction<.8:
                return {"status":"STOP_BACKGROUND_COVERAGE_BELOW_PREDECLARED_80PCT",
                        "available_geographic_background_fraction":available_fraction}
            return _pca_env_2d(a,b,o,p)
        area=min(audit["native"]["complete_pixel_background_fraction_of_all_WCVP_available_area"],
                 audit["introduced"]["complete_pixel_background_fraction_of_all_WCVP_available_area"])
        subject["overall_pca_env"]=analyse_if_supported(
            native_en,introduced_en,sn,si,area)
        for continent,groups in si.groupby("tdwg1",sort=True):
            key=CONTINENTS[int(continent)]
            region_bg=introduced_en.loc[introduced_en.TDWG_level1.eq(int(continent))]
            coverage=audit.get("introduced:"+key,{})
            frac=coverage.get("fraction_original_background_tile_coverage_area_weighted",0.)
            subject["by_introduced_continent"][key]={
                "source_wild_climate_complete_original_sightings":len(groups),
                "original_background_points":len(region_bg),
                "available_background_geographic_fraction":frac,
                "climate_median_shift":climatic_shift_summary(sn,groups)
                    if len(groups)>=100 else {"status":"STOP_OCCUPANCY_TOO_SMALL"},
                "pca_env":analyse_if_supported(native_en,region_bg,sn,groups,
                     min(frac,audit["native"]["complete_pixel_background_fraction_of_all_WCVP_available_area"]))
            }
        results[name]=subject
    return {
        "schema_version":"v4_two_Cirsium_full_TDWG_available_background_PCA_env_20261021_v1",
        "status":"ORIGINAL_TIFF_SHA_VERIFIED_ALL_FIVE_COMPLETE_NATIVE_INTRODUCED_AVAILABLE_ENVIRONMENT",
        "focal_source_species":results,
        "total_wild_source_observations_before_observer_cell_balance":EXPECTED_ORIGINAL_WILD_TWO,
        "n_background_TDWG_land_mesh_points":len(land),
        "n_background_all_five_climates_available":len(background),
        "background_mesh_degrees":.5,
        "background_area_weight":"cos(latitude)",
        "PCA_env_method":"original CHELSA 2.1 full-grid ~1km pixels sampled on 0.5-degree TDWG source-availability mesh; 80x80 joint PC grid, smoothed presence density divided by WCVP available-environment density",
        "not_exact_ecospat_or_MaxEnt":True,
        "without_image_orientation_trait_selection":True,
        "no_genetic_adaptation_or_evolution_inferred":True,
        "no_original_photo_observer_coordinates_in_public_report":True,
        "mean_annual_climate_not_floral_rain_event":True,
    }


def execute(original:Path,authority_archive:Path,source_meta:Path,
            output:Path,work:Path,key:str,
            *,measurements:Path|None=None,measurements_key:str="")->dict:
    if measurements is None:
        wild,auth=prepare_original_wild_sources(original,authority_archive,key)
        photo=None
    else:
        wild,auth,all_wild=prepare_original_wild_sources(
            original,authority_archive,key,keep_unbalanced=True)
        photo=prepare_verified_photo_angles(
            all_wild,measurements,measurements_key)
        del all_wild
    land=background_mesh(auth)
    if photo is None:
        wild,land,audit=apply_original_five_CHELSA(
            wild,land,source_meta=source_meta,work=work)
    else:
        wild,land,audit,photo=apply_original_five_CHELSA(
            wild,land,source_meta=source_meta,work=work,photo_angles=photo)
    report=_collect_background_result(wild,land,auth)
    report["source_TIFF_provenance_per_variable"]=audit
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_two_species_full_background_PCA_env.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
        encoding="utf-8")
    if photo is not None:
        new=analyse_photo_orientation(photo,bootstrap=500)
        new["five_original_CHELSA_SHA256_by_variable"]={
            k:v["source_SHA256"] for k,v in audit.items()
        }
        new["frozen_photo_measurement_archive_sha256"]=MEASUREMENT_SHA
        (output/"S2_two_species_full_rain_head_response.json").write_text(
            json.dumps(new,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
            encoding="utf-8")
        # Outcome-blind source niche threshold; only then matched angle distributions.
        # Existing original TIFFs/measurements are already in-memory: no new acquisition.
        wet_report=save_wet_edge(
            wild,photo,output,
            source_climate_sha={k:v["source_SHA256"] for k,v in audit.items()})
        print(json.dumps({
            "wet_edge_source_complete":wet_report["source_wild_complete_climate_n"],
            "wet_edge_angle_complete":wet_report["orientation_original_image_balanced_climate_complete_n"],
            "status":wet_report["status"],
        },ensure_ascii=False),flush=True)
        print(json.dumps({
            "annual_BIO12_primary_head_image_response":{
                name:{
                    continent:{
                        label:case["status"] for label,case in
                        region["hypotheses"].items()
                    } for continent,region in result["regions"].items()
                } for name,result in new["focal"].items()
            }
        },ensure_ascii=False),flush=True)
    print(json.dumps({
        "status":report["status"],
        "native_background_coverage":{
            name:z["background"]["native"]["complete_pixel_background_fraction_of_all_WCVP_available_area"]
            for name,z in report["focal_source_species"].items()
        },
        "introduced_background_coverage":{
            name:z["background"]["introduced"]["complete_pixel_background_fraction_of_all_WCVP_available_area"]
            for name,z in report["focal_source_species"].items()
        },
        "species":{
            name:{
                "D":z["overall_pca_env"].get(
                    "Schoener_D_available_environment_corrected_density"),
                "status":z["overall_pca_env"]["status"],
            } for name,z in report["focal_source_species"].items()
        }
    },ensure_ascii=False),flush=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("original","authority_archive","source_meta","output","work"):
        p.add_argument("--"+name.replace("_","-"),required=True,type=Path)
    p.add_argument("--measurements",type=Path)
    a=p.parse_args()
    execute(a.original,a.authority_archive,a.source_meta,a.output,a.work,
            os.environ.get("AZAMI_METADATA_KEY",""),
            measurements=a.measurements,
            measurements_key=os.environ.get("AZAMI_MEASUREMENTS_KEY",""))


if __name__=="__main__":
    main()

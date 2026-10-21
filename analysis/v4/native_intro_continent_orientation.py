"""Did Cirsium vulgare capitula change IMAGE orientation after introduction?

Direct, previously unreported descriptive comparison of the SAME original
angle-measuring algorithm between source-native observations and separately
classified introduced continents. A within-continental comparison of
bio01-warm or bio15-high-seasonality pixels is explicitly separate and
requires original CHELSA coverage. This cannot establish genuine gravitational
nodding, longitudinal phenotype change, protection of pollen, rapid evolution,
or selection.

Private input from the authenticated 665115-source workflow ONLY; original
record, user, photo IDs and geographic coordinates MUST NEVER be exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_introduced_climate_transfer import (
    load_verified_focal, dated_and_balanced
)
from analysis.v4.native_cached_climate import (
    physical_climate_for_focal, VARIABLES
)

SEED=20261021
BOOTSTRAP=1000
MIN_DESCRIPTIVE_PHOTOS=30
MIN_DESCRIPTIVE_CELLS=5
MIN_DESCRIPTIVE_OBSERVERS=10
MIN_INTERVAL_CELLS=15
CONTINENTS={1:"Europe",2:"Africa",3:"Asia temperate",4:"Asia tropical",
            5:"Australasia",6:"Pacific",7:"North America",8:"South America",
            9:"Antarctic"}
FOCAL={"Cirsium vulgare","Cirsium arvense"}
TARGET_HYPOTHESES={
    "North America":{"chelsa_bio01":"warmer_than_native_q95"},
    "South America":{"chelsa_bio15":"more_precip_seasonal_than_native_q95"},
    "Africa":{"chelsa_bio15":"more_precip_seasonal_than_native_q95"}
}


def cell_weighted_metrics(frame:pd.DataFrame)->dict:
    if frame.empty:
        return {"status":"NO_PHOTO_MEASUREMENTS","original_image_photos":0}
    if (frame.obs_id.isna().any() or frame.obs_id.duplicated().any() or
            frame.user_id.isna().any() or frame._cell.isna().any()):
        raise ValueError("Private original source photograph identity incomplete")
    response=frame.orientation_angle_degrees_median.to_numpy(float)
    if not np.isfinite(response).all() or not ((response>=0)&(response<=180)).all():
        raise ValueError("Invalid original image-vertical angle")
    per_cell=frame.groupby("_cell",sort=True).orientation_angle_degrees_median.mean()
    n=len(frame); k=int(len(per_cell)); users=int(frame.user_id.nunique())
    research=int(frame.quality_grade.eq("research").sum())
    result={
        "source_photo_orientation_measurements":n,
        "original_distinct_observers":users,
        "original_5degree_equal_sin_cells":k,
        "source_research_grade_fraction":research/n,
        "orientation_photo_observation_equal_mean_degrees":float(response.mean()),
        "orientation_photo_observation_median_degrees":float(np.median(response)),
        "orientation_cell_equal_mean_degrees":float(per_cell.mean()),
        "orientation_cell_equal_median_of_cell_means_degrees":float(per_cell.median()),
        "orientation_cell_mean_sd_degrees":float(per_cell.std(ddof=0)),
        "photo_orientation_quantiles_q10_q50_q90":[float(x) for x in np.quantile(response,[.1,.5,.9])],
        "scope":"EXIF-corrected image-vertical head axis; NOT validated gravity-based nodding or pose-independent anatomy",
    }
    if n<MIN_DESCRIPTIVE_PHOTOS or k<MIN_DESCRIPTIVE_CELLS or users<MIN_DESCRIPTIVE_OBSERVERS:
        result["status"]="INSUFFICIENT_SOURCE_OBSERVERS_OR_SPATIAL_CELLS"
    else:
        result["status"]="DESCRIPTIVE_SOURCE_IMAGE_ANGLE_SUPPORTED"
    return result


def _source_cells(frame:pd.DataFrame)->np.ndarray:
    if frame._cell.nunique()<MIN_INTERVAL_CELLS:
        raise ValueError("At least 15 independent source cells needed for uncertainty")
    return frame.groupby("_cell",sort=True).orientation_angle_degrees_median.mean().to_numpy(float)


def _paired_cell_bootstrap(a:np.ndarray,b:np.ndarray,tag:str,
                           bootstrap:int=BOOTSTRAP)->dict:
    if len(a)<MIN_INTERVAL_CELLS or len(b)<MIN_INTERVAL_CELLS:
        return {"status":"STOP_INSUFFICIENT_SOURCE_CELLS_FOR_INTERVAL",
                "native_cells":len(a),"introduced_cells":len(b),
                "interval":None}
    salt=int.from_bytes(hashlib.sha256(tag.encode("utf-8")).digest()[:8],"little")
    rng=np.random.default_rng(salt ^ SEED)
    draw_a=rng.integers(0,len(a),size=(bootstrap,len(a)))
    draw_b=rng.integers(0,len(b),size=(bootstrap,len(b)))
    boot=b[draw_b].mean(axis=1)-a[draw_a].mean(axis=1)
    if not np.isfinite(boot).all():
        raise ValueError("Invalid original-cell-angle bootstrap")
    return {
        "status":"ORIGINAL_SOURCE_CELL_BOOTSTRAP_DESCRIPTIVE",
        "introduced_minus_native_cell_mean_image_angle_deg":float(b.mean()-a.mean()),
        "bootstrap_site_cell_replicates":bootstrap,
        "site_cell_bootstrap_q025_q50_q975":[float(q) for q in np.quantile(boot,[.025,.5,.975])],
        "bootstrap_intervals_do_not_identify_technical_camera_or_genetic_change":True
    }


def _climate_technical_coverage(frame:pd.DataFrame,
                                variables:dict[str,np.ndarray]) -> dict:
    x=frame.copy()
    for v,arr in variables.items():
        if len(arr)!=len(x):
            raise ValueError("Source climate arrays disagree with original observation order")
        x[v]=arr
    if set(VARIABLES)!=set(variables):
        raise ValueError("Must use all five previously fixed source-climate covariates")
    result={}
    for (species,range_state,continent),g in x.groupby(
        ["taxon_name","range_state","tdwg1"],dropna=False,sort=True
    ):
        if pd.isna(continent):
            raise ValueError("Unknown source continent status")
        name=CONTINENTS[int(continent)]
        status={
            "all_original_photo_angle_observations":len(g),
            "five_original_CHELSA_pixels_complete":int(
                np.isfinite(g[list(VARIABLES)].to_numpy(float)).all(axis=1).sum()),
            "available_original_pixels_by_variable":{
                v:int(np.isfinite(g[v].to_numpy(float)).sum()) for v in VARIABLES
            },
            "does_not_impute_missing_original_pixels":True,
        }
        result[species+"|"+range_state+"|"+name]=status
    return result,x


def hypothesis_within_continent(
        native:pd.DataFrame,introduced:pd.DataFrame,
        *,climate_variable:str,scenario:str) -> dict:
    if climate_variable not in {"chelsa_bio01","chelsa_bio15"}:
        raise ValueError("Only literature-registered warm and precipitation-seasonality axes")
    n=native.loc[np.isfinite(native[climate_variable])].copy()
    i=introduced.loc[np.isfinite(introduced[climate_variable])].copy()
    out={
        "scenario":scenario,
        "source_photos_with_native_reference_climate":len(n),
        "introduced_photos_with_original_climate_pixel":len(i),
        "rainfall_related_literature_target":"warm North America or high BIO15 seasonality Africa/South America",
        "angle_reference":"image vertical, not gravity",
    }
    if (len(n)<100 or n._cell.nunique()<15 or
            len(i)<2*MIN_DESCRIPTIVE_PHOTOS or i._cell.nunique()<2*MIN_DESCRIPTIVE_CELLS):
        out.update(status="STOP_NO_ORIGINAL_CLIMATE_COVERAGE_OR_INDEPENDENT_CELLS",
                   within_continent_angle_change=None)
        return out
    native_cell_counts=n.groupby("_cell").obs_id.transform("size").to_numpy(float)
    w=1/native_cell_counts
    order=np.argsort(n[climate_variable].to_numpy(float),kind="mergesort")
    values=n[climate_variable].to_numpy(float)[order]
    weights=w[order]
    cum=(np.cumsum(weights)-weights/2)/weights.sum()
    threshold=float(np.interp(.95,cum,values))
    higher=i.loc[i[climate_variable]>threshold].copy()
    analogue=i.loc[i[climate_variable]<=threshold].copy()
    a=cell_weighted_metrics(higher)
    b=cell_weighted_metrics(analogue)
    out.update({
        "reference_native_cell_equal_q95":threshold,
        "introduced_beyond_native_q95":a,
        "introduced_at_or_below_native_q95":b
    })
    if any(z["status"]!="DESCRIPTIVE_SOURCE_IMAGE_ANGLE_SUPPORTED" for z in (a,b)):
        out.update(status="STOP_SOURCE_COMPARISON_NOT_IDENTIFIABLE",
                   within_continent_angle_change=None)
        return out
    dif=float(
        a["orientation_cell_equal_mean_degrees"]-
        b["orientation_cell_equal_mean_degrees"]
    )
    difference={
        "novel_warm_or_precip_seasonal_minus_local_other_cell_mean_image_angle":dif,
        "image_angle_difference_not_gravitational_head_nodding":True,
    }
    if (a["original_5degree_equal_sin_cells"]>=MIN_INTERVAL_CELLS and
        b["original_5degree_equal_sin_cells"]>=MIN_INTERVAL_CELLS):
        difference["cell_cluster_bootstrap"]=_paired_cell_bootstrap(
            _source_cells(analogue),_source_cells(higher),
            tag="new_niche_"+scenario+"_"+climate_variable
        )
        out["status"]="SUPPORTED_DESCRIPTIVE_IMAGE_ANGLE_CLIMATE_REGION_CONTRAST"
    else:
        difference["cell_cluster_bootstrap"]={"status":"STOP_LESS_THAN_15_CELLS_PER_EXPOSURE"}
        out["status"]="DESCRIPTIVE_ONLY_INSUFFICIENT_INDEPENDENT_CELLS_FOR_INTERVAL"
    out["within_continent_angle_change"]=difference
    return out


def analyse_regions(
        private_source:pd.DataFrame,*,climate:dict[str,np.ndarray]|None=None
) -> dict:
    if not set(private_source.taxon_name)==FOCAL:
        raise ValueError("Two source-verified taxon names required without synonym pooling")
    if set(private_source.range_state)!={"native","introduced"}:
        raise ValueError("Exact WCVP source status needed")
    original=private_source.copy()
    if len(original)!=173737:
        raise ValueError("Original two-species global source observation census changed")
    balanced=dated_and_balanced(original)
    if balanced.empty or balanced.obs_id.duplicated().any():
        raise ValueError("No valid original observer×cell and time records")
    payload={}
    for name,g in balanced.groupby("taxon_name",sort=True):
        native=g.loc[g.range_state.eq("native")]
        ref=cell_weighted_metrics(native)
        a=_source_cells(native)
        regions={}
        for continent,g2 in g.loc[g.range_state.eq("introduced")].groupby("tdwg1"):
            if pd.isna(continent) or int(continent) not in CONTINENTS:
                raise ValueError("Original WCVP continent not classified")
            title=CONTINENTS[int(continent)]
            measures=cell_weighted_metrics(g2)
            if measures["status"]=="DESCRIPTIVE_SOURCE_IMAGE_ANGLE_SUPPORTED":
                measures["vs_native_cell_mean_angle"]=_paired_cell_bootstrap(
                    a,
                    _source_cells(g2) if g2._cell.nunique()>=MIN_INTERVAL_CELLS else
                    g2.groupby("_cell").orientation_angle_degrees_median.mean().to_numpy(float),
                    tag=name+"|"+title
                )
                measures["vs_native_descriptive_difference_deg"]=float(
                    measures["orientation_cell_equal_mean_degrees"]-
                    ref["orientation_cell_equal_mean_degrees"]
                )
            else:
                measures["vs_native_descriptive_difference_deg"]=None
                measures["vs_native_cell_mean_angle"]={
                    "status":"STOP_UNSUPPORTED_SMALL_SOURCE_REGION"}
            regions[title]=measures
        payload[name]={
            "native_all_regions":ref,
            "introduced_by_TDWG_level1":regions,
            "observed_introduced_regions":len(regions),
            "basis":"matched original image-angle archive; one source observation per taxon×original photographer×5° spatial cell",
        }
    output={
        "schema_version":"v4_S2_original_image_angle_by_introduced_continent_20261021_v1",
        "status":"OBSERVED_IMAGE_RELATIVE_HEAD_ANGLE_BY_NATIVE_AND_INTRODUCED_CONTINENT",
        "full_source_original_cirsium_observations":665115,
        "two_species_archived_eligible_image_angle_observations":len(original),
        "two_species_after_original_observer_cell_date_balancing":len(balanced),
        "source_taxa":payload,
        "primary_interpretation":"Different photo-derived angle distributions across introduced continents are NOT direct longitudinal change or botanical nodding, even when cell-cluster intervals do not overlap zero.",
        "original_observer_spatial_cell_balancing_applied":True,
        "no_private_original_identifiers_in_public_result":True,
        "same_head_repeatability_not_validated":True,
        "introduced_species_population_history_not_established":True,
    }
    if climate is not None:
        coverage,with_climate=_climate_technical_coverage(original,climate)
        with_climate=dated_and_balanced(with_climate)
        comparisons={}
        for name,g in with_climate.groupby("taxon_name",sort=True):
            native=g.loc[g.range_state.eq("native")]
            introduced=g.loc[g.range_state.eq("introduced")]
            for name_region,clim in TARGET_HYPOTHESES.items():
                gregion=introduced.loc[introduced.tdwg1.eq(
                    next(k for k,v in CONTINENTS.items() if v==name_region)
                )]
                for v,label in clim.items():
                    comparisons[name+"|"+name_region+"|"+v]=hypothesis_within_continent(
                        native,gregion,climate_variable=v,scenario=label
                    )
        output["source_pixel_climate_coverage_by_region"]=coverage
        output["climate_novelty_x_photo_angle_by_introduced_region"]=comparisons
    else:
        output["source_pixel_climate_coverage_by_region"]={"status":"NOT_COMPUTED"}
        output["climate_novelty_x_photo_angle_by_introduced_region"]={
            "status":"NOT_COMPUTED_WITHOUT_ORIGINAL_CHELSA_PIXELS"}
    return output


def execute(private_focal:Path,source_receipt:Path,
            output:Path,climate_root:Path|None=None,
            metadata:Path|None=None)->dict:
    original,receipt=load_verified_focal(private_focal,source_receipt)
    climate=None
    if climate_root is not None:
        if metadata is None:
            raise ValueError("CHELSA source metadata required for original pixels")
        climate,audit=physical_climate_for_focal(original,climate_root,metadata)
    result=analyse_regions(original,climate=climate)
    result["source_status_sha256"]=hashlib.sha256(source_receipt.read_bytes()).hexdigest()
    if climate_root is not None:
        result["five_original_CHELSA_grid_sha256"]={
            name:audit[name]["source_raster_sha256"] for name in VARIABLES
        }
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_original_head_angle_by_continent.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
        encoding="utf-8"
    )
    print(json.dumps({
        "status":result["status"],
        "vulgare":{
            k:{
                "photographs":v["source_photo_orientation_measurements"],
                "cells":v["original_5degree_equal_sin_cells"],
                "delta_image_degrees":v["vs_native_descriptive_difference_deg"],
            } for k,v in result["source_taxa"]["Cirsium vulgare"]
                ["introduced_by_TDWG_level1"].items()
        }
    },ensure_ascii=False),flush=True)
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("private_focal","source_receipt","output"):
        p.add_argument("--"+key.replace("_","-"),type=Path,required=True)
    p.add_argument("--climate-root",type=Path)
    p.add_argument("--metadata",type=Path)
    a=p.parse_args()
    execute(a.private_focal,a.source_receipt,a.output,a.climate_root,a.metadata)


if __name__=="__main__":
    main()

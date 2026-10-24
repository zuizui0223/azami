"""Does new wet / seasonally rainy climate occupancy covary with head PHOTO angle?

This is a NEXT analysis following completed full available-background PCA-env,
NOT a replacement for the phenotypically independent ecological-niche analysis.

Within each existing original source-qualified C. vulgare / C. arvense
introduced continent, compare exactly the same source-archive image-axis
head angles on a climate-only split:
  A) original native image-eligible climate BIO12 or BIO15 95th percentile,
  B) introduced continent's own climate-only bottom-vs-top tercile.
Native reference is NOT a new head-angle phenotype threshold. All BIO12/BIO15
values are sampled from original CHELSA v2.1 source TIFFs after SHA256
validation, for the first time including Africa/South America raster pixels.

Original photos, user IDs, WCVP coordinates, raw TIFFs and decrypted metadata
are private in the Actions runner. Only source counts and aggregate comparisons
leave the runner. A documented angle difference is not a gravity-referenced
nodding phenotype, temporal evolutionary response, adaptation or causal rain
shelter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_introduced_climate_transfer import (
    dated_and_balanced,load_verified_focal
)
from analysis.v4.native_intro_continent_orientation import (
    FOCAL, CONTINENTS, _paired_cell_bootstrap, _source_cells,
    cell_weighted_metrics, hypothesis_within_continent,
    MIN_INTERVAL_CELLS, MIN_DESCRIPTIVE_CELLS, MIN_DESCRIPTIVE_PHOTOS,
    MIN_DESCRIPTIVE_OBSERVERS
)
from analysis.v4.native_climate_trait_coupling import weighted_quantiles

PRIMARY_AXIS="chelsa_bio12"  # Annual rain: verified S1 between-taxon association
SECONDARY_AXIS="chelsa_bio15"  # Uncertain S1 tendency: 95% interval includes zero
AXES=(PRIMARY_AXIS,SECONDARY_AXIS)
SOURCE_ORIGINAL_ALL_HEADS=173737
TERCILE_BREAKS=(1/3,2/3)
PREREGISTERED_NATIVE_REFERENCE=.95


def _area_weighted_site_quantiles(group:pd.DataFrame,
                                  variable:str)->tuple[float,float]:
    if group.empty:
        raise ValueError("Cannot define introduced climate tercile from no original photos")
    n=group.groupby("_cell").obs_id.transform("size").to_numpy(float)
    w=1/n
    return tuple(weighted_quantiles(group[variable].to_numpy(float),w,TERCILE_BREAKS))


def climatic_tercile_photo_contrast(source:pd.DataFrame,variable:str,
                                   *,taxon:str,region:str)->dict:
    if variable not in AXES:
        raise ValueError("Rainfall/seasonality axes only")
    base={
        "species":taxon,"introduced_continent":region,"climate_variable":variable,
        "all_original_source_image_angle_rows":len(source),
        "source_cells":int(source._cell.nunique()) if len(source) else 0,
        "status":"STOP_INSUFFICIENT_CLIMATE_ORIGINAL_SOURCE_SUPPORT",
        "photo_angle_tercile_contrast":None
    }
    if len(source)<2*MIN_DESCRIPTIVE_PHOTOS or source._cell.nunique()<MIN_DESCRIPTIVE_CELLS:
        return base
    if source[variable].isna().any() or not np.isfinite(source[variable].to_numpy(float)).all():
        raise ValueError("Original CHELSA missing: never substitute neighbour or zero")
    low,high=_area_weighted_site_quantiles(source,variable)
    lesser=source.loc[source[variable]<low]
    greater=source.loc[source[variable]>high]
    a=cell_weighted_metrics(lesser)
    b=cell_weighted_metrics(greater)
    base.update({
        "introduced_climate_tercile_lower_threshold":float(low),
        "introduced_climate_tercile_upper_threshold":float(high),
        "climate_terciles_are_outcome_blind":True,
        "lower_climate_third":a,"upper_climate_third":b,
    })
    if a["status"]!="DESCRIPTIVE_SOURCE_IMAGE_ANGLE_SUPPORTED" or b["status"]!="DESCRIPTIVE_SOURCE_IMAGE_ANGLE_SUPPORTED":
        return base
    contrast={
        "upper_minus_lower_image_angle_degrees":
            float(b["orientation_cell_equal_mean_degrees"]-
                  a["orientation_cell_equal_mean_degrees"]),
        "not_gravity_angle_or_within_genet_change":True,
    }
    if min(a["original_5degree_equal_sin_cells"],
           b["original_5degree_equal_sin_cells"]) >= MIN_INTERVAL_CELLS:
        contrast["original_site_cell_bootstrap"]=_paired_cell_bootstrap(
            _source_cells(lesser),_source_cells(greater),
            tag="introduced_tercile_"+taxon+"_"+region+"_"+variable
        )
        base["status"]="DESCRIPTIVE_INTRODUCED_CLIMATE_TERCILE_IMAGE_ANGLE_CONTRAST"
    else:
        contrast["original_site_cell_bootstrap"]={
            "status":"STOP_BELOW_15_SOURCE_CELLS_PER_CLIMATE_TERCILE",
        }
        base["status"]="DESCRIPTIVE_ONLY_TOO_FEW_INDEPENDENT_CELLS_FOR_UNCERTAINTY"
    base["photo_angle_tercile_contrast"]=contrast
    return base


def run_angle_seasonality(source:pd.DataFrame,original_meta:dict)->dict:
    if len(source)!=SOURCE_ORIGINAL_ALL_HEADS or set(source.taxon_name)!=FOCAL:
        raise ValueError("Original two species 173737 head-image measurement cohort changed")
    if set(source.range_state)!={"native","introduced"}:
        raise ValueError("Original WCVP native/introduced status changed")
    if source.obs_id.isna().any() or not source.obs_id.is_unique or source.user_id.isna().any():
        raise ValueError("Original head observation/photographer validation failed")
    if not np.isfinite(source[list(AXES)].to_numpy(float)).all():
        raise ValueError("Original BIO12 and BIO15 pixels missing after source download")
    eligible=dated_and_balanced(source)
    if eligible.obs_id.duplicated().any() or not set(AXES)<=set(eligible):
        raise ValueError("Original observer×cell and climate support changed")
    results={}
    for species,group in eligible.groupby("taxon_name",sort=True):
        native=group.loc[group.range_state.eq("native")]
        intro=group.loc[group.range_state.eq("introduced")]
        if len(native)<100 or native._cell.nunique()<15:
            raise ValueError("Native reference poorly sampled: fail closed")
        regions={}
        for continent,regional in intro.groupby("tdwg1",sort=True):
            if pd.isna(continent) or int(continent) not in CONTINENTS:
                raise ValueError("Unknown WCVP original introduced region")
            name=CONTINENTS[int(continent)]
            one={
                "source_original_image_orientation_rows":len(regional),
                "source_unique_observer_cell_5deg":int(regional._cell.nunique()),
                "axes":{},
            }
            for variable in AXES:
                # A. Strict novel climate > native 95% (one-variable only).
                extreme=hypothesis_within_continent(
                    native,regional,climate_variable=variable,
                    scenario="native_q95_"+variable+"_"+species+"_"+name
                )
                # B. Replicated gradient within the introduced region independent
                # of native-climate rectangle shared-support requirements.
                tert=climatic_tercile_photo_contrast(
                    regional,variable,taxon=species,region=name)
                one["axes"][variable]={
                    "above_native_q95_vs_other_within_same_introduced_continent":extreme,
                    "upper_vs_lower_introduced_climate_tercile":tert,
                }
            regions[name]=one
        results[species]={
            "native_source_head_image_angle_observations":len(native),
            "native_source_5deg_cells":int(native._cell.nunique()),
            "introduced_continents":regions,
        }
    return {
        "schema_version":"v4_S2_two_Cirsium_continent_actual_rain_and_seasonality_head_image_angle_20261021_v1",
        "status":"EXECUTED_ORIGINAL_TWO_CHELSA_RAIN_VARIABLES_VS_ORIGINAL_IMAGE_HEAD_ANGLES",
        "original_wild_source_observations_all_Cirsium":665115,
        "source_original_eligible_image_head_angles":SOURCE_ORIGINAL_ALL_HEADS,
        "observer_cell_balanced_image_head_angles":len(eligible),
        "range_status_WCVP":"Wild un-obscured unique species accepted native vs explicit introduced original TDWG3",
        "climate_variables":list(AXES),
        "primary_registered_climate_effect":{
            "axis":PRIMARY_AXIS,
            "why":"S1 among-species association +0.4382 image degrees / 100mm with 95pct interval +0.2058 to +0.6706; secondary precipitation seasonality was not an independently resolved effect",
        },
        "secondary_exploratory_climate_effect":{
            "axis":SECONDARY_AXIS,
            "S1_between_deg_per_plus10_units":0.7880825140436949,
            "S1_between_95pct_interval":[-0.2662791617924732,1.8424441898798627],
            "S1_between_interval_includes_zero":True,
            "S1_within_deg_per_plus10_units":-0.15339124857684425,
            "S1_within_95pct_interval":[-0.5267153435666101,0.2199328464129216],
            "S1_within_interval_includes_zero":True,
            "status":"SECONDARY_SUPPORTING_EXPLORATORY_NOT_CONFIRMED"
        },
        "climate_window":"CHELSA v2.1 original 1981–2010 normals, NOT rainfall on observation date or during anthesis",
        "climate_threshold_selection":"Source-native (within species) 95th percentile and within introduced continent climate terciles, derived only from climate and cell weights",
        "distinct_focal_species":results,
        "full_original_CHELSA_source_SHA_and_grid":original_meta,
        "source_image_angle_reference":"EXIF upright image vertical (0=image-up,90=image-horizontal,180=image-down); NOT botanical gravitational angle",
        "private_observation_user_geo_or_image_data_published":False,
        "conditional_science_limits":[
            "Extreme new-rain and seasonal reference is one marginal climate axis, not full PCA-env novelty nor multivariate climate matching.",
            "Upper/lower introduced-climate tercile comparisons are geographical cross sections, not evolution of the same lineage or an experimental rain shield.",
            "Photo orientation depends on camera pitch/pose, head crop, flowering stage and uploaded image selection.",
            "Photographic independence and source-introduction founding-genotype correspondence are unverified.",
            "Cell bootstrap interval crossing zero is not an equivalence test or proof of evolutionary fixity.",
            "The two climate normals do not identify wetting of pollen, nectar, lifetime fitness or adaptive genetic change.",
        ]
    }


def execute(focal:Path,original_status:Path,metadata:Path,
            work:Path,output:Path)->dict:
    # Rasterio is ONLY required for the original-TIFF Actions execution. CI
    # source-selection and estimand tests run without adding dependencies to
    # the SHA-frozen GEB V3 integrity environment.
    from analysis.v4.native_pca_env_fullcoverage import (
        acquire_verified_source, sample_original_pixels
    )
    frame,source=load_verified_focal(focal,original_status)
    grid=json.loads(metadata.read_text(encoding="utf-8"))
    provenance={}
    work.mkdir(parents=True,exist_ok=True)
    for variable in AXES:
        original,receipt=acquire_verified_source(variable,grid[variable],work)
        try:
            raw=sample_original_pixels(original,frame)
            values=raw*float(grid[variable]["scale"])+float(grid[variable]["offset"])
            if not np.isfinite(values).all():
                raise ValueError("Missing BIO12/BIO15 source pixel on a wild, exact-coordinate, eligible head: STOP")
            frame[variable]=values
            provenance[variable]=receipt
        finally:
            original.unlink(missing_ok=True)
    result=run_angle_seasonality(frame,provenance)
    result["source_original_status_SHA256"]=hashlib.sha256(original_status.read_bytes()).hexdigest()
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_introduced_rain_seasonality_vs_original_head_image_angle.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "source_angle_observations":result["source_original_eligible_image_head_angles"],
        "balanced_angle_observations":result["observer_cell_balanced_image_head_angles"],
        "vulgare":{
            c:{
                ax:{
                    "extreme":z["above_native_q95_vs_other_within_same_introduced_continent"]["status"],
                    "tercile":z["upper_vs_lower_introduced_climate_tercile"]["status"]
                } for ax,z in v["axes"].items()
            } for c,v in result["distinct_focal_species"]["Cirsium vulgare"]["introduced_continents"].items()
        }
    },ensure_ascii=False),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("focal","original_status","metadata","work","output"):
        p.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    execute(a.focal,a.original_status,a.metadata,a.work,a.output)


if __name__=="__main__":
    main()

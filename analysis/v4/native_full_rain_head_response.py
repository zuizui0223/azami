"""Direct test: does a wetter newly occupied climate correspond to head PHOTO angle?

One study, two SOURCE-SUPPORTED Cirsium species, native versus each introduced
continent. Primary exposure is annual precipitation BIO12, secondary BIO15
precipitation seasonality. The climate-only native threshold is frozen at
5th/95th percentiles, cell-equal weighted, BEFORE consulting any head angles.

Model endpoints: original image-vertical 0° up / 90° sideways / 180° down;
NOT a gravity-referenced anatomical angle. Exposures are original-pixel
CHELSA 1981-2010 normals, not rain during anthesis.

Source originals/photo/observer/location IDs stay on the private runner.
Public summaries contain taxon-continent aggregates, no source point data.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_introduced_climate_transfer import (
    dated_and_balanced, cell_weighted_fit, standardization
)
from analysis.v4.native_climate_trait_coupling import (
    weighted_quantiles, original_cell_weights
)
from analysis.v4.native_intro_continent_orientation import CONTINENTS

SPECIES=("Cirsium arvense","Cirsium vulgare")
PRIMARY="chelsa_bio12"
SECONDARY="chelsa_bio15"
MIN_CATEGORY_RECORDS=30
MIN_CATEGORY_OBSERVERS=10
MIN_DESCRIPTIVE_CELLS=5
MIN_INTERVAL_CELLS=15
BOOTSTRAPS=500
SEED=20261021


def climate_reference(native:pd.DataFrame)->dict:
    if len(native)<100 or native["_cell"].nunique()<20:
        raise ValueError("No spatially replicated original native floral-angle source")
    if not set(VARIABLES)<=set(native):
        raise ValueError("Frozen five CHELSA predictors missing")
    weight=original_cell_weights(native)
    result={}
    for name in (PRIMARY,SECONDARY):
        v=native[name].to_numpy(float)
        if not np.isfinite(v).all():
            raise ValueError("Original native climate includes nonfinite values")
        lower,upper=weighted_quantiles(v,weight)
        if not upper>lower:
            raise ValueError("Native climatic envelope degenerate")
        result[name]={"native_cell_equal_q05":lower,"native_cell_equal_q95":upper,
                      "source":"original native original-observer×cell filtered head-angle sample",
                      "angle_ignored_in_threshold":True}
    return result


def _cell_average(group:pd.DataFrame)->dict:
    n=len(group)
    cells=int(group["_cell"].nunique())
    observers=int(group.user_id.nunique())
    supported=n>=MIN_CATEGORY_RECORDS and cells>=MIN_DESCRIPTIVE_CELLS and observers>=MIN_CATEGORY_OBSERVERS
    summary={"original_observations":n,"equal_area_5deg_cells":cells,
             "original_observers":observers,
             "status":"DESCRIPTIVE_SUPPORTED" if supported else "STOP_INSUFFICIENT_RECORDS_CELLS_OR_OBSERVERS"}
    if supported:
        y=group.groupby("_cell",sort=True).orientation_angle_degrees_median.mean()
        summary["image_angle_cell_equal_mean_deg"]=float(y.mean())
        summary["BIO12_cell_equal_mean_mm"]=float(group.groupby("_cell").chelsa_bio12.mean().mean())
        summary["BIO15_cell_equal_mean"]=float(group.groupby("_cell").chelsa_bio15.mean().mean())
    else:
        summary["image_angle_cell_equal_mean_deg"]=None
    return summary


def _bootstrap_shared_cells(
    wet:pd.DataFrame,analogue:pd.DataFrame, *,
    name:str,replicates:int=BOOTSTRAPS
)->dict:
    """Resample distinct source cells, preserving wet/analogue overlap by cell."""
    w=wet.groupby("_cell").orientation_angle_degrees_median.mean()
    a=analogue.groupby("_cell").orientation_angle_degrees_median.mean()
    if min(len(w),len(a))<MIN_INTERVAL_CELLS:
        return {"status":"STOP_LESS_THAN_15_CELLS_EACH_GROUP",
                "wet_cells":len(w),"analogue_cells":len(a),
                "difference":float(w.mean()-a.mean()),
                "interval":None}
    cells=sorted(set(w.index).union(a.index))
    n=len(cells)
    presence=np.zeros((n,2),float)
    total=np.zeros((n,2),float)
    ix={cell:k for k,cell in enumerate(cells)}
    for key,values,col in ((w,w,0),(a,a,1)):
        for cell,value in values.items():
            k=ix[cell]
            presence[k,col]=1.
            total[k,col]=float(value)
    seed=int.from_bytes(hashlib.sha256(name.encode()).digest()[:8],"little")^SEED
    rng=np.random.default_rng(seed)
    valid=0
    differences=[]
    for _ in range(replicates):
        counts=np.bincount(rng.integers(n,size=n),minlength=n)
        n_group=counts@presence
        if np.any(n_group<=0):
            continue
        value=counts@total
        differences.append(float(value[0]/n_group[0]-value[1]/n_group[1]))
        valid+=1
    if valid<.9*replicates:
        return {"status":"STOP_CLUSTER_BOOTSTRAP_INSUFFICIENT_VALID_DRAWS",
                "valid_draws":valid,"difference":float(w.mean()-a.mean()),
                "interval":None}
    q=np.quantile(differences,[.025,.5,.975])
    return {
        "status":"EXPLORATORY_REGIONAL_CELL_BOOTSTRAP",
        "difference_wet_minus_analogue_image_degrees":float(w.mean()-a.mean()),
        "site_cell_bootstrap_q025_q50_q975":[float(z) for z in q],
        "draws_valid":valid,
        "pointwise_interval_not_causal_or_gravity_calibrated":True
    }


def compare_within_continent(region:pd.DataFrame,reference:dict, *,
                             taxon:str,continent:str,
                             bootstrap:int=BOOTSTRAPS)->dict:
    """Compare wetter/higher-BIO15 regional sites with SAME-continent analogues."""
    result={
        "status":"ORIGINAL_SOURCE_IMAGE_ANGLE_AND_CLIMATE_COMPARISON",
        "region":continent,
        "source_original_observer_cell_balanced_records":len(region),
        "source_region_5deg_cells":int(region._cell.nunique()),
        "native_reference":reference,
        "response_axis":"EXIF image-up=0°, horizontal=90°, image-down=180°, not true gravity angle",
        "rain_classes_chosen_without_image_angle":True,
        "hypotheses":{},
    }
    for climate,label in ((PRIMARY,"annual_precipitation"),
                           (SECONDARY,"precipitation_seasonality")):
        lo=reference[climate]["native_cell_equal_q05"]
        hi=reference[climate]["native_cell_equal_q95"]
        wet=region.loc[region[climate]>hi].copy()
        analogue=region.loc[region[climate].between(lo,hi,inclusive="both")].copy()
        dry=region.loc[region[climate]<lo].copy()
        categories={
            "above_native_q95":_cell_average(wet),
            "within_native_q05_q95":_cell_average(analogue),
            "below_native_q05":_cell_average(dry),
        }
        comparison={
            "status":"NOT_ESTIMABLE_CATEGORY_SAMPLE",
            "threshold_native_cell_equal_q95":hi,
            "threshold_native_cell_equal_q05":lo,
            "categories":categories,
            "higher_than_native_q95_minus_native_analogue_same_continent":None,
            "camera_angle_never_used_in_class_selection":True,
        }
        if (categories["above_native_q95"]["status"]=="DESCRIPTIVE_SUPPORTED"
           and categories["within_native_q05_q95"]["status"]=="DESCRIPTIVE_SUPPORTED"):
            comparison["status"]="REGIONAL_DESCRIPTIVE_CONTRAST_SUPPORTED"
            comparison["higher_than_native_q95_minus_native_analogue_same_continent"]=(
                _bootstrap_shared_cells(wet,analogue,
                  name=taxon+"|"+continent+"|"+climate,replicates=bootstrap)
            )
        result["hypotheses"][label]=comparison
    return result


def analyse_photo_orientation(wild_photo:pd.DataFrame,
                              *,bootstrap:int=BOOTSTRAPS)->dict:
    required={"obs_id","taxon_name","range_state","user_id","_cell","tdwg1",
              "orientation_angle_degrees_median","observed_on",*VARIABLES}
    if not required<=set(wild_photo):
        raise ValueError("Original photo/observer/original-pixel climate values incomplete")
    if wild_photo.obs_id.duplicated().any() or wild_photo.user_id.isna().any():
        raise ValueError("Source photo observation not unique or missing human original observer")
    if set(wild_photo.taxon_name)!=set(SPECIES):
        raise ValueError("Only the two exhaustively supported original source taxa")
    if set(wild_photo.range_state)!={"native","introduced"}:
        raise ValueError("Original WCVP range classes lost")
    good=(np.isfinite(wild_photo[list(VARIABLES)].to_numpy(float)).all(axis=1)
          &np.isfinite(wild_photo.orientation_angle_degrees_median.to_numpy(float))
          &wild_photo.orientation_angle_degrees_median.between(0,180).to_numpy())
    filtered=wild_photo.loc[good].copy()
    counted=dated_and_balanced(filtered)
    result={}
    for name,frame in counted.groupby("taxon_name",sort=True):
        native=frame.loc[frame.range_state.eq("native")].copy()
        introduced=frame.loc[frame.range_state.eq("introduced")].copy()
        threshold=climate_reference(native)
        climates=standardization(frame)
        output={
            "native_source_observation_count":len(native),
            "introduced_source_observation_count":len(introduced),
            "native_5deg_cell_count":int(native._cell.nunique()),
            "introduced_5deg_cell_count":int(introduced._cell.nunique()),
            "native_reference":threshold,
            "regions":{},
            "original_five_exposure_between_region_support":"annual BIO12 primary and BIO15 secondary, all source values retained without raster substitution",
        }
        for continent,region in introduced.groupby("tdwg1",sort=True):
            if pd.isna(continent) or int(continent) not in CONTINENTS:
                raise ValueError("Original introduced TDWG continent invalid")
            group=CONTINENTS[int(continent)]
            regional=compare_within_continent(region,threshold,taxon=name,
                                             continent=group,bootstrap=bootstrap)
            # A five-covariate linear slope is secondary to the direct
            # predeclared wet/dry vs analogue regional comparison.
            if len(region)>=100 and region._cell.nunique()>=15:
                fitted=cell_weighted_fit(region,climates)
                regional["secondary_continent_five_climate_rain_slope"]=fitted
            else:
                regional["secondary_continent_five_climate_rain_slope"]={
                    "status":"STOP_INSUFFICIENT_SOURCE_CELLS_FOR_CONDITIONAL_SLOPE"}
            output["regions"][group]=regional
        result[name]=output
    return {
        "schema_version":"v4_S2_verified_original_full_CHELSA_head_angle_rain_novelty_20261021_v1",
        "status":"DIRECT_CLIMATE_NOVELTY_AND_IMAGE_ORIENTATION_FINISHED_EXPLORATORY",
        "source_wild_photo_eligible_before_date_cell_dedup":len(wild_photo),
        "source_original_climate_and_date_eligible":len(filtered),
        "source_after_original_human_cell_balance":len(counted),
        "focal":result,
        "scientific_limits":[
            "Annual precipitation BIO12 is the primary hypothesis; BIO15 precipitation seasonality is secondary and was not independently supported in original global S1.",
            "The reference climate envelope is conditioned on source photographic angle-eligible sightings, not a full ecological niche.",
            "Wetter-than-native q95 and analogue comparisons are matched within each introduced continent, not the same genets before and after invasion.",
            "Weather conditions during flowering are unmeasured; original CHELSA 1981-2010 values are climate normals.",
            "Image vertical is not actual gravitational head inclination; photographic camera bias, individual identity, flowering stage and source genotype remain unresolved.",
            "Same-continent five-climate conditional linear slopes are exploratory; correlations among BIO1, BIO4, BIO12, BIO15, RSDS, month and sampling geography remain.",
            "Non-significant intervals are NOT evidence of zero response, niche conservatism, genetic fixation or rapid plasticity/evolution.",
            "Bootstrap resamples source 5-degree cells, not individual plant genet or independent introductions."
        ],
        "no_original_location_observation_photographer_identifiers_in_public_result":True
    }

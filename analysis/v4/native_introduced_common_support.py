"""S2 climate transfer falsification on prespecified COMMON SUPPORT.

The previously completed native→introduced comparison used all original-pixel
observations and reported negative out-of-domain prediction skill without a
common climate-support restriction. This supplementary analysis applies
*only* outcome-blind, source-climate 5%-95% interval overlap for BIO12 alone
and for ALL FIVE frozen predictors, separately for C. arvense and C. vulgare.

BIO12-only overlap is a 1-D exposure availability diagnostic. It does NOT
eliminate extrapolation in remaining four climates. The all-five rectangle
is a conservative MARGINAL overlap screen, not proof of true multivariate
density overlap; no genetic, adaptation, fitness or evolution claim follows.
Missing cache pixels stay missing, and empty/rank-deficient strata FAIL CLOSED.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_introduced_climate_transfer import (
    SOURCE_REPORT_SHA256, cell_weighted_fit, climate_common_support,
    dated_and_balanced, design, evaluate_focal, load_verified_focal,
    standardization, transfer_against_native_baseline,
)
from analysis.v4.native_cached_climate import physical_climate_for_focal, digest_file

REFERENCE_RUN=39539209624
REFERENCE_ARTIFACT=12317764638
REFERENCE_RESULTS_SHA256="a95e0eaa085cc4cbd55a8302f5d527e26e3800761f7bb043ae9e65bfecdbdb26"
MIN_CELLS=15
MIN_OBSERVERS=10
MIN_ROWS=30
QUANTILES=(.05,.95)
SECTOR_NAMES={1:"Europe",2:"Africa",3:"Asia temperate",4:"Asia tropical",
              5:"Australasia",6:"Pacific",7:"North America",8:"South America",
              9:"Antarctic"}


def overlap_rectangle(native:pd.DataFrame,introduced:pd.DataFrame,
                      variables:tuple[str,...])->dict:
    """Outcome-blind source-climate marginal 90% intersection.

    The two range-specific distributions have already passed original source
    WCVP established-status filters. No angle is consulted for thresholds,
    conditional on this previously photographic-eligible input frame.
    """
    if not variables or len(set(variables))!=len(variables):
        raise ValueError("A unique non-empty prespecified climate vector is required")
    if not set(variables)<=set(native) or not set(variables)<=set(introduced):
        raise ValueError("Frozen climate predictors unavailable")
    result={}
    for name in variables:
        n=native[name].to_numpy(float);i=introduced[name].to_numpy(float)
        if not np.isfinite(n).all() or not np.isfinite(i).all() or not len(n) or not len(i):
            raise ValueError("Do not compute overlap from missing climate or zero source rows")
        nq=np.quantile(n,QUANTILES);iq=np.quantile(i,QUANTILES)
        lo=max(float(nq[0]),float(iq[0]))
        hi=min(float(nq[1]),float(iq[1]))
        result[name]={
            "native_q05_q95":[float(nq[0]),float(nq[1])],
            "introduced_q05_q95":[float(iq[0]),float(iq[1])],
            "intersection":[lo,hi] if hi>lo else None
        }
    bounds_valid=all(v["intersection"] is not None for v in result.values())
    def within(data):
        if not bounds_valid:
            return np.zeros(len(data),bool)
        mask=np.ones(len(data),bool)
        for var,v in result.items():
            lo,hi=v["intersection"]
            mask &= data[var].to_numpy(float)>=lo
            mask &= data[var].to_numpy(float)<=hi
        return mask
    ma=within(native);mb=within(introduced)
    return {
        "predictors":list(variables),
        "limits":result,
        "climate_intersections_nonempty":bounds_valid,
        "native_source_rows":len(native),
        "introduced_source_rows":len(introduced),
        "native_rows_in_overlap":int(ma.sum()),
        "introduced_rows_in_overlap":int(mb.sum()),
        "native_fraction_in_overlap":float(ma.mean()),
        "introduced_fraction_in_overlap":float(mb.mean()),
        "no_response_value_entered_overlap_selection":True,
        "native_mask":ma,
        "introduced_mask":mb,
    }


def report_support(rectangle:dict)->dict:
    return {k:v for k,v in rectangle.items() if k not in {"native_mask","introduced_mask"}}


def _needs_sources(s:pd.DataFrame)->dict:
    return {
        "n_original_observer_cell_balanced":int(len(s)),
        "unique_original_observers":int(s.user_id.nunique()),
        "unique_original_5degree_cells":int(s._cell.nunique()),
        "original_TDWG_level1_regions":int(s.tdwg1.nunique()),
    }


def compare_on_support(native:pd.DataFrame,introduced:pd.DataFrame,
                       scales:dict,rectangle:dict,
                       *,compare_models:bool=True)->dict:
    n=native.loc[rectangle["native_mask"]].copy()
    i=introduced.loc[rectangle["introduced_mask"]].copy()
    out={
        "climate_overlap":report_support(rectangle),
        "native_support":_needs_sources(n),
        "introduced_support":_needs_sources(i)
    }
    viable=all(len(a)>=MIN_ROWS and a._cell.nunique()>=MIN_CELLS
               and a.user_id.nunique()>=MIN_OBSERVERS for a in (n,i))
    if not viable:
        out.update(
            status="STOP_INSUFFICIENT_CLIMATE_COMMON_SUPPORT",
            slope_difference=None,
            native_trained_transfer=None,
            reason="At least one range has too few observer-cell rows, independent 5-degree cells or observers after climate-only overlap"
        )
        return out
    if not compare_models:
        return dict(out,status="support_sufficient_not_refitted")
    slopes={name:cell_weighted_fit(s,scales) for name,s in (("native",n),("introduced",i))}
    out["refitted_models"]=slopes
    if any(x["status"]!="estimated_five_exposure_cell_cluster_WLS" for x in slopes.values()):
        out.update(
            status="STOP_FIVE_CLIMATE_DESIGN_RANK_OR_CONDITION",
            slope_difference=None,native_trained_transfer=None
        )
        return out
    p=slopes["native"];q=slopes["introduced"]
    delta=(q["annual_precip_degrees_per_100mm"]
           -p["annual_precip_degrees_per_100mm"])
    delta_se=np.hypot(q["annual_precip_robust_se_per_100mm"],
                      p["annual_precip_robust_se_per_100mm"])
    out.update(
        status="ESTIMATED_CONDITIONAL_MARGINAL_COMMON_SUPPORT_SENSITIVITY",
        slope_difference={
            "introduced_minus_native_BIO12_degrees_per_100mm":float(delta),
            "descriptive_cell_cluster_95pct_interval":[
                float(delta-1.96*delta_se),float(delta+1.96*delta_se)
            ],
            "independence_of_origins_or_photography_not_asserted":True,
            "not_equivalence_testing":True
        },
        native_trained_transfer=transfer_against_native_baseline(n,i,scales),
        support_does_not_remove_camera_technical_error=True
    )
    return out


def continent_sensitivity(native:pd.DataFrame,introduced:pd.DataFrame,
                          scales:dict)->dict:
    """Declare introduced continent-specific non-transfer before outcome tests."""
    out={}
    for continent,part in introduced.groupby("tdwg1",sort=True):
        if pd.isna(continent):
            raise ValueError("Unknown original introduced territory parent continent")
        key=SECTOR_NAMES.get(int(continent))
        if not key:
            raise ValueError("Unregistered TDWG level 1 continent")
        rect=overlap_rectangle(native,part,tuple(VARIABLES))
        # Preserve every territory's observed support and failure, even when
        # the original pixel cache has no rows in that continent.
        out[key]=compare_on_support(native,part,scales,rect)
        out[key]["introduced_region_original_observation_count"]=int(len(part))
    return out


def re_audit(private_focal:pd.DataFrame,report:dict)->dict:
    if set(private_focal.taxon_name)!={"Cirsium arvense","Cirsium vulgare"}:
        raise ValueError("Only preregistered focal species allowed")
    previous=report["focal_species"]
    if set(previous)!=set(private_focal.taxon_name):
        raise ValueError("Original fitted species set differs from same source")
    response_key="orientation_angle_degrees_median"
    results={}
    for name,frame in private_focal.groupby("taxon_name",sort=True):
        frame=frame.copy()
        frame["_complete"]=np.isfinite(frame[list(VARIABLES)]).all(axis=1)
        original=previous[name]
        if not frame.loc[frame.range_state.eq("native")].shape[0]==original["range_support"]["native"]["source_orientation_observations"]:
            raise ValueError("Native source and original fitted analysis disagree")
        if not frame.loc[frame.range_state.eq("introduced")].shape[0]==original["range_support"]["introduced"]["source_orientation_observations"]:
            raise ValueError("Introduced source and original fitted analysis disagree")
        f=frame.loc[frame._complete].copy()
        scales=standardization(f)
        base=dated_and_balanced(f)
        native=base.loc[base.range_state.eq("native")]
        introduced=base.loc[base.range_state.eq("introduced")]
        original_overlap=previous[name]["joint_climatic_support"]
        if len(native)!=original["range_support"]["native"]["source_original_observer_cell_balanced_rows"]:
            raise ValueError("Native balanced cohort changed")
        if len(introduced)!=original["range_support"]["introduced"]["source_original_observer_cell_balanced_rows"]:
            raise ValueError("Introduced balanced cohort changed")
        rain=overlap_rectangle(native,introduced,("chelsa_bio12",))
        joint=overlap_rectangle(native,introduced,tuple(VARIABLES))
        original_beta_n=previous[name]["estimated_within_range_slope"]["native"]["annual_precip_degrees_per_100mm"]
        original_beta_i=previous[name]["estimated_within_range_slope"]["introduced"]["annual_precip_degrees_per_100mm"]
        if abs(cell_weighted_fit(native,scales)["annual_precip_degrees_per_100mm"]-original_beta_n)>1e-8:
            raise ValueError("Original native five-climate regression not reproduced")
        if abs(cell_weighted_fit(introduced,scales)["annual_precip_degrees_per_100mm"]-original_beta_i)>1e-8:
            raise ValueError("Original introduced five-climate regression not reproduced")
        for var in VARIABLES:
            bounds=original_overlap[var]["joint_90pct_central_overlap"]
            got=overlap_rectangle(native,introduced,(var,))["limits"][var]["intersection"]
            if bounds is None and got is None:
                continue
            if bounds is None or got is None or not np.allclose(bounds,got,atol=1e-9):
                raise ValueError("Original univariate overlap distribution differs")
        results[name]={
            "previous_all_exposure_transfer":{
                "native_slope":original_beta_n,
                "introduced_slope":original_beta_i,
                "native_to_introduced_prediction_skill_vs_native_constant":previous[name]["original_native_trained_introduced_transfer"]["native_trained_model_to_introduced_skill_vs_native_constant"],
                "out_of_range_extrapolation_not_excluded":True
            },
            "one_exposure_BIO12_5_95pct_overlap_diagnostic":compare_on_support(native,introduced,scales,rain),
            "all_five_exposure_joint_rectangle_5_95pct_overlap":compare_on_support(native,introduced,scales,joint),
            "introduced_continent_sensitivity":continent_sensitivity(native,introduced,scales),
            "design_scale":"Exactly original per-species SD and the same BIO1 BIO4 BIO12 BIO15 RSDS plus month/year",
        }
    return {
        "schema_version":"v4_S2_native_introduced_common_support_falsification_20261021_v1",
        "status":"SOURCE_BACKED_OUTCOME_BLIND_COMMON_SUPPORT_AND_CONDITIONAL_TRANSFERS",
        "source_original_transfer_action":REFERENCE_RUN,
        "source_original_transfer_artifact_id":REFERENCE_ARTIFACT,
        "source_original_transfer_json_sha256":REFERENCE_RESULTS_SHA256,
        "source_status_json_sha256":SOURCE_REPORT_SHA256,
        "focal_species":results,
        "reported_units":"camera-image degrees per 100 mm annual precipitation, not physical head tilt",
        "limitations":[
            "Climatic bounds are 5th–95th intervals computed on existing orientation-eligible observer-balanced photos: this is NOT an unbiased flowering population.",
            "Individual climate ranges alone do not guarantee density common support in 5D; joint rectangular intersection is necessary but not sufficient.",
            "BIO12-only region does not control extrapolation of BIO1/BIO4/BIO15/RSDS; it is a rainfall support diagnostic, not proof of transfer.",
            "Different source continents have independent possible histories, no introduction year or genotype inferred from these pictures.",
            "A narrow confidence interval still cannot distinguish evolutionary fixation from plasticity, sampling, technical photo angle and founders.",
            "Conditional cell-cluster intervals do not remove phylogenetic, camera angle, individual or geographic provenance dependence.",
            "Model is run on recorded original photo response, not gravity-referenced or experimentally manipulated floral heads.",
        ],
        "no_photo_or_location_or_user_ids_published":True
    }


def execute(focal:Path,source_receipt:Path,baseline:Path,
            climate_root:Path,metadata:Path,output:Path):
    # The previous raw source and frozen 5-climate values must be regenerated
    # on the SAME private runner; a run-specific JSON hash verifies estimands.
    from analysis.v4.native_introduced_climate_transfer import (
        load_verified_focal
    )
    f,_=load_verified_focal(focal,source_receipt)
    vals,audit=physical_climate_for_focal(f,climate_root,metadata)
    for var,data in vals.items():
        f[var]=data
    payload=baseline.read_bytes()
    if hashlib.sha256(payload).hexdigest()!=REFERENCE_RESULTS_SHA256:
        raise ValueError("Previous native-introduced numerical result hash changed")
    original=json.loads(payload)
    result=re_audit(f,original)
    result["verified_CHELSA_source_grid_hashes"]={
        name:audit[name]["source_raster_sha256"] for name in VARIABLES
    }
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_native_introduced_common_support.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    brief={
        taxon:{
            "one_rain_bio12_only_support":q["one_exposure_BIO12_5_95pct_overlap_diagnostic"]["status"],
            "all_five_climate_support":q["all_five_exposure_joint_rectangle_5_95pct_overlap"]["status"],
            "all_five_retained":q["all_five_exposure_joint_rectangle_5_95pct_overlap"]["climate_overlap"]
                ["native_rows_in_overlap"]
        } for taxon,q in result["focal_species"].items()
    }
    print(json.dumps({"status":result["status"],"focal":brief},ensure_ascii=False),flush=True)
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    for field in ("focal","source_receipt","baseline","climate_root","metadata","output"):
        p.add_argument("--"+field.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    execute(a.focal,a.source_receipt,a.baseline,a.climate_root,a.metadata,a.output)


if __name__=="__main__":
    main()

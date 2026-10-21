"""Pair climatic novelty with capitulum photo orientation on the SAME records.

A climate-only native envelope assigns each introduced record to one of four
predeclared mutually exclusive groups:
  (1) wetter than native BIO12 q95, (2) drier than native BIO12 q05,
  (3) other-climate novel but native-analogue rain, (4) analogue all five.
All bounds use native photo-eligible source values, with equal 5-degree-cell
weights and NO orientation response. This is an observational *realized
occupancy vs photographic trait* comparison, NOT evolution, phenotypic change
through time, plasticity or causal rain shielding.

Outcomes: within introduced records compare actual photo angle and predicted
native five-climate+month/year angle for wet/dry/other-novel vs analogue,
using cell-cluster intervals. Never publish photo/user/site/observation IDs.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_introduced_climate_transfer import (
    cell_weighted_fit, design, MIN_CELLS, LIMIT_CONDITION,
)

CLASSES=("rain_wetter_than_native","rain_drier_than_native",
         "other_climate_novel","native_climate_analogue")
REF="native_climate_analogue"
MIN_CATEGORY_ROWS=30
MIN_CATEGORY_CELLS=5
MIN_CATEGORY_OBSERVERS=10
SEED=20261021


def weighted_quantiles(values:np.ndarray,weights:np.ndarray,quantiles=(.05,.95))->list[float]:
    values=np.asarray(values,float)
    weights=np.asarray(weights,float)
    if (len(values)==0 or len(weights)!=len(values) or
        not np.isfinite(values).all() or not np.isfinite(weights).all() or
        np.any(weights<=0)):
        raise ValueError("Native climatic reference must have finite cell weights")
    order=np.argsort(values,kind="mergesort")
    csum=(np.cumsum(weights[order])-.5*weights[order])/weights.sum()
    return [float(np.interp(q,csum,values[order],
                             left=values[order][0],right=values[order][-1]))
            for q in quantiles]


def original_cell_weights(frame:pd.DataFrame)->np.ndarray:
    if "_cell" not in frame or frame["_cell"].isna().any():
        raise ValueError("Original 5-degree cell assignments not present")
    n=frame.groupby("_cell",sort=False).obs_id.transform("size").to_numpy(float)
    if not np.isfinite(n).all() or (n<=0).any():
        raise ValueError("Original cell weights changed")
    return 1/n


def climate_only_groups(native:pd.DataFrame,introduced:pd.DataFrame)->tuple[np.ndarray,dict]:
    if not set(VARIABLES)<=set(native) or not set(VARIABLES)<=set(introduced):
        raise ValueError("Five frozen climate axes unavailable")
    if len(native)<100 or len(introduced)<100:
        raise ValueError("Too few original orientation-eligible climate records")
    weights=original_cell_weights(native)
    bounds={}
    outside=np.zeros(len(introduced),bool)
    for v in VARIABLES:
        vn=native[v].to_numpy(float);vi=introduced[v].to_numpy(float)
        if not np.isfinite(vn).all() or not np.isfinite(vi).all():
            raise ValueError("Missing original CHELSA pixel values")
        low,high=weighted_quantiles(vn,weights)
        if low>=high:
            raise ValueError("Degenerate native climate envelope")
        bounds[v]={"native_cell_equal_q05":low,"native_cell_equal_q95":high}
        outside|=(vi<low)|(vi>high)
    rain=introduced["chelsa_bio12"].to_numpy(float)
    lo=bounds["chelsa_bio12"]["native_cell_equal_q05"]
    hi=bounds["chelsa_bio12"]["native_cell_equal_q95"]
    labels=np.full(len(introduced),REF,dtype=object)
    labels[outside]="other_climate_novel"
    labels[rain<lo]="rain_drier_than_native"
    labels[rain>hi]="rain_wetter_than_native"
    if not set(labels)<=set(CLASSES):
        raise AssertionError("Four climate-only groups did not partition introduced source")
    return labels,{
        "climate_axes":list(VARIABLES),
        "thresholds":"Native geographic-cell equal-weighted fifth and ninety-fifth percentiles",
        "native_reference_bounds":bounds,
        "rule_precedence":["rain wetter","rain drier",
                           "other climate novel but rainfall in native bounds",
                           "analogue all five native marginal intervals"],
        "climate_group_computed_without_orientation_values":True,
        "not_proof_of_joint_climatic_support_or_fundamental_niche":True,
    }


def native_camera_model(native:pd.DataFrame,scales:dict)->tuple[np.ndarray,dict]:
    """Same 5-predictor cell-equal WLS used by the S2 initial transfer fit."""
    fitted=cell_weighted_fit(native,scales)
    if fitted["status"]!="estimated_five_exposure_cell_cluster_WLS":
        raise ValueError("Native climate-photo-angle fit not identifiable")
    x=design(native,scales)
    w=original_cell_weights(native)
    xtwx=x.T@(x*w[:,None])
    if np.linalg.matrix_rank(xtwx)<x.shape[1] or np.linalg.cond(xtwx)>LIMIT_CONDITION:
        raise ValueError("Native original-angle model lost its design rank")
    beta=np.linalg.solve(xtwx,x.T@(w*native.orientation_angle_degrees_median.to_numpy(float)))
    annual=float(beta[list(VARIABLES).index("chelsa_bio12")+1]
                 *100/scales["chelsa_bio12"]["sd"])
    if abs(annual-fitted["annual_precip_degrees_per_100mm"])>1e-8:
        raise ValueError("Native fitted climate effect no longer matches frozen model")
    return beta,fitted


def _group_cell_summary(frame:pd.DataFrame,category:str)->dict:
    eligible=frame.loc[frame["_novelty"].eq(category)]
    n=len(eligible)
    cells=int(eligible["_cell"].nunique())
    observers=int(eligible.user_id.nunique())
    out={
        "photo_observations":n,"original_5degree_cells":cells,
        "distinct_original_observers":observers,
        "satisfies_predeclared_descriptive_support":bool(
            n>=MIN_CATEGORY_ROWS and cells>=MIN_CATEGORY_CELLS
            and observers>=MIN_CATEGORY_OBSERVERS),
    }
    if not out["satisfies_predeclared_descriptive_support"]:
        out["status"]="NOT_ESTIMABLE_INDEPENDENT_CELLS_OR_OBSERVERS"
        out["observed_mean_image_angle_degrees"]=None
        out["native_model_predicted_mean_image_angle_degrees"]=None
        out["observed_minus_native_model_predicted_degrees"]=None
        out["observed_mean_BIO12_mm_per_year"]=None
        return out
    means=eligible.groupby("_cell",sort=True)[
        ["_y","_pred","_resid","chelsa_bio12"]
    ].mean()
    out.update({
        "status":"CELL_EQUAL_DESCRIPTIVE_OBSERVED_ANGLE_AND_NATIVE_PREDICTION",
        "observed_mean_image_angle_degrees":float(means["_y"].mean()),
        "native_model_predicted_mean_image_angle_degrees":float(means["_pred"].mean()),
        "observed_minus_native_model_predicted_degrees":float(means["_resid"].mean()),
        "observed_mean_BIO12_mm_per_year":float(means["chelsa_bio12"].mean()),
    })
    return out


def _cell_cluster_contrasts(frame:pd.DataFrame,groups:dict,*,bootstrap:int)->dict:
    if bootstrap<100:
        raise ValueError("At least 100 prespecified cell bootstrap resamples required")
    bycell=frame.groupby(["_cell","_novelty"],sort=True)[["_y","_pred","_resid"]].mean()
    unique=sorted(frame._cell.unique())
    index={cell:i for i,cell in enumerate(unique)}
    cells=len(unique)
    presence=np.zeros((cells,len(CLASSES)))
    metrics={m:np.zeros((cells,len(CLASSES))) for m in ("_y","_pred","_resid")}
    for (cell,label),r in bycell.iterrows():
        j=CLASSES.index(label)
        k=index[cell]
        presence[k,j]=1
        for key in metrics:
            metrics[key][k,j]=float(r[key])
    rng=np.random.default_rng(SEED)
    # Bootstrap cells with replacement, preserving multiple climate groups
    # occurring in one original geographical cell.
    draw=rng.multinomial(cells,[1/cells]*cells,size=bootstrap).astype(float)
    counts=draw@presence
    out={}
    refix=CLASSES.index(REF)
    for label in CLASSES:
        if label==REF:
            continue
        if not (groups[label]["satisfies_predeclared_descriptive_support"] and
                groups[REF]["satisfies_predeclared_descriptive_support"]):
            out[label]={"status":"STOP_INSUFFICIENT_CATEGORY_OR_REFERENCE_SUPPORT"}
            continue
        j=CLASSES.index(label)
        shared=(counts[:,j]>0)&(counts[:,refix]>0)
        if int(shared.sum())<.80*bootstrap:
            out[label]={"status":"STOP_CLUSTER_RESAMPLE_TOO_FEW_VALID_SOURCE_CELLS",
                        "bootstrap_valid_replicates":int(shared.sum())}
            continue
        result={
            "status":"INTRODUCED_CLIMATE_CONTRAST_ESTIMATED_CELL_CLUSTER_CONDITIONAL",
            "bootstrap_clusters":"original distinct geographic 5-degree cells",
            "bootstrap_replicates_requested":bootstrap,
            "bootstrap_replicates_valid":int(shared.sum()),
            "reference_class":REF,
            "contrast_class":label,
        }
        for key,label2 in (
            ("_y","observed_image_angle"),
            ("_pred","native_model_expected_image_angle"),
            ("_resid","observed_minus_expected_residual")):
            block=draw@metrics[key]
            deltas=(block[shared,j]/counts[shared,j]
                    -block[shared,refix]/counts[shared,refix])
            point=(groups[label][{
                "_y":"observed_mean_image_angle_degrees",
                "_pred":"native_model_predicted_mean_image_angle_degrees",
                "_resid":"observed_minus_native_model_predicted_degrees"}[key]]
                -groups[REF][{
                "_y":"observed_mean_image_angle_degrees",
                "_pred":"native_model_predicted_mean_image_angle_degrees",
                "_resid":"observed_minus_native_model_predicted_degrees"}[key]])
            result[label2+"_difference_degrees"]=float(point)
            result[label2+"_cell_cluster_q025_q50_q975"]=[
                float(z) for z in np.quantile(deltas,[.025,.5,.975])]
        # A rain shelter direction claim is permitted only conditionally and
        # only if actual observed (not model expected) angle contrast has
        # a cell bootstrap interval wholly in the registered direction.
        q=result["observed_image_angle_cell_cluster_q025_q50_q975"]
        direction=1 if label=="rain_wetter_than_native" else -1 if label=="rain_drier_than_native" else 0
        result["image_orientation_direction_consistent_with_rain_prediction"]=(
            bool(q[0]>0) if direction>0 else bool(q[2]<0) if direction<0 else None
        )
        result["not_gravity_or_functional_protection_evidence"]=True
        out[label]=result
    return out


def analyse_niche_trait_pair(native:pd.DataFrame,introduced:pd.DataFrame,
                             scales:dict,*,bootstrap:int=500)->dict:
    """New climatic regime and photo trait change evaluated together per species."""
    if not set(native.range_state)=={"native"} or not set(introduced.range_state)=={"introduced"}:
        raise ValueError("Focal source introduced/native WCVP status changed")
    if len(native)<100 or len(introduced)<100:
        raise ValueError("Source sample insufficient")
    for d in (native,introduced):
        if not d.obs_id.is_unique or d.user_id.isna().any():
            raise ValueError("Non-independent missing original observer/observation")
        if not d.orientation_angle_degrees_median.between(0,180).all():
            raise ValueError("Nonfinite or invalid archived photo orientation")
        if not d._cell.notna().all():
            raise ValueError("No true source 5-degree region assignments")
    labels,classification=climate_only_groups(native,introduced)
    beta,fit=native_camera_model(native,scales)
    introduced=introduced.copy()
    introduced["_novelty"]=labels
    introduced["_y"]=introduced.orientation_angle_degrees_median.astype(float)
    introduced["_pred"]=design(introduced,scales)@beta
    introduced["_resid"]=introduced["_y"]-introduced["_pred"]
    if not np.isfinite(introduced[["_y","_pred","_resid"]]).all().all():
        raise ValueError("Original camera response or projection nonfinite")
    groups={label:_group_cell_summary(introduced,label) for label in CLASSES}
    if sum(groups[label]["photo_observations"] for label in CLASSES)!=len(introduced):
        raise AssertionError("Climate bins no longer partition original introduced photos")
    contrast=_cell_cluster_contrasts(introduced,groups,bootstrap=bootstrap)
    return {
        "status":"EXECUTED_CLIMATE_NOVELTY_X_PHOTO_ANGLE_SCREEN_NOT_CAUSAL_TRAIT_EVOLUTION",
        "source_eligible_native_observations":len(native),
        "source_eligible_introduced_observations":len(introduced),
        "native_source_distinct_cells":int(native._cell.nunique()),
        "introduced_source_distinct_cells":int(introduced._cell.nunique()),
        "native_original_five_exposure_BIO12_slope_deg_per_100mm":
            fit["annual_precip_degrees_per_100mm"],
        "native_photo_model_prespecified_exposures":list(VARIABLES),
        "native_photo_model_adjustments":["observation_month","observation_year"],
        "introduced_four_climate_categories":groups,
        "categories_vs_climate_analogue":contrast,
        "native_climate_reference":classification,
        "technical_and_causal_boundaries":[
            "Climate group assignment uses only source climate, never observed photo angle.",
            "Difference in introduced photo angle is between geographic groups and not longitudinal change of the same plant/head.",
            "Source camera tilt and flower stage may differ by region; angle is image-vertical rather than gravity-vertical.",
            "A native-only fitted climate photo model extrapolates outside the native BIO12/other covariate envelope.",
            "Image measuring/dropout selection, origin genotype and founder effects can explain group contrasts.",
            "Cell bootstrap cannot identify independent genets, ecological protection, local adaptation or evolutionary transition time.",
            "No inference that an absent wet/dry category has a true zero ecological effect."
        ]
    }

"""Direct S2 comparison of original continuous head-photo angles: native vs ALL introduced.

This is the main S2 descriptive range contrast, not the q75/q95 wet-edge or
6-band decomposition. No angle re-binning: weighted empirical quantiles,
central spreads, extreme tails, and weighted 1-D Earth mover's distance.
Confidence sensitivity uses original geographically clustered 5-degree cells,
not photograph/individual bootstrap. No original IDs or points leave runner.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
from analysis.v4.native_climate_trait_coupling import (
    original_cell_weights, weighted_quantiles
)

SPECIES = ("Cirsium arvense", "Cirsium vulgare")
ANGLE = "orientation_angle_degrees_median"
QUANTILES = (.05,.10,.25,.50,.75,.90,.95)
QUANTILE_LABELS = ("Q05","Q10","Q25","Q50","Q75","Q90","Q95")
TAIL_BOUNDS = (30.,60.,90.,120.,150.)
TAIL_LABELS = ("ge30","ge60","ge90","ge120","ge150")
BOOTS = 299
SEED = 20261021
EXPECTED_ORIGINAL_PHOTO_ELIGIBLE = 173737
EXPECTED_BALANCED_COUNTS = {
    "Cirsium arvense":{"native":21870,"introduced":25657},
    "Cirsium vulgare":{"native":16246,"introduced":30287}
}


def _validate_frame(frame:pd.DataFrame) -> None:
    required = {"obs_id","user_id","taxon_name","range_state","latitude",
                "longitude","observed_on",ANGLE}
    if not required <= set(frame):
        raise ValueError("Authentic original image-angle source columns missing")
    if not frame.obs_id.is_unique or frame.user_id.isna().any():
        raise ValueError("Original source observations not unique or observer missing")
    if set(frame.taxon_name) != set(SPECIES) or set(frame.range_state) != {
        "native","introduced"
    }:
        raise ValueError("Original two species native/explicit introduced source changed")
    angle = frame[ANGLE].to_numpy(float)
    if not np.isfinite(angle).all() or np.any((angle<0)|(angle>180)):
        raise ValueError("Invalid archived image-vertical angle")


def _summary(frame:pd.DataFrame, weight:np.ndarray) -> dict:
    y=frame[ANGLE].to_numpy(float)
    weights=np.asarray(weight,float)
    if (not len(y) or len(weights)!=len(y) or
        not np.isfinite(weights).all() or np.any(weights<0) or weights.sum()<=0):
        raise ValueError("Invalid positive weighted angle sample")
    mask=weights>0
    if not mask.any():
        raise ValueError("No support remains after geographical weighting")
    yy=y[mask]
    ww=weights[mask]/weights[mask].sum()
    q=weighted_quantiles(yy,ww,quantiles=QUANTILES)
    mean=float(ww@yy)
    variance=float(ww@np.square(yy-mean))
    ge={k:float(ww[yy>=threshold].sum())
        for k,threshold in zip(TAIL_LABELS,TAIL_BOUNDS)}
    bins = [
        float(ww[(yy>=lo)&(yy<hi)].sum())
        for lo,hi in ((0.,30.),(30.,60.),(60.,90.),(90.,120.),(120.,150.))
    ]+[float(ww[yy>=150].sum())]
    if not np.isclose(sum(bins),1,atol=1e-12):
        raise ValueError("Camera image-angle bins do not sum to one")
    return {
        "quantiles_degrees":dict(zip(QUANTILE_LABELS,[float(z) for z in q])),
        "mean_degrees":mean,
        "sd_degrees":float(np.sqrt(max(variance,0.))),
        "IQR_Q75_minus_Q25_degrees":float(q[4]-q[2]),
        "central_80pct_Q90_minus_Q10_degrees":float(q[5]-q[1]),
        "central_90pct_Q95_minus_Q05_degrees":float(q[6]-q[0]),
        "threshold_exceedance_fraction":ge,
        "nonoverlapping_angle_regions":dict(zip(
            ("0_to_lt30","30_to_lt60","60_to_lt90",
             "90_to_lt120","120_to_lt150","150_to_180"),bins))
    }


def _delta(native:dict, introduced:dict) ->dict:
    return {
        "quantiles_deg":{key:introduced["quantiles_degrees"][key]
                           -native["quantiles_degrees"][key]
                         for key in QUANTILE_LABELS},
        "mean_deg":introduced["mean_degrees"]-native["mean_degrees"],
        "sd_deg":introduced["sd_degrees"]-native["sd_degrees"],
        "IQR_deg":introduced["IQR_Q75_minus_Q25_degrees"]
                  -native["IQR_Q75_minus_Q25_degrees"],
        "central80_width_deg":introduced["central_80pct_Q90_minus_Q10_degrees"]
                              -native["central_80pct_Q90_minus_Q10_degrees"],
        "central90_width_deg":introduced["central_90pct_Q95_minus_Q05_degrees"]
                              -native["central_90pct_Q95_minus_Q05_degrees"],
        "threshold_exceedance_fraction":{
            key:introduced["threshold_exceedance_fraction"][key]
                 -native["threshold_exceedance_fraction"][key]
            for key in TAIL_LABELS
        }
    }


def _continuous_w1(left:pd.DataFrame,right:pd.DataFrame,
                   wl:np.ndarray,wr:np.ndarray)->float:
    al=left[ANGLE].to_numpy(float)
    ar=right[ANGLE].to_numpy(float)
    a=wl>0
    b=wr>0
    return float(wasserstein_distance(al[a],ar[b],
                   u_weights=wl[a],v_weights=wr[b]))


def compare(native:pd.DataFrame,introduced:pd.DataFrame,
            *,seed:int,replicates:int=BOOTS)->dict:
    """Full source phototype from two WCVP range states, no rain strata."""
    if min(len(native),len(introduced))<100:
        return {"status":"STOP_TOO_FEW_IMAGES"}
    if min(native._cell.nunique(),introduced._cell.nunique())<15:
        return {"status":"STOP_TOO_FEW_SPATIAL_BLOCKS"}
    wn=original_cell_weights(native)
    wi=original_cell_weights(introduced)
    s0=_summary(native,wn)
    s1=_summary(introduced,wi)
    observed=_delta(s0,s1)
    observed["W1_continuous_original_angles_degrees"]=_continuous_w1(
        native,introduced,wn,wi)
    # Apply exactly the SAME 5-degree-cell resample to any cell present in
    # both geographical establishment groups (coarse cells can overlap).
    cell_union, index=np.unique(
        np.concatenate((native._cell.to_numpy(int),
                        introduced._cell.to_numpy(int))),
        return_inverse=True
    )
    j=len(native)
    ni=index[:j]
    ii=index[j:]
    rng=np.random.default_rng(seed)
    stats=[]
    for k in range(replicates):
        multiplicity=np.bincount(
            rng.integers(len(cell_union),size=len(cell_union)),
            minlength=len(cell_union))
        wa=wn*multiplicity[ni]
        wb=wi*multiplicity[ii]
        if (wa.sum()<=0 or wb.sum()<=0 or
            np.count_nonzero(np.bincount(ni,weights=wa))<10 or
            np.count_nonzero(np.bincount(ii,weights=wb))<10):
            continue
        a=_summary(native,wa)
        b=_summary(introduced,wb)
        d=_delta(a,b)
        d["W1_continuous_original_angles_degrees"]=_continuous_w1(
            native,introduced,wa,wb)
        stats.append(d)
    if len(stats)<.90*replicates:
        return {"status":"STOP_UNSTABLE_BLOCK_BOOTSTRAP",
                "replicates_valid":len(stats)}
    def interval(values):
        return [float(x) for x in np.quantile(values,(.025,.50,.975))]
    ci={
        "quantiles_deg":{key:interval([x["quantiles_deg"][key] for x in stats])
                         for key in QUANTILE_LABELS},
        "IQR_deg":interval([x["IQR_deg"] for x in stats]),
        "central80_width_deg":interval([x["central80_width_deg"] for x in stats]),
        "central90_width_deg":interval([x["central90_width_deg"] for x in stats]),
        "mean_deg":interval([x["mean_deg"] for x in stats]),
        "sd_deg":interval([x["sd_deg"] for x in stats]),
        "threshold_exceedance_fraction":{
            key:interval([x["threshold_exceedance_fraction"][key] for x in stats])
            for key in TAIL_LABELS},
        "W1_continuous_original_angles_degrees":
            interval([x["W1_continuous_original_angles_degrees"] for x in stats])
    }
    return {
        "status":"EXPLORATORY_SOURCE_IMAGE_CONTINUOUS_ANGLE_DISTRIBUTION_COMPARED",
        "native":s0,
        "all_introduced":s1,
        "introduced_minus_native":observed,
        "fixed_original_cell_bootstrap_sensitivity_95pct":ci,
        "original_native_angle_observations":len(native),
        "original_introduced_angle_observations":len(introduced),
        "native_5deg_cells":int(native._cell.nunique()),
        "all_introduced_5deg_cells":int(introduced._cell.nunique()),
        "native_original_observers":int(native.user_id.nunique()),
        "introduced_original_observers":int(introduced.user_id.nunique()),
        "spatial_bootstrap":{
            "unit":"union original 5degree geographical cell",
            "draws_requested":replicates,
            "valid_draws":len(stats),
            "not_true_plant_genet_or_image_camera_error_interval":True,
            "W1_nonnegative_interval_not_test_of_distribution_equality":True,
        }
    }


def evaluate(photo:pd.DataFrame,*,bootstrap:int=BOOTS)->dict:
    _validate_frame(photo)
    if len(photo)!=EXPECTED_ORIGINAL_PHOTO_ELIGIBLE:
        raise ValueError("Archived two-species photo-eligible source changed")
    balanced=dated_and_balanced(photo)
    counts={
        name:{s:int(len(g)) for s,g in frame.groupby("range_state")}
        for name,frame in balanced.groupby("taxon_name")
    }
    if counts!=EXPECTED_BALANCED_COUNTS:
        raise ValueError("Original photo-angle observer×cell range cohort changed")
    result={}
    for name in SPECIES:
        frame=balanced.loc[balanced.taxon_name.eq(name)]
        a=frame.loc[frame.range_state.eq("native")]
        b=frame.loc[frame.range_state.eq("introduced")]
        species_seed=SEED^int.from_bytes(hashlib.sha256(name.encode()).digest()[:8],"little")
        result[name]=compare(a,b,seed=species_seed,replicates=bootstrap)
    return {
        "schema":"azami_S2_native_vs_all_introduced_full_continuous_image_angle_distribution_v1",
        "status":"ORIGINAL_SOURCE_TWO_SPECIES_DIRECT_7_QUANTILES_TAILS_SPREAD_AND_W1",
        "source_original_eligible_angles":len(photo),
        "balanced_observer_cell_angle_observations":len(balanced),
        "comparison":"native vs ALL explicitly WCVP introduced; not continent matched and no rain-banding",
        "angle_reference":"EXIF-adjusted image vertical (0 up, 90 horizontal, 180 down), not actual gravity",
        "quantile_probabilities":list(QUANTILES),
        "distance_is_degrees_on_original_continuous_angle_not_binned":True,
        "species":result,
        "public_aggregate_no_original_obs_photographer_or_coordinates":True,
        "limits":[
            "Two original WCVP source populations, not before/after individual plants",
            "Photographer, plant developmental phase and camera tilt are not separately identified",
            "High observed within-species residual variance is not proven genetic polymorphism",
            "Cell sensitivity intervals conditional on measured image-angle and sample, not genetic evolution",
            "W1 always nonnegative; reported bootstrap interval not a significance test",
            "Rainfall occupancy medians arise from a different original non-trait-selected cohort"
        ]
    }


def save(photo:pd.DataFrame,output:Path,*,bootstrap:int=BOOTS)->dict:
    result=evaluate(photo,bootstrap=bootstrap)
    output.mkdir(parents=True,exist_ok=True)
    target=output/"S2_native_vs_all_introduced_full_continuous_orientation.json"
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                      encoding="utf-8")
    return result

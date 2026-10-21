"""S2: start with observed full distributions, not preselected difference indices.

Compare WCVP-native versus ALL explicitly introduced occurrences for two Cirsium
species. Climate BIO12 uses original non-trait-selected source observations.
Image orientation uses original angle-qualified observations with identical
5-degree cell/photographer balancing. Store entire empirical distribution
curves on a shared, pooled, equal-range-probability *display* grid.

The display grid is NOT a chosen response metric, rain threshold, niche width,
ecological hypothesis test or significance test. The user decides which
features matter AFTER inspecting these unconditional, uncalibrated curves.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_climate_trait_coupling import original_cell_weights, weighted_quantiles
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced

SPECIES=("Cirsium arvense","Cirsium vulgare")
RAIN="chelsa_bio12"
ANGLE="orientation_angle_degrees_median"
DISPLAY_RESOLUTION=401
SPATIAL_BOOTSTRAPS=399
MIN_VALID_FRACTION=.90


def _validate(original:pd.DataFrame,photo:pd.DataFrame)->None:
    for frame in (original,photo):
        needed={"obs_id","user_id","taxon_name","range_state","_cell"}
        if not needed<=set(frame):
            raise ValueError("Missing original source identity for full distribution")
        if frame.obs_id.duplicated().any() or frame.user_id.isna().any():
            raise ValueError("Original source identity not unique or missing observer")
        if set(frame.taxon_name)!=set(SPECIES):
            raise ValueError("Unexpected original Cirsium source species")
        if set(frame.range_state)!={"native","introduced"}:
            raise ValueError("Invalid original WCVP native/introduced classification")
    if not set(VARIABLES)<=set(original) or ANGLE not in photo:
        raise ValueError("Missing original CHELSA source climate or photo orientation")
    if len(original)<1000 or len(photo)<1000:
        raise ValueError("Full original wild and photo source not available")


def _weights(frame:pd.DataFrame)->np.ndarray:
    w=original_cell_weights(frame)
    return w/w.sum()


def _pooled_shared_grid(a:np.ndarray,b:np.ndarray,
                        wa:np.ndarray,wb:np.ndarray)->np.ndarray:
    """Pooled mixture ranks supply plotting knots, never bin/score boundaries."""
    pooled=np.concatenate([a,b])
    mixture_weight=np.concatenate([wa/2,wb/2])
    grid=weighted_quantiles(
        pooled,mixture_weight,
        quantiles=np.linspace(0,1,DISPLAY_RESOLUTION)
    )
    # Include every source min/max and never truncate distribution tails.
    grid=np.unique(np.r_[min(a.min(),b.min()),grid,max(a.max(),b.max())])
    if len(grid)<2 or not np.isfinite(grid).all():
        raise ValueError("Degenerate pooled empirical distribution")
    return grid


def _ecdf(values:np.ndarray,weights:np.ndarray,grid:np.ndarray)->np.ndarray:
    order=np.argsort(values,kind="stable")
    ranks=np.searchsorted(values[order],grid,side="right")
    prefix=np.r_[0.,np.cumsum(weights[order])]
    cdf=prefix[ranks]
    cdf[-1]=1.0
    return cdf


def _cell_empirical_curves(frame:pd.DataFrame,column:str,grid:np.ndarray
                           )->tuple[np.ndarray,np.ndarray]:
    """Each geographic cell contributes total mass 1 before area balancing.

    Preaggregate full empirical CDF per original 5-degree cell and use exact
    same display knots for all spatial draws. No rain/angle cutpoints chosen.
    """
    ids,inverse=np.unique(frame._cell.to_numpy(int),return_inverse=True)
    if len(ids)<15:
        raise ValueError("Insufficient independent regional cells for uncertainty")
    counts=np.bincount(inverse)
    values=frame[column].to_numpy(float)
    indices=np.searchsorted(grid,values,side="left")
    if np.any(indices>=len(grid)):
        raise ValueError("Pooled full-distribution grid excluded source extrema")
    mass=np.zeros((len(ids),len(grid)),dtype=float)
    np.add.at(mass,(inverse,indices),1.0/counts[inverse])
    curves=np.cumsum(mass,axis=1)
    if not np.allclose(curves[:,-1],1.,atol=1e-10):
        raise ValueError("Original source cell mass lost in the ECDF")
    return ids,curves


def _geographic_uncertainty(native:pd.DataFrame,introduced:pd.DataFrame,
                            column:str,grid:np.ndarray,
                            native_curve:np.ndarray,intro_curve:np.ndarray,
                            *,seed_label:str,replicates:int=SPATIAL_BOOTSTRAPS)->dict:
    """Exploratory simultaneous whole-ECDF uncertainty, no biological score.

    Resample the UNION of native and introduced 5-degree cells in one draw:
    locations appearing in both groups get the same geographic multiplicity.
    The 401 pooled-rank display grid stays fixed. Both input distributions
    were made before looking at response differences.
    """
    if replicates<50:
        raise ValueError("Insufficient full-curve geographic resamples")
    n_cells,n_curves=_cell_empirical_curves(native,column,grid)
    i_cells,i_curves=_cell_empirical_curves(introduced,column,grid)
    union=np.union1d(n_cells,i_cells)
    rng=np.random.default_rng(
        int.from_bytes(hashlib.sha256(seed_label.encode()).digest()[:8],"little")
        ^20261021)
    counts=np.stack([
        np.bincount(rng.integers(len(union),size=len(union)),
                    minlength=len(union))
        for _ in range(replicates)
    ]).astype(float)
    n_count=counts[:,np.searchsorted(union,n_cells)]
    i_count=counts[:,np.searchsorted(union,i_cells)]
    # Require geographically dispersed contributors in both groups. Preserve
    # source cells, never replace original photos with pseudo-replicated rows.
    keep=(((n_count>0).sum(axis=1)>=10)&
          ((i_count>0).sum(axis=1)>=10))
    n_count=n_count[keep]
    i_count=i_count[keep]
    if len(n_count)<replicates*MIN_VALID_FRACTION:
        return {"status":"STOP_INSUFFICIENT_INDEPENDENT_GEOGRAPHIC_BOOTSTRAPS",
                "requested":replicates,"valid":int(len(n_count))}
    bootstrap_native=(n_count@n_curves)/n_count.sum(axis=1,keepdims=True)
    bootstrap_intro=(i_count@i_curves)/i_count.sum(axis=1,keepdims=True)
    observed_delta=intro_curve-native_curve
    bootstrap_delta=bootstrap_intro-bootstrap_native
    deviations=bootstrap_delta-observed_delta[None,:]
    # Unstudentised supremum over the WHOLE axis; this is a method for
    # plotting uncertainty, NOT a selected evolutionary effect-size index.
    radius=float(np.quantile(np.max(np.abs(deviations),axis=1),.95))
    simultaneous_low=np.maximum(-1,observed_delta-radius)
    simultaneous_high=np.minimum(1,observed_delta+radius)
    pointwise=np.quantile(bootstrap_delta,[.025,.975],axis=0)
    return {
        "status":"EXPLORATORY_FULL_ECDF_GEOGRAPHIC_CELL_UNCERTAINTY",
        "axis":"same entire source-pooled display grid as native/introduced ECDF",
        "introduced_minus_native_ECDF":observed_delta.tolist(),
        "simultaneous_95pct_spatial_sensitivity_lower":simultaneous_low.tolist(),
        "simultaneous_95pct_spatial_sensitivity_upper":simultaneous_high.tolist(),
        "pointwise_95pct_spatial_sensitivity_lower":pointwise[0].tolist(),
        "pointwise_95pct_spatial_sensitivity_upper":pointwise[1].tolist(),
        "resamples_requested":replicates,
        "resamples_valid":int(len(n_count)),
        "native_5degree_cells":int(len(n_cells)),
        "introduced_5degree_cells":int(len(i_cells)),
        "zero_ECDF_delta_means_equal_cumulative_probability_not_equal_data":True,
        "negative_ECDF_delta_means_introduced_more_rightward_at_that_x":True,
        "simultaneous_confidence_envelope_exploratory_conditional_on_source":True,
        "no_causal_evolution_or_camera_gravity_inference":True,
        "no_biological_tail_or_width_metric_selected":True,
    }


def _curve(native:pd.DataFrame,intro:pd.DataFrame,column:str,*,unit:str,seed_label:str)->dict:
    if min(len(native),len(intro))<100:
        raise ValueError("Insufficient two-range source support")
    left=native[column].to_numpy(float)
    right=intro[column].to_numpy(float)
    if (not np.isfinite(left).all() or not np.isfinite(right).all()
        or (column==ANGLE and
            (np.any((left<0)|(left>180)) or np.any((right<0)|(right>180))))):
        raise ValueError("Invalid source original-pixel rain or image angle")
    wa=_weights(native);wb=_weights(intro)
    grid=_pooled_shared_grid(left,right,wa,wb)
    fn=_ecdf(left,wa,grid)
    fi=_ecdf(right,wb,grid)
    if (np.diff(fn)<-1e-9).any() or (np.diff(fi)<-1e-9).any():
        raise ValueError("Empirical distribution not monotone")
    cell_uncertainty=_geographic_uncertainty(
        native,intro,column,grid,fn,fi,seed_label=seed_label)
    return {
        "whole_curve_spatial_uncertainty":cell_uncertainty,
        "unit":unit,
        "display_only_unselected_common_knots":grid.tolist(),
        "native_empirical_cdf":fn.tolist(),
        "all_introduced_empirical_cdf":fi.tolist(),
        "native_original_observations":int(len(native)),
        "all_introduced_original_observations":int(len(intro)),
        "native_geographic_cells":int(native._cell.nunique()),
        "all_introduced_geographic_cells":int(intro._cell.nunique()),
        "native_original_observers":int(native.user_id.nunique()),
        "all_introduced_original_observers":int(intro.user_id.nunique()),
        "native_and_introduced_geographic_cells_equal_weight":True,
        "no_selected_quantile_spread_or_width_estimand":True,
    }


def explore(wild:pd.DataFrame,photo:pd.DataFrame,
            source_sha:dict|None=None)->dict:
    _validate(wild,photo)
    full=wild.loc[np.isfinite(wild[list(VARIABLES)].to_numpy(float)).all(axis=1)].copy()
    angles=photo.loc[photo[ANGLE].between(0,180) &
         np.isfinite(photo[ANGLE].to_numpy(float))].copy()
    balanced=dated_and_balanced(angles)
    result={}
    for name in SPECIES:
        sw=full.loc[full.taxon_name.eq(name)]
        sp=balanced.loc[balanced.taxon_name.eq(name)]
        result[name]={
            "rainfall_environment_full_distribution":
                _curve(sw.loc[sw.range_state.eq("native")],
                       sw.loc[sw.range_state.eq("introduced")],RAIN,unit="mm_per_year",seed_label=name+"|rain"),
            "image_angle_full_distribution":
                _curve(sp.loc[sp.range_state.eq("native")],
                       sp.loc[sp.range_state.eq("introduced")],ANGLE,unit="image_degrees",seed_label=name+"|angle"),
        }
    return {
        "schema":"S2_exploratory_full_distribution_no_chosen_response_metrics_v1",
        "status":"ORIGINAL_NATIVE_ALL_INTRODUCED_FULL_EMPIRICAL_CURVES",
        "question":"What differences appear in the entire rainfall and image-angle distributions before selecting effects or summary indices?",
        "compared_groups":"Same two species separately, original WCVP-native versus ALL explicit introduced ranges",
        "no_selected_outcome_threshold_or_effect_score":True,
        "complete_ECDF_geographic_uncertainty_without_selected_metric":True,
        "display_grid_definition":"Equal pooled rank spacing per axis (401 requested ranks), plus full minimum and maximum; duplicates removed; no response metric derived from display knots",
        "original_CHELSA_source_sha256":source_sha or {},
        "species":result,
        "privacy":{"observation_id_exported":False,"photographer_id_exported":False,
                   "source_location_exported":False,"original_photo_exported":False},
        "limits":[
            "Original observed distribution and camera-relative image-angle distribution have different measurement-eligibility denominators",
            "No causal rainfall/trait relationship or climate adaptation inferred",
            "This is descriptive ECDF only: no range-width or extreme threshold was prespecified or tested",
            "Image vertical is NOT true gravity, no camera pose correction",
            "ECDF plotting knots are selected for graphical resolution, never treated as independent ecological samples",
        ]
    }


def save(wild:pd.DataFrame,photo:pd.DataFrame,output:Path,
         *,source_sha:dict|None=None)->dict:
    result=explore(wild,photo,source_sha=source_sha)
    output.mkdir(parents=True,exist_ok=True)
    (output/"S2_unselected_full_rain_and_head_angle_distributions.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
        encoding="utf-8")
    return result

"""S2: geographically robust entire conditional capitulum-angle CDF surface.

Question: On the SAME rainfall exposure grid, how do all original photo-angle
distributions compare between the native and ALL explicitly introduced ranges?
Preserve the source-verified, phenotype-blind rainfall display grid and the
complete 0..180 degree axis. Do not select a q90, image-angle threshold, width,
climate band, or an ecological effect after seeing the result.

Native and introduced 5-degree geographic cells are jointly resampled. An
exploratory, unstudentised simultaneous envelope covers the WHOLE common
rainfall-bin × image-angle-CDF surface. Unsupported climate levels are missing
(null), never interpolated or interpreted as no biological difference.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_climate_trait_coupling import original_cell_weights
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced

ANGLE = "orientation_angle_degrees_median"
RAIN = "chelsa_bio12"
SPECIES = ("Cirsium arvense", "Cirsium vulgare")
EXPECTED_PHOTO = 173737
EXPECTED_BALANCED = 94060
DRAWS = 399
MIN_CELLS = 5
MIN_OBSERVATIONS = 30
MIN_COMMON_RAIN_BINS = 5
MIN_VALID_DRAW_FRACTION = .90


def _cell_grid(frame:pd.DataFrame,rain_edges:np.ndarray,
               angle_edges:np.ndarray)->tuple[np.ndarray,np.ndarray]:
    """Per-cell histogram with every cell's total original weight equal to 1."""
    if frame.empty or frame.obs_id.duplicated().any():
        raise ValueError("Missing unmodified source photo observations")
    angles=frame[ANGLE].to_numpy(float)
    rain=frame[RAIN].to_numpy(float)
    if (not np.isfinite(angles).all() or not np.isfinite(rain).all()
        or np.any(angles<0) or np.any(angles>180)):
        raise ValueError("Invalid original image-vertical orientation or climate")
    cell,inv=np.unique(frame._cell.to_numpy(int),return_inverse=True)
    n_by_cell=np.bincount(inv)
    rb=np.clip(np.searchsorted(rain_edges,rain,side="right")-1,0,len(rain_edges)-2)
    ab=np.clip(np.searchsorted(angle_edges,angles,side="right")-1,0,len(angle_edges)-2)
    result=np.zeros((len(cell),len(rain_edges)-1,len(angle_edges)-1),float)
    np.add.at(result,(inv,rb,ab),1./n_by_cell[inv])
    if not np.allclose(result.sum(axis=(1,2)),1.,atol=1e-9):
        raise ValueError("Original equal-cell image weights not conserved")
    return cell,result


def _full_cdf(mass:np.ndarray)->np.ndarray:
    """Native/introduced range-specific P(angle <= x | source rain display bin)."""
    denom=mass.sum(axis=-1,keepdims=True)
    probs=np.divide(mass,denom,out=np.zeros_like(mass),
                    where=denom>0)
    return probs.cumsum(axis=-1)


def _nullable(matrix:np.ndarray,supported:np.ndarray)->list:
    """Entire curve on the registered, outcome-blind raster display grid."""
    return [row.tolist() if bool(ok) else [None]*matrix.shape[-1]
            for row,ok in zip(matrix,supported)]


def _aligned_source_joint(photo:pd.DataFrame,joint:dict)->pd.DataFrame:
    if (joint.get("status")!="ORIGINAL_S2_JOINT_RAIN_AND_IMAGE_ANGLE_ALL_INTRODUCED"
        or not joint.get("all_introduced_remains_primary",False)
        or joint.get("two_species_original_angle_eligible")!=EXPECTED_PHOTO
        or any(joint.get("privacy",{}).values())
        or set(joint.get("species",{}))!=set(SPECIES)):
        raise ValueError("Source-verified prior whole-distribution joint result required")
    if len(photo)!=EXPECTED_PHOTO:
        raise ValueError("Original source photo archive cardinality changed")
    fields={"obs_id","user_id","taxon_name","range_state","_cell",RAIN,ANGLE}
    if not fields<=set(photo) or photo.obs_id.duplicated().any():
        raise ValueError("Original photo identity columns missing")
    keep=photo[ANGLE].between(0,180)&np.isfinite(photo[RAIN].to_numpy(float))
    balanced=dated_and_balanced(photo.loc[keep])
    if len(balanced)!=EXPECTED_BALANCED:
        raise ValueError("Original authenticated balanced photos changed")
    return balanced


def _species_surface(frame:pd.DataFrame,joint_species:dict,
                     *,label:str,draws:int=DRAWS)->dict:
    rain_edges=np.asarray(joint_species["display_climate_BIO12_edges_mm_year"],float)
    angle_edges=np.asarray(joint_species["display_angle_edges_image_degrees"],float)
    if (len(rain_edges)<10 or len(angle_edges)<10
        or np.any(np.diff(rain_edges)<=0) or np.any(np.diff(angle_edges)<=0)):
        raise ValueError("Original phenotype-independent graph axes unavailable")
    part={}
    sample={}
    for state,receipt in (("native",joint_species["native"]),
                          ("introduced",joint_species["all_introduced"])):
        sample[state]=frame.loc[frame.range_state.eq(state)]
        if (len(sample[state])!=receipt["original_photo_observations"]
            or sample[state]._cell.nunique()!=receipt["original_geographic_cells"]):
            raise ValueError("Original same-photo cohort no longer matches joint receipt")
        ids,grid=_cell_grid(sample[state],rain_edges,angle_edges)
        aggregate=grid.sum(axis=0)
        reported=np.asarray(receipt["original_photo_joint_mass"],float)
        if (reported.shape!=aggregate.shape
            or not np.allclose(aggregate/aggregate.sum(),reported,atol=1e-8)):
            raise ValueError("Photo×climate cell grid inconsistent with original joint distribution")
        rain_obs=np.asarray(receipt["photo_n_per_display_rain_bin"],int)
        rain_cells=np.asarray(receipt["photo_cells_per_display_rain_bin"],int)
        active=(grid.sum(axis=2)>0).sum(axis=0)
        if not np.array_equal(active,rain_cells):
            raise ValueError("Original 5degree cell support disagrees with source receipt")
        part[state]=(ids,grid,aggregate,rain_obs,rain_cells)
    supported=np.logical_and(
        np.logical_and(part["native"][3]>=MIN_OBSERVATIONS,
                       part["native"][4]>=MIN_CELLS),
        np.logical_and(part["introduced"][3]>=MIN_OBSERVATIONS,
                       part["introduced"][4]>=MIN_CELLS))
    common=np.flatnonzero(supported)
    if len(common)<MIN_COMMON_RAIN_BINS:
        return {"status":"STOP_INSUFFICIENT_COMMON_RAINFALL_SUPPORT",
                "supported_rain_plot_bins":len(common)}
    observed_native=_full_cdf(part["native"][2])
    observed_intro=_full_cdf(part["introduced"][2])
    diff=observed_intro-observed_native

    if draws<50:
        raise ValueError("Full joint spatial resampling requires >=50 replicates")
    combined=np.union1d(part["native"][0],part["introduced"][0])
    seed=int.from_bytes(hashlib.sha256(label.encode()).digest()[:8],"little")^20261021
    rng=np.random.default_rng(seed)
    multiplier=rng.multinomial(
        len(combined),np.full(len(combined),1/len(combined)),size=draws).astype(float)
    boot={}
    for state in ("native","introduced"):
        cells,grid=part[state][:2]
        cell_mult=multiplier[:,np.searchsorted(combined,cells)]
        # All photo angles and rainfall intervals, no chosen trait tail.
        replicate_hists=np.einsum("kc,cij->kij",cell_mult,grid,optimize=True)
        boot[state]=replicate_hists
    n_mass=boot["native"].sum(axis=-1)
    i_mass=boot["introduced"].sum(axis=-1)
    valid=np.logical_and(
        (n_mass[:,supported]>0).all(axis=1),
        (i_mass[:,supported]>0).all(axis=1))
    if valid.sum()<draws*MIN_VALID_DRAW_FRACTION:
        return {"status":"STOP_INSUFFICIENT_CELL_BOOTSTRAP_SUPPORT",
                "requested_draws":draws,
                "valid_draws":int(valid.sum()),
                "common_original_rainfall_bins":int(len(common))}
    c_native=_full_cdf(boot["native"][valid])
    c_intro=_full_cdf(boot["introduced"][valid])
    replicates=c_intro-c_native
    # Treat every supported rainfall plot bin and angle cut jointly.
    deviations=replicates[:,supported,:]-diff[None,supported,:]
    largest=np.max(np.abs(deviations),axis=(1,2))
    radius=float(np.quantile(largest,.95))
    upper=np.minimum(1.,diff+radius)
    lower=np.maximum(-1.,diff-radius)
    positive=np.logical_and(supported[:,None],lower>0)
    negative=np.logical_and(supported[:,None],upper<0)
    pointwise=np.quantile(replicates,[.025,.975],axis=0)
    return {
        "status":"EXPLORATORY_SIMULTANEOUS_SAME_RAIN_FULL_HEAD_CDF_CELL_ENVELOPE",
        "rainfall_grid_from_original_wild_not_angle_selected":True,
        "full_image_angle_axis_not_preselected_threshold":True,
        "native_vs_ALL_introduced_only":True,
        "source_climate_display_edges_mm_year":rain_edges.tolist(),
        "image_angle_CDF_upper_edges_degrees":angle_edges[1:].tolist(),
        "rainfall_display_bins_originally_supported_both_ranges":
            supported.tolist(),
        "common_rainfall_bins_originally_supported":int(len(common)),
        "observed_conditional_head_angle_CDF_introduced_minus_native":
            _nullable(diff,supported),
        "simultaneous_95pct_spatial_sensitivity_lower":
            _nullable(lower,supported),
        "simultaneous_95pct_spatial_sensitivity_upper":
            _nullable(upper,supported),
        "pointwise_95pct_sensitivity_lower":
            _nullable(pointwise[0],supported),
        "pointwise_95pct_sensitivity_upper":
            _nullable(pointwise[1],supported),
        "native_original_photo_n_by_rain_display_bin":part["native"][3].tolist(),
        "introduced_original_photo_n_by_rain_display_bin":
            part["introduced"][3].tolist(),
        "native_original_cells_by_rain_display_bin":part["native"][4].tolist(),
        "introduced_original_cells_by_rain_display_bin":
            part["introduced"][4].tolist(),
        "simultaneous_interval_diagnostic":{
            "region":"entire original supported rain × image CDF display surface",
            "absolute_supremum_deviation_95percent":radius,
            "all_supported_surface_points_with_lower_above_zero":int(positive.sum()),
            "all_supported_surface_points_with_upper_below_zero":int(negative.sum()),
            "raw_cell_draws_requested":draws,
            "draws_valid_for_all_common_rain_bins":int(valid.sum()),
            "spatial_resampling_unit":"union original 5degree geographic cells",
            "conditional_not_camera_adjusted_or_genotype_replication":True,
            "not_a_reject_null_test_for_each_rain_bin":True,
        },
    }


def evaluate(photo:pd.DataFrame,joint:dict,*,draws:int=DRAWS)->dict:
    frame=_aligned_source_joint(photo,joint)
    species={}
    for name in SPECIES:
        s=frame.loc[frame.taxon_name.eq(name)]
        species[name]=_species_surface(s,joint["species"][name],
                                       label=name+"|S2_conditional_curve",draws=draws)
    return {
        "schema":"azami_S2_entire_conditional_image_orientation_CDF_spatial_20261022_v1",
        "status":"SOURCE_VERIFIED_SAME_RAIN_UNSELECTED_ANGLE_FULL_CDF_CELL_SENSITIVITY",
        "n_source_photo_angle_eligible":EXPECTED_PHOTO,
        "n_photographer_cell_dated_balanced":len(frame),
        "source_original_five_CHELSA_SHA256":joint.get("source_raster_SHA256",{}),
        "species":species,
        "range_definitions":"same exact WCVP native and ALL introduced as original S2",
        "no_chosen_width_quantile_extreme_tail_effect_metric":True,
        "camera_orientation_is_not_true_head_gravity":True,
        "source_privacy":{"obs_ids":False,"user_ids":False,"coordinates":False,"photos":False},
        "limitations":[
            "Conditional curves use photo-eligible sample, not all wild source occurrences",
            "5degree geographic cells are not independently introduced genotypes",
            "Same-rain comparison does not control for temperature, season, photo framing or history",
            "Simultaneous band is exploratory and conditional on 5degree-cell sampling",
            "BIO12 annual climatic normal is not the rain during flowering"
        ],
    }


def save(photo:pd.DataFrame,joint:dict,output:Path,*,draws:int=DRAWS)->dict:
    out=evaluate(photo,joint,draws=draws)
    output.mkdir(parents=True,exist_ok=True)
    path=output/"S2_full_conditional_rain_angle_ECDF_cell_sensitivity.json"
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                    encoding="utf-8")
    return out

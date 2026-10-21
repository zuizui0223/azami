"""S2: source-defined novel wet-edge -> conditional *image-angle distribution*.

Climate envelope and nearest-neighbour links are made without reading image
ANGLE VALUES. The full WCVP source (not trait-selected) supplies the native
rainfall reference; introduced comparisons stay within species and continent.
Original observation IDs/user IDs/coordinates never leave the private runner.

Exploratory image relative outcomes only, not a gravity-referenced posture,
independent plant phenotype, plasticity, or adaptive evolution.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wasserstein_distance

from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.native_climate_trait_coupling import (
    original_cell_weights, weighted_quantiles
)
from analysis.v4.native_introduced_climate_transfer import dated_and_balanced
from analysis.v4.native_intro_continent_orientation import CONTINENTS

SPECIES = ("Cirsium arvense", "Cirsium vulgare")
RAIN = "chelsa_bio12"
OTHER = tuple(x for x in VARIABLES if x != RAIN)
COL_ANGLE = "orientation_angle_degrees_median"
PREDECLARED_NATIVE_QUANTILES = (0.50, 0.75, 0.90, 0.95, 0.975)
REF_LOWER = 0.75
REF_UPPER = 0.95
MIN_ROWS = 30
MIN_INDEPENDENT_CELLS = 15
MIN_ORIGINAL_OBSERVERS = 10
MIN_MATCHED_ESS = 30
MAX_STANDARDISED_DISTANCE = 1.25
MAX_FEATURE_DISTANCE = 2.0
MAX_STANDARDISED_MEAN_IMBALANCE = 0.25
NEIGHBOURS = 3
SEED = 20261021


def _require_source(frame: pd.DataFrame, *, need_angle: bool) -> None:
    required = {"obs_id", "taxon_name", "range_state", "user_id",
                "latitude", "longitude", "tdwg1", "_cell",
                "observed_on", "quality_grade", *VARIABLES}
    if need_angle:
        required.add(COL_ANGLE)
    if not required <= set(frame):
        raise ValueError("Source identity, climate or image qualification fields missing")
    if frame.obs_id.duplicated().any() or frame.user_id.isna().any():
        raise ValueError("Repeated observation or missing original photographer")
    if not set(frame.taxon_name) <= set(SPECIES):
        raise ValueError("Unexpected focal species")
    if not set(frame.range_state) <= {"native", "introduced"}:
        raise ValueError("Unverified WCVP establishment state")


def _finite(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.loc[np.isfinite(frame[list(VARIABLES)].to_numpy(float)).all(axis=1)].copy()


def native_envelope(source_native: pd.DataFrame) -> dict:
    """Angle-blind, balanced source wildlife only; never the measured subset."""
    if len(source_native) < 200 or source_native._cell.nunique() < 20:
        raise ValueError("Not enough native ORIGINAL wild climate support")
    weights = original_cell_weights(source_native)
    q = weighted_quantiles(
        source_native[RAIN].to_numpy(float), weights,
        quantiles=PREDECLARED_NATIVE_QUANTILES
    )
    if any(a >= b for a,b in zip(q[:-1],q[1:])):
        raise ValueError("Degenerate predeclared native rainfall quantiles")
    bounds={}
    for column in OTHER:
        low,high=weighted_quantiles(
            source_native[column].to_numpy(float), weights)
        if low >= high:
            raise ValueError("Native nonrain climatic envelope degenerate")
        bounds[column]={"q05":low,"q95":high}
    z_reference = np.column_stack([
        source_native[column].to_numpy(float) for column in OTHER
    ])
    mean=z_reference.mean(axis=0)
    sd=z_reference.std(axis=0)
    if not np.isfinite(sd).all() or (sd<=1e-10).any():
        raise ValueError("Native source scale missing for nonrain matching")
    return {
        "BIO12_native_source_q50_q75_q90_q95_q975_mm":list(q),
        "native_nonrain_marginal_bounds":bounds,
        "nonrain_source_means":dict(zip(OTHER,mean.tolist())),
        "nonrain_source_sds":dict(zip(OTHER,sd.tolist())),
        "source_n":int(len(source_native)),
        "source_cells":int(source_native._cell.nunique()),
        "definition":"Angle-blind cell-equal native source, not trait-selected; q95 wet-edge, q75..q95 wet analogue",
    }


def _source_band(frame: pd.DataFrame, q75: float, q95: float) -> np.ndarray:
    rain=frame[RAIN].to_numpy(float)
    return np.select(
        [rain>q95, (rain>=q75)&(rain<=q95)],
        ["beyond_native_wet_edge", "native_known_wet_band"],
        default="drier_than_native_q75"
    ).astype(str)


def _source_support(frame: pd.DataFrame, photo_ids: set[str],
                    q75: float, q95: float) -> dict:
    """Selection denominators: never infer a photographed plant from absence."""
    out={}
    for region,unit in frame.groupby("range_state",sort=True):
        def counts(f):
            return {
                "source_wild_observations":int(len(f)),
                "original_observer_cell_units":int(f[["_cell","user_id"]].drop_duplicates().shape[0]),
                "original_observers":int(f.user_id.nunique()),
                "coarse_cells":int(f._cell.nunique()),
                "original_obs_with_eligible_angle":int(f.obs_id.isin(photo_ids).sum()),
                "quality_grade_missing_source":int(f.quality_grade.isna().sum()),
            }
        out[region]={"all":counts(unit)}
        categories=_source_band(unit,q75,q95)
        for category in ("beyond_native_wet_edge","native_known_wet_band",
                         "drier_than_native_q75"):
            out[region][category]=counts(unit.loc[categories==category])
    return out


def _features(frame: pd.DataFrame, envelope: dict) -> np.ndarray:
    """Only OTHER climate, timing; deliberately excludes BIO12 and angle."""
    columns=[]
    for name in OTHER:
        mean=envelope["nonrain_source_means"][name]
        sd=envelope["nonrain_source_sds"][name]
        columns.append((frame[name].to_numpy(float)-mean)/sd)
    m=frame["_month"].to_numpy(float)
    year=frame["_year"].to_numpy(float)
    columns += [np.sin(2*np.pi*m/12)/0.75,
                np.cos(2*np.pi*m/12)/0.75, (year-2015)/8.0]
    values=np.column_stack(columns)
    if not np.isfinite(values).all():
        raise ValueError("Nonrain/date features invalid")
    return values


def _counts(frame: pd.DataFrame) -> dict:
    return {"n":int(len(frame)), "cells":int(frame._cell.nunique()),
            "observers":int(frame.user_id.nunique())}


def _cell_weights(frame: pd.DataFrame) -> np.ndarray:
    return original_cell_weights(frame)


def _ess(w: np.ndarray) -> float:
    v=np.asarray(w,float)
    if v.size==0 or v.sum()<=0:
        return 0.
    return float((v.sum()**2)/np.square(v).sum())


def match_without_angles(wet: pd.DataFrame, reference: pd.DataFrame,
                         envelope: dict) -> dict:
    """Three nearest comparators per case within quality stratum, no rainfall or angles.

    Observational nearest-neighbour overlap, NOT 1:1 pairing of plants or a
    causal propensity score. Caliper is RMS on 7 predeclared standardised
    nonrain/season/year features, with per-coordinate upper bound.
    """
    xw=_features(wet,envelope)
    xr=_features(reference,envelope)
    case_rows=[];control_rows=[];distances=[]
    grade_w=wet.quality_grade.fillna("unknown").astype(str).to_numpy()
    grade_r=reference.quality_grade.fillna("unknown").astype(str).to_numpy()
    for grade in np.unique(grade_w):
        widx=np.flatnonzero(grade_w==grade)
        ridx=np.flatnonzero(grade_r==grade)
        if len(ridx)<NEIGHBOURS or not len(widx):
            continue
        tree=cKDTree(xr[ridx])
        d, near=tree.query(xw[widx],k=NEIGHBOURS)
        for row,indices in enumerate(np.atleast_2d(near)):
            for j in indices:
                ctrl=int(ridx[int(j)])
                resid=np.abs(xw[widx[row]]-xr[ctrl])
                rms=float(np.sqrt(np.mean(resid**2)))
                if rms<=MAX_STANDARDISED_DISTANCE and float(resid.max())<=MAX_FEATURE_DISTANCE:
                    case_rows.append(int(widx[row]))
                    control_rows.append(ctrl)
                    distances.append(rms)
    if not case_rows:
        return {"status":"NO_NONRAIN_SEASONAL_ANALOGUE","wet_indices":np.array([],int),
                "control_indices":np.array([],int)}
    case_rows=np.asarray(case_rows,int)
    control_rows=np.asarray(control_rows,int)
    counts=np.bincount(case_rows,minlength=len(wet))
    case_unique=np.flatnonzero(counts>0)
    # Each admitted wet source observation has unit weight within its 5-degree
    # cell. Share that weight across its own matched controls.
    chosen=wet.iloc[case_unique]
    case_base=_cell_weights(chosen)
    weights_lookup=np.zeros(len(wet),float)
    weights_lookup[case_unique]=case_base
    edge_weights=weights_lookup[case_rows]/counts[case_rows]
    ctrl_weights=np.bincount(control_rows,weights=edge_weights,
                             minlength=len(reference))
    ctrl_unique=np.flatnonzero(ctrl_weights>0)
    wt_case=case_base/case_base.sum()
    wt_ctrl=ctrl_weights[ctrl_unique]/ctrl_weights[ctrl_unique].sum()
    feature_difference=(wt_case@xw[case_unique])-(wt_ctrl@xr[ctrl_unique])
    return {
        "status":"ANGLE_BLIND_MATCH_COMPLETE",
        "wet_indices":case_unique,"control_indices":ctrl_unique,
        "wet_weights":wt_case,"control_weights":wt_ctrl,
        "matched_pairs":int(len(case_rows)),
        "maximum_absolute_feature_mean_difference_native_scaled":
            float(np.abs(feature_difference).max()),
        "feature_mean_differences_native_scaled":
            dict(zip([*OTHER,"season_sine","season_cosine","observation_year"],
                     feature_difference.tolist())),
        "max_matched_nonrain_season_year_RMS":
            float(np.max(distances)),
        "matched_wet_original_cells":int(chosen._cell.nunique()),
        "matched_control_original_cells":int(reference.iloc[ctrl_unique]._cell.nunique()),
        "unique_controls":int(len(ctrl_unique)),
        "wet_ESS":_ess(wt_case),"control_ESS":_ess(wt_ctrl),
    }


def geographic_block_uncertainty(wet:pd.DataFrame, ref:pd.DataFrame,
                                 matched:dict, *, replicates:int=400)->dict:
    """A sensitivity interval from joint 5-degree cell resampling.

    Keep outcome-blind matching links fixed; resample the UNION of geographic
    cells shared by wet and reference. This captures geographic leverage in
    the observed photo sample but NOT matching uncertainty, camera error,
    separate introduction events or any causal estimand.
    """
    w=wet.iloc[matched["wet_indices"]]
    r=ref.iloc[matched["control_indices"]]
    yw=w[COL_ANGLE].to_numpy(float)
    yr=r[COL_ANGLE].to_numpy(float)
    ww=matched["wet_weights"]
    wr=matched["control_weights"]
    cells=np.unique(np.concatenate([
        w._cell.to_numpy(int),r._cell.to_numpy(int)
    ]))
    wcell=np.searchsorted(cells,w._cell.to_numpy(int))
    rcell=np.searchsorted(cells,r._cell.to_numpy(int))
    rng=np.random.default_rng(SEED)
    values=[]
    for _ in range(replicates):
        picks=rng.integers(len(cells),size=len(cells))
        counts=np.bincount(picks,minlength=len(cells))
        aw=ww*counts[wcell]
        ar=wr*counts[rcell]
        ok1=aw>0
        ok2=ar>0
        if ok1.sum()<MIN_ROWS or ok2.sum()<MIN_ROWS:
            continue
        if (len(np.unique(wcell[ok1]))<8
            or len(np.unique(rcell[ok2]))<8):
            continue
        aw=aw[ok1]/aw[ok1].sum()
        ar=ar[ok2]/ar[ok2].sum()
        qw=weighted_quantiles(yw[ok1],aw,(.10,.50,.90))
        qr=weighted_quantiles(yr[ok2],ar,(.10,.50,.90))
        values.append([*(float(a-b) for a,b in zip(qw,qr)),
                       float(wasserstein_distance(yw[ok1],yr[ok2],
                                                  u_weights=aw,v_weights=ar))])
    if len(values)<int(replicates*.90):
        return {
            "status":"NOT_ESTIMABLE_TOO_FEW_VALID_GEOGRAPHIC_BLOCK_DRAWS",
            "valid_draws":len(values),"requested_draws":replicates
        }
    draws=np.asarray(values,float)
    intervals=np.quantile(draws,[.025,.975],axis=0)
    return {
        "status":"GEOGRAPHIC_CELL_SENSITIVITY_INTERVAL_NOT_CAUSAL_OR_CAMERA_VALIDATED",
        "method":"Joint 5-degree source-cell bootstrap; matched neighbour links held fixed",
        "geographic_blocks":int(len(cells)),
        "valid_draws":len(values),
        "Q10_Q50_Q90_delta_95pct_intervals":[
            [float(intervals[0,i]),float(intervals[1,i])] for i in range(3)],
        "W1_95pct_sensitivity_interval":[float(intervals[0,3]),float(intervals[1,3])],
        "Q90_positive_fraction_of_draws":float((draws[:,2]>0).mean()),
        "not_adjusted_for_camera_and_measurement_error":True,
        "matching_and_selection_variance_not_estimated":True,
        "not_an_independent_population_confidence_interval":True,
    }


def conditional_distribution(wet:pd.DataFrame, ref:pd.DataFrame,
                             matched:dict)->dict:
    """Point estimates; independent camera calibration and cluster CI pending."""
    if matched["status"]!="ANGLE_BLIND_MATCH_COMPLETE":
        return {"status":matched["status"]}
    iw=matched["wet_indices"];ir=matched["control_indices"]
    a=wet.iloc[iw];b=ref.iloc[ir]
    limits={
        "wet":_counts(a), "reference":_counts(b),
        "wet_effective_n":matched["wet_ESS"],
        "reference_effective_n":matched["control_ESS"],
        "nonrain_covariate_max_abs_scaled_imbalance":
            matched["maximum_absolute_feature_mean_difference_native_scaled"],
    }
    valid=(
        min(len(a),len(b))>=MIN_ROWS
        and min(a._cell.nunique(),b._cell.nunique())>=MIN_INDEPENDENT_CELLS
        and min(a.user_id.nunique(),b.user_id.nunique())>=MIN_ORIGINAL_OBSERVERS
        and min(matched["wet_ESS"],matched["control_ESS"])>=MIN_MATCHED_ESS
        and matched["maximum_absolute_feature_mean_difference_native_scaled"]
              <=MAX_STANDARDISED_MEAN_IMBALANCE
    )
    if not valid:
        return {"status":"INSUFFICIENT_MATCHED_SUPPORT_OR_COVARIATE_BALANCE",
                "support":limits,
                "no_biological_zero_inferred":True}
    yw=a[COL_ANGLE].to_numpy(float)
    yr=b[COL_ANGLE].to_numpy(float)
    ww=matched["wet_weights"]
    wr=matched["control_weights"]
    if (not np.isfinite(yw).all() or not np.isfinite(yr).all()
        or ((yw<0)|(yw>180)).any() or ((yr<0)|(yr>180)).any()):
        raise ValueError("Invalid original image angle after source-matching")
    quantiles=(.10,.50,.90)
    qa=weighted_quantiles(yw,ww,quantiles)
    qb=weighted_quantiles(yr,wr,quantiles)
    result = {
        "status":"EXPLORATORY_MATCHED_DISTRIBUTION_ESTIMATED_NOT_INDEPENDENT_CONFIRMED",
        "support":limits,
        "comparison":"Same-species and same introduced continent; novel wet-edge versus native q75-q95 rain analogue; nearest nonrain/season/year climate and equal quality grade",
        "response":"EXIF-corrected image vertical degrees; 0 up, 90 horizontal, 180 down",
        "wet_Q10_Q50_Q90":qa,
        "reference_Q10_Q50_Q90":qb,
        "wet_minus_reference_Q10_Q50_Q90":[float(a-b) for a,b in zip(qa,qb)],
        "wet_minus_reference_interquantile_width":
            float((qa[2]-qa[0])-(qb[2]-qb[0])),
        "weighted_one_dimensional_wasserstein_W1_degrees":
            float(wasserstein_distance(yw,yr,u_weights=ww,v_weights=wr)),
        "wet_image_horizontal_or_down_fraction":
            float(ww[yw>=90].sum()),
        "reference_image_horizontal_or_down_fraction":
            float(wr[yr>=90].sum()),
        "geographic_block_uncertainty":geographic_block_uncertainty(wet,ref,matched),
        "no_gravity_calibration_or_causal_interpretation":True,
    }
    # Geographic uncertainty stays separate from raw point estimates and
    # treats the frozen climate-only matches as fixed; camera and matching
    # model uncertainty are not incorporated.
    from analysis.v4.native_wet_edge_uncertainty import geographic_block_bootstrap
    result["geographic_cell_uncertainty"] = geographic_block_bootstrap(
        a,b,ww,wr,
        label=str(a.taxon_name.iloc[0])+"|"+str(a.tdwg1.iloc[0])+"|BIO12-native-q95",
    )
    return result


def evaluate(source_wild:pd.DataFrame, photo:pd.DataFrame,
             *,source_climate_sha:dict|None=None) -> dict:
    """Run A support ledger then B *only* where angle-blind gates pass."""
    _require_source(source_wild,need_angle=False)
    _require_source(photo,need_angle=True)
    if len(source_wild)<1000 or len(photo)<1000:
        raise ValueError("Not full pre-phototype two-species source")
    src=_finite(source_wild)
    measured=_finite(photo)
    measured=measured.loc[
        np.isfinite(measured[COL_ANGLE].to_numpy(float))
        & measured[COL_ANGLE].between(0,180)
    ].copy()
    measured=dated_and_balanced(measured)
    measured["_source_grade"]=measured.quality_grade.fillna("unknown").astype(str)
    if len(measured)<100:
        raise ValueError("Original photo response missing after date/QC gate")
    output={}
    photo_ids=set(photo.obs_id.astype(str))
    for species in SPECIES:
        n=src.loc[src.taxon_name.eq(species)&src.range_state.eq("native")]
        en=native_envelope(n)
        q50,q75,q90,q95,q975=en["BIO12_native_source_q50_q75_q90_q95_q975_mm"]
        sp=src.loc[src.taxon_name.eq(species)]
        subj={
            "native_original_unselected_climate_reference":{
                "BIO12_q50_q75_q90_q95_q975_mm":[q50,q75,q90,q95,q975],
                "n_source":en["source_n"],"n_cells":en["source_cells"]
            },
            "angle_blind_source_coverage":_source_support(sp,photo_ids,q75,q95),
            "introduced_by_continent":{},
        }
        ps=measured.loc[measured.taxon_name.eq(species) &
                        measured.range_state.eq("introduced")].copy()
        for tdwg,reg in ps.groupby("tdwg1",sort=True):
            if pd.isna(tdwg) or int(tdwg) not in CONTINENTS:
                raise ValueError("Original TDWG continent not authoritative")
            name=CONTINENTS[int(tdwg)]
            subset=sp.loc[sp.range_state.eq("introduced") &
                          sp.tdwg1.eq(tdwg)]
            w_src=subset.loc[subset[RAIN]>q95]
            ref_src=subset.loc[subset[RAIN].between(q75,q95)]
            wetlands=reg.loc[reg[RAIN]>q95]
            comparators=reg.loc[reg[RAIN].between(q75,q95)]
            entry={
                "source_wet_edge":{"n":len(w_src),"cells":int(w_src._cell.nunique())},
                "source_native_wet_band":{"n":len(ref_src),"cells":int(ref_src._cell.nunique())},
                "photographed_wet_edge":_counts(wetlands),
                "photographed_reference":_counts(comparators),
                "angle_prevalence_not_plant_prevalence":True,
                "comparison":"BIO12>native q95 versus native q75..q95 within one introduced continent; other four exposures/time matched, not BIO12"
            }
            # Gate A: do not use angles for matching or admission decisions.
            if (len(wetlands)>=MIN_ROWS and len(comparators)>=MIN_ROWS
                and wetlands._cell.nunique()>=MIN_INDEPENDENT_CELLS
                and comparators._cell.nunique()>=MIN_INDEPENDENT_CELLS
                and wetlands.user_id.nunique()>=MIN_ORIGINAL_OBSERVERS
                and comparators.user_id.nunique()>=MIN_ORIGINAL_OBSERVERS):
                mat=match_without_angles(wetlands,comparators,en)
                entry["matching"]={key:value for key,value in mat.items()
                                   if not isinstance(value,np.ndarray)}
                entry["conditional_angle_distribution"]=conditional_distribution(
                    wetlands,comparators,mat)
            else:
                entry["matching"]={"status":"STOP_INSUFFICIENT_INDEPENDENT_WET_EDGE_OR_ANALOGUE_CELLS"}
                entry["conditional_angle_distribution"]={
                    "status":"NOT_IDENTIFIED_NOT_ZERO"}
            # Further native-Z novelty is a *separate* descriptive gate;
            # a within-introduction analogue does not prove source-climate identity.
            within_native_nonrain=np.ones(len(wetlands),dtype=bool)
            for key,b in en["native_nonrain_marginal_bounds"].items():
                v=wetlands[key].to_numpy(float)
                within_native_nonrain&=(v>=b["q05"])&(v<=b["q95"])
            entry["wet_edge_within_native_nonrain_marginal_envelope"]={
                "n":int(within_native_nonrain.sum()),
                "n_outside":int((~within_native_nonrain).sum()),
                "not_a_joint_five_dimensional_niche_test":True,
            }
            subj["introduced_by_continent"][name]=entry
        output[species]=subj
    return {
        "schema_version":"S2_wet_edge_angle_blind_source_and_matched_distribution_20261021_v1",
        "status":"SOURCE_VERIFIED_WET_EDGE_PREFLIGHT_AND_CONDITIONAL_EXPLORATORY_WHERE_ADMITTED",
        "source_two_species_full_original_wild_observations":297622,
        "source_orientation_eligible_before_date_and_cell_balancing":173737,
        "source_wild_cell_redistributed_n":int(len(source_wild)),
        "source_wild_complete_climate_n":int(len(src)),
        "orientation_original_image_balanced_climate_complete_n":int(len(measured)),
        "registered_threshold":"source (non-photo-selected) native BIO12 cell-weighted q95; q75..q95 comparator",
        "calipers":{"nonrain_season_year_rms":MAX_STANDARDISED_DISTANCE,
                    "max_individual_scaled_axis":MAX_FEATURE_DISTANCE,
                    "max_post_match_scaled_mean_diff":MAX_STANDARDISED_MEAN_IMBALANCE,
                    "k":NEIGHBOURS},
        "focal":output,
        "source_five_original_CHELSA_SHA256":source_climate_sha or {},
        "privacy":{"any_individual_source_identifier_exported":False,
                   "any_original_geo_coordinate_exported":False,
                   "any_original_photo_exported":False},
        "limits":[
            "WCVP/TDWG native/introduced is not confirmed genetic provenance",
            "Photo vertical is not gravity, camera attitude unvalidated",
            "Image assessability varies with geography, source effort, seasons and head development",
            "KNN comparison is observational; no balance or source cells => no result",
            "Pointwise matched percentiles and W1 have NO cell/photographer-block confidence intervals yet",
            "Niche q95 measures sampled source occupancy, not fundamental niche limit",
            "Within introduced wet-edge is not the same individual before/after colonization",
        ],
    }


def save(source_wild:pd.DataFrame,photo:pd.DataFrame,output:Path,
         source_climate_sha:dict|None=None)->dict:
    result=evaluate(source_wild,photo,source_climate_sha=source_climate_sha)
    output.mkdir(parents=True,exist_ok=True)
    target=output/"S2_wet_edge_angle_blind_distribution.json"
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                 allow_nan=False)+"\n",encoding="utf-8")
    return result

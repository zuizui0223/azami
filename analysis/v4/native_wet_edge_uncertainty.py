"""Cluster-conditional uncertainty for S2 matched image-angle distributions.

Resample geographical 5-degree cells from the UNION of already matched wet
and control records, keeping the same cell draw for both arms. Matches were
selected without angles. This is a *fixed-match, geographical-sampling*
sensitivity, not inference over camera calibration, plant genets, or matching
uncertainty. No geographic IDs or original identifiers are returned.
"""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance
from analysis.v4.native_climate_trait_coupling import weighted_quantiles

DEFAULT_REPLICATES = 399
MIN_VALID_FRACTION = .90
QUANTILES = (.10,.50,.90)
ANGLE_COLUMN = "orientation_angle_degrees_median"

def _weighted_summary(yw: np.ndarray, yr: np.ndarray,
                      ww: np.ndarray, wr: np.ndarray) -> dict:
    left=ww>0
    right=wr>0
    if not left.any() or not right.any():
        raise ValueError("Empty geographic bootstrap arm")
    y1=yw[left]
    y0=yr[right]
    w1=ww[left]/ww[left].sum()
    w0=wr[right]/wr[right].sum()
    q1=weighted_quantiles(y1,w1,QUANTILES)
    q0=weighted_quantiles(y0,w0,QUANTILES)
    qdiff=np.asarray(q1)-np.asarray(q0)
    return {
        "qdiff":qdiff,
        "qwidth":float(qdiff[2]-qdiff[0]),
        "w1":float(wasserstein_distance(y1,y0,u_weights=w1,v_weights=w0)),
        "horizontal_or_down_fraction_delta":
            float(w1[y1>=90].sum()-w0[y0>=90].sum()),
    }


def geographic_block_bootstrap(
    wet:pd.DataFrame,
    reference:pd.DataFrame,
    wet_weights:np.ndarray,
    reference_weights:np.ndarray,
    *,
    label:str,
    replicates:int=DEFAULT_REPLICATES
)->dict:
    """Shared-cell geographic block resampling of two fixed matched samples.

    Resampling a 5-degree cell jointly in both arms acknowledges same-cell
    observations are spatially linked. It does not restore uncertainty due
    to matching controls or original photographer/camera measurement error.
    """
    if replicates<30:
        raise ValueError("At least 30 fixed-source spatial replicates required")
    for sample,w in ((wet,wet_weights),(reference,reference_weights)):
        if "_cell" not in sample or ANGLE_COLUMN not in sample:
            raise ValueError("Geographic block or original image angle missing")
        y=sample[ANGLE_COLUMN].to_numpy(float)
        if (len(sample)!=len(w) or len(sample)==0 or
            not np.isfinite(y).all() or
            not np.isfinite(w).all() or np.any(w<=0) or
            np.any((y<0)|(y>180))):
            raise ValueError("Bad fixed-match bootstrap input")
    wet_cells=wet["_cell"].to_numpy()
    ref_cells=reference["_cell"].to_numpy()
    union,wi,ri=np.unique(
        np.concatenate((wet_cells,ref_cells)),return_inverse=True
    )
    na=len(wet)
    ia=wi[:na]
    ib=wi[na:]
    if np.unique(wet_cells).size<15 or np.unique(ref_cells).size<15:
        return {"status":"STOP_INSUFFICIENT_INDEPENDENT_SPATIAL_BLOCKS"}
    yw=wet[ANGLE_COLUMN].to_numpy(float)
    yr=reference[ANGLE_COLUMN].to_numpy(float)
    wa=np.asarray(wet_weights,float)
    wb=np.asarray(reference_weights,float)
    wa/=wa.sum()
    wb/=wb.sum()
    seed=int.from_bytes(hashlib.sha256(label.encode("utf-8")).digest()[:8],
                        "little") ^ 20261021
    rng=np.random.default_rng(seed)
    qd=[]
    qw=[]
    wdist=[]
    frac=[]
    for _ in range(replicates):
        # All recorded cells, not individual photographs, are resampled.
        n=np.bincount(rng.integers(len(union),size=len(union)),minlength=len(union))
        a=wa*n[ia]
        b=wb*n[ib]
        if a.sum()<=0 or b.sum()<=0:
            continue
        s=_weighted_summary(yw,yr,a,b)
        qd.append(s["qdiff"])
        qw.append(s["qwidth"])
        wdist.append(s["w1"])
        frac.append(s["horizontal_or_down_fraction_delta"])
    if len(qd)<MIN_VALID_FRACTION*replicates:
        return {"status":"STOP_TOO_FEW_VALID_GEO_BLOCK_DRAWS",
                "valid":len(qd),"requested":replicates}
    def interval(a):
        return [float(x) for x in np.quantile(np.asarray(a,float),
                                               (.025,.50,.975))]
    delta=np.asarray(qd,float)
    return {
        "status":"EXPLORATORY_FIXED_MATCH_GEOGRAPHIC_BLOCK_INTERVALS",
        "replicates_requested":replicates,
        "replicates_valid":len(qd),
        "unique_union_5degree_cells":int(len(union)),
        "wet_5degree_cells":int(np.unique(wet_cells).size),
        "control_5degree_cells":int(np.unique(ref_cells).size),
        "delta_Q10_Q50_Q90_2p5_50_97p5":[interval(delta[:,j]) for j in range(3)],
        "delta_interquantile_width_2p5_50_97p5":interval(qw),
        "W1_distance_2p5_50_97p5":interval(wdist),
        "delta_horizontal_or_down_fraction_2p5_50_97p5":interval(frac),
        "fixed_angle_blind_match":True,
        "geographic_cell_cluster_resampled_jointly_across_arms":True,
        "no_match_selection_or_camera_error_bootstrap":True,
        "confidence_intervals_unadjusted_exploratory":True,
        "no_original_coordinates_or_cell_identifiers_exported":True,
    }

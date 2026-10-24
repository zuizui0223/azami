"""Limited original V4 external-involucre phenotypic LEGO bridge to aza3.

Retains only image-derived bract PROJECTION and spread as candidate proxies,
not authentic anatomical free phyllary posture, spine length, stickiness, or
functional access. Different common cohorts must never be pooled with the 8
core-construct 28-edge V4 analysis or labelled V3's 33/36.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.v4.native_integration import association
from analysis.v4.s1_s2_lego_triple import maxent_from_observed

COLUMNS=("orientation_image_vertical_angle","corolla_lab_chroma",
         "bract_projection_p95","bract_spread_fraction")
FEATURES={
 "presentation":[COLUMNS[0]],
 "colour_chroma":[COLUMNS[1]],
 "candidate_phyllary_envelope":[COLUMNS[2],COLUMNS[3]],
}
TRIAD=(COLUMNS[0],COLUMNS[1],COLUMNS[2])
MIN_TAXA=20
MIN_N_PER_TAXON=5


def valid_source(native, env):
    missing=[c for c in COLUMNS if c not in native or c+"__eligible" not in native]
    if missing:
        return None,{"status":"SOURCE_ENDPOINTS_NOT_MEASURED","missing_endpoint_or_QC":missing}
    keep=np.ones(len(native),bool)
    for c in COLUMNS:
        keep &= np.isfinite(pd.to_numeric(native[c],errors="coerce").to_numpy(float))
        keep &= native[c+"__eligible"].fillna(False).astype(bool).to_numpy()
    original= native.loc[keep].copy()
    measured_n=len(original)
    measured_taxa=original.accepted_key.nunique()
    keep_env=np.isfinite(original[list(env)].to_numpy(float)).all(axis=1)
    complete=original.loc[keep_env].copy()
    counts=complete.groupby("accepted_key").size()
    complete=complete.loc[complete.accepted_key.isin(
        counts[counts>=MIN_N_PER_TAXON].index)].copy()
    receipt={"total_bract_angle_colour_same_observation_QC":measured_n,
             "taxa_QC":int(measured_taxa),
             "n_after_exact_12_climate_and_min5":len(complete),
             "taxa_after_min5":int(complete.accepted_key.nunique()),
             "candidate_spine_and_phyllary_proxy_only":True}
    if complete.accepted_key.nunique()<MIN_TAXA:
        return None,{"status":"INSUFFICIENT_COMMON_BRACHT_TAXA_FOR_SCALE_TEST",
                      **receipt}
    return complete,{"status":"SOURCE_SUPPORT_SATISFIED",**receipt}


def normalize(frame):
    traits=frame[list(COLUMNS)].to_numpy(float)
    scale=traits.std(axis=0)
    if min(scale)<1e-9 or not np.isfinite(scale).all():
        raise ValueError("INVALID_BRACHT_SCALE")
    return (traits-traits.mean(axis=0))/scale


def group_covariances(values, groups):
    labels,inv=np.unique(np.asarray(groups).astype(str),return_inverse=True)
    counts=np.bincount(inv)
    means=np.zeros((len(labels),values.shape[1]),float)
    np.add.at(means,inv,values)
    means/=counts[:,None]
    within_residual=values-means[inv]
    within=(within_residual/counts[inv,None]).T@within_residual/len(labels)
    between_c=means-means.mean(axis=0)
    between=between_c.T@between_c/len(labels)
    return within,between,means,inv,counts


def bract_edges(within,among):
    index={name:[COLUMNS.index(col) for col in cols] for name,cols in FEATURES.items()}
    names=list(FEATURES)
    edges=[]
    for i,left in enumerate(names):
        for right in names[i+1:]:
            w=association(within,index[left],index[right])
            a=association(among,index[left],index[right])
            edges.append({"left":left,"right":right,
                          "within_RV":float(w[0]),"among_RV":float(a[0]),
                          "delta_RV":float(a[0]-w[0]),
                          "not_genetic_modularity":True})
    return edges


def triple_configuration(values,groups):
    original=values[:,[COLUMNS.index(c) for c in TRIAD]]
    within,between,means,inv,counts=group_covariances(original,groups)
    bits=(original-means[inv]>0).astype(int)
    code=bits[:,0]*4+bits[:,1]*2+bits[:,2]
    by_taxon=np.zeros((len(counts),8),float)
    np.add.at(by_taxon,(inv,code),1.)
    by_taxon/=counts[:,None]
    within_freq=by_taxon.mean(axis=0)
    means_bits=(means-np.median(means,axis=0)>0).astype(int)
    among_freq=np.bincount(means_bits[:,0]*4+means_bits[:,1]*2+
                            means_bits[:,2],minlength=8).astype(float)/len(counts)
    return {"triad":list(TRIAD),"n_taxa":len(counts),
            "within":maxent_from_observed(within_freq),
            "among":maxent_from_observed(among_freq),
            "triple_is_visible_phyllary_projection_NOT_genetic_G123":True}


def evaluate(native,all_env,scenario_selected):
    cohort,receipt=valid_source(native,all_env)
    if cohort is None:
        return receipt
    v=normalize(cohort)
    groups=cohort.accepted_key
    w,a,*_=group_covariances(v,groups)
    result={"status":"SOURCE_LINKED_BRACHT_PROXY_EXPLORATORY",
            **receipt,"raw_bract_pair_relations":bract_edges(w,a),
            "raw_visible_triad_beyond_pairwise":triple_configuration(v,groups),
            "adjusted_scenarios":{},
            "distinct_population_from_eight_core_constructs":True}
    from analysis.v4.s1_s2_extreme_lego import residualize
    for name in ("annual_BIO12_elevation","warmest_quarter_BIO18_elevation"):
        model=scenario_selected[name]
        if model["status"]!="SELECTED":
            result["adjusted_scenarios"][name]={"status":model["status"]}
            continue
        predictors=model["selected"]
        try:
            adj,*_=residualize(cohort,v,predictors)
            within,among,*_=group_covariances(adj,groups)
            result["adjusted_scenarios"][name]={
                "status":"DESCRIPTIVE",
                "predictors":predictors,
                "pair_relations":bract_edges(within,among),
                "three_way_visible_state":triple_configuration(adj,groups),
                "n_head_observations":len(cohort)}
        except (ValueError,np.linalg.LinAlgError) as exc:
            result["adjusted_scenarios"][name]={"status":"NOT_IDENTIFIABLE",
                                                 "reason":str(exc)[:180]}
    return result

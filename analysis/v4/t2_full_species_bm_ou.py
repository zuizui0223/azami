"""T2 census-wide BM/OU/WN screen: every available S1 angle tip, plus independent nuclear controls.

Primary census is ALL 191 concept-clean S1 photo-angle tips mapped to each of
52 NON-INDEPENDENT placement hypotheses. They are NOT a Cirsium nuclear tree.
NA63 published nuclear climate niche (all 63) and same measured S1 NA42 (all
42 overlapping) are independently-rooted regional comparisons.
This does not infer independent nodding origins or adaptive optimum shifts.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_signal import (
    verify_source, taxon_sets, model_design,
)
from analysis.v4.s1_continuous_phylogenetic_history import (
    audited_pruned_tree, shared_path_matrix,
)
from analysis.v4.t2_public_niche_bm_ou import _covariance
from analysis.v4.s1_continuous_na_nuclear_signal import source_rows
from analysis.v4.t2_pgls_placement_screen import (
    download_author_states, read_author_nuclear_tree,
    AUTHOR_NUCLEAR_TREE_GIT_BLOB, AUTHOR_NUCLEAR_TREE_PATH, SOURCE_DOI,
)
from analysis.v4.t2_phylogenetic_climate_null import (
    fetch_author, read_published_climate, ENV_PATH, ENV_GIT_BLOB,
)

FIXED_S1 = ("chelsa_bio01", "chelsa_bio04", "chelsa_bio12",
            "chelsa_bio15", "chelsa_rsds_mean")
PRIMARY_S1_N = 191
PLACEMENT_N = 52
NUCLEAR_N = 63
NUCLEAR_PHOTO_N = 42
NUCLEAR_AXES = ("BIO1", "BIO4", "BIO12", "BIO15")


def fit_three_models(y, X, depth, mrca, nugget=0.05, grid_points=45):
    """Gaussian ML with a shared design X: WN, BM, and nested one-optimum OU.

    OU root is assumed at its stationary optimum, root-state uncertainty is
    NOT propagated. Nugget is a fixed covariance sensitivity, not learned
    measurement error. Non-ultrametric trees are accepted as phylograms, and
    branch units are retained (not reinterpreted as dates).
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    depth = np.asarray(depth, dtype=float)
    mrca = np.asarray(mrca, dtype=float)
    n, p = X.shape
    if (y.shape != (n,) or depth.shape != (n,) or mrca.shape != (n, n)
            or n < max(9, p+5) or not np.isfinite(y).all()
            or not np.isfinite(X).all() or not np.isfinite(mrca).all()
            or np.any(depth <= 0) or not (0 <= nugget <= 1)
            or np.linalg.matrix_rank(X) != p or np.std(y) < 1e-10):
        raise ValueError("Unidentified shared-tip BM/OU/WN dataset")
    yy = (y-np.mean(y))/np.std(y, ddof=1)
    def nll_from_cov(C):
        try:
            ch = cho_factor(C, lower=True, check_finite=False)
            ciX = cho_solve(ch, X, check_finite=False)
            ciy = cho_solve(ch, yy, check_finite=False)
            A = X.T @ ciX
            b = X.T @ ciy
            beta = np.linalg.solve(A,b)
            z = yy-X@beta
            sigma2 = float(z @ cho_solve(ch,z,check_finite=False)/n)
            if not np.isfinite(sigma2) or sigma2 < 1e-12:
                return float("inf")
            logdet = 2*np.log(np.diag(ch[0])).sum()
            return float(.5*(n*(np.log(2*np.pi)+1+np.log(sigma2))+logdet))
        except (ValueError, np.linalg.LinAlgError):
            return float("inf")
    # Independent tips are a diagnostic null, not an ancestral model.
    nll_wn = nll_from_cov(np.eye(n))
    def fun(alpha):
        return nll_from_cov(_covariance(depth,mrca,alpha,nugget))
    nll_bm = fun(0.)
    grid = np.r_[0.,np.geomspace(1e-4,250,grid_points)]
    values = np.asarray([fun(a) for a in grid],dtype=float)
    if not np.isfinite(values).any() or not np.isfinite(nll_bm):
        raise ValueError("No finite Brownian or OU likelihood")
    idxs = np.argsort(values)[:4]
    best = [(float(values[i]),float(grid[i])) for i in idxs]
    for i in idxs:
        if i < 2 or i >= len(grid)-1:
            continue
        result = minimize_scalar(lambda x:fun(float(np.exp(x))),
                  bounds=(np.log(grid[i-1]),np.log(grid[i+1])),
                  method="bounded",options={"xatol":5e-4})
        if result.success and np.isfinite(result.fun):
            best.append((float(result.fun),float(np.exp(result.x))))
    best.append((float(nll_bm),0.0))
    nll_ou, alpha = min(best, key=lambda x:x[0])
    if nll_ou > nll_bm + 1e-8:
        raise ValueError("OU not nested in BM on the same measured tips")
    def aicc(nll,k):
        if n <= k+1:
            raise ValueError("AICc undefined at given sample size")
        return float(2*nll+2*k+(2*k*(k+1))/(n-k-1))
    k = p+1 # p regression parameters + diffusion/normal variance
    aics = {"WN":aicc(nll_wn,k), "BM":aicc(nll_bm,k),
            "OU1":aicc(nll_ou,k+1)}
    rank = sorted(aics,key=aics.get)
    return {
        "n_tips":n, "n_regression_coefficients":p,
        "nugget_assumption":nugget,
        "loglik":{"WN":-float(nll_wn),"BM":-float(nll_bm),
                  "OU1":-float(nll_ou)},
        "aicc":aics, "preferred":rank[0],
        "delta_ou_minus_bm_aicc":aics["OU1"]-aics["BM"],
        "delta_wn_minus_bm_aicc":aics["WN"]-aics["BM"],
        "OU_alpha_times_median_root_tip":float(alpha),
        "OU_optimum_at_upper_search_boundary":bool(alpha >= 249),
        "single_optimum_not_multi_regime":True
    }


def fit_on_tree(tree,rows,outcome,adjusted,nugget=.05):
    ordered,X = model_design(rows,adjusted)
    values = np.asarray([float(r[outcome]) for r in ordered])
    names = [r["tip"] for r in ordered]
    pruned,audit = audited_pruned_tree(tree,names)
    C = shared_path_matrix(pruned,names)[0]
    result = fit_three_models(values,X,np.diag(C),C,nugget=nugget)
    result["zero_length_branches_regularized"] = audit["zero_length_edges_regularised"]
    return result


def analyze_s1(source,trees):
    rows=verify_source(source)
    sets=taxon_sets(rows)
    if ({k:len(v) for k,v in sets.items()} != {
           "original194":194,"collision_free191":191,
           "collision_free_at_least5":171,"collision_free_at_least10":143
        } or len(trees)!=PLACEMENT_N):
        raise ValueError("S1 all-species census or 52-placement source changed")
    primary=sets["collision_free191"]
    results=[]
    for i,(tree_label,tree) in enumerate(trees):
        for adjusted in (False,True):
            fitted=fit_on_tree(tree,primary,"image_mean_deg",adjusted)
            results.append({"tree":tree_label,"tree_index":i,
                            "cohort":"ALL_191_UNAMBIGUOUS_TREE_MAPPED",
                            "response":"S1_mean_EXIF_image_vertical_angle_deg",
                            "covariates":"five_S1_fixed_climates" if adjusted else "intercept_only",
                            **fitted})
    # Confirm effect of three alias collisions and image depth on SAME frozen
    # leading placement (never silently promote these to the primary analysis).
    sensitivity=[]
    for cohort_name in ("original194","collision_free_at_least10"):
        fitted=fit_on_tree(trees[0][1],sets[cohort_name],
                           "image_mean_deg",False)
        sensitivity.append({"cohort":cohort_name,"tree":trees[0][0],
                            "response":"S1_mean","covariates":"intercept_only",
                            **fitted})
    for metric in ("image_mean_deg","image_median_deg"):
        fitted=fit_on_tree(trees[0][1],primary,metric,False,nugget=.25)
        sensitivity.append({"cohort":"collision_free191",
                            "tree":trees[0][0],
                            "response":metric,"covariates":"intercept_only",
                            **fitted})
    grouped={}
    for c in ("intercept_only","five_S1_fixed_climates"):
        group=[r for r in results if r["covariates"]==c]
        grouped[c]={
            "n_nonindependent_placements":len(group),
            "model_winners":dict(Counter(x["preferred"] for x in group)),
            "ou_minus_bm_aicc_range":[min(x["delta_ou_minus_bm_aicc"] for x in group),
                                     max(x["delta_ou_minus_bm_aicc"] for x in group)],
            "n_ou_alpha_at_upper_search_boundary":sum(
                  x["OU_optimum_at_upper_search_boundary"] for x in group),
        }
    return {
        "source":{"S1_native_total_photo_eligible":170066,
                  "S1_native_accepted_keys":244,
                  "mapped_to_historical_placements":194,
                  "concept_clean_primary":191,
                  "unmapped_to_placements":50,
                  "mapped_but_concept_collided":3,
                  "primary_includes_every_eligible_noncollided_tree_tip":True,
                  "source_artifact_sha256":hashlib.sha256(Path(source).read_bytes()).hexdigest()},
        "primary_52_placement_summary":grouped,
        "primary_each_placement":results,
        "same_first_tree_measurement_sampling_sensitivities":sensitivity,
        "phylogeny_status":"GRAFTED_GBOTB_VPhyloMaker_NOT_INDEPENDENT_CIRSIUM_NUCLEAR",
    }


def analyze_independent_na63():
    author_states,trait_sha=download_author_states()
    tree_raw,tree_sha=fetch_author(AUTHOR_NUCLEAR_TREE_PATH,AUTHOR_NUCLEAR_TREE_GIT_BLOB)
    climate_raw,climate_sha=fetch_author(ENV_PATH,ENV_GIT_BLOB)
    climate=read_published_climate(climate_raw)
    tree=read_author_nuclear_tree(tree_raw,author_states)
    names=sorted(author_states)
    if len(names)!=NUCLEAR_N or len(climate)!=NUCLEAR_N:
        raise ValueError("63-tip nuclear author sample contract changed")
    phy, _=audited_pruned_tree(tree,names)
    C=shared_path_matrix(phy,names)[0]
    results={}
    for axis in NUCLEAR_AXES:
        yy=[climate[t][axis] for t in names]
        result=fit_three_models(yy,np.ones((len(names),1)),np.diag(C),C)
        results[axis]=result
    return {
        "n_taxa":NUCLEAR_N,
        "n_botanical_nodding":sum(author_states.values()),
        "independent_published_nuclear_tree":True,
        "source_doi":SOURCE_DOI,
        "source_sha256":{"traits":trait_sha,"tree":tree_sha,"climate":climate_sha},
        "climate_axis_units":"published-source units; BIO1 may be scaled; do not assume degrees Celsius",
        "results":results,
        "branch_lengths":"original author ultrametric tree; numerical scale not automatically calendar time",
        "trait_type":"present_day_species_environment_not_biological_head_orientation",
    }


def analyze_nuclear_s1_overlap(joint, tree):
    rows=source_rows(joint)
    if len(rows)!=NUCLEAR_PHOTO_N:
        raise ValueError("All 42 exact North American S1 photo taxa required")
    # Same 42 concepts and same PHOTO ANGLE on one original published nuclear tree.
    clean=[{"tip":r["tip"],"image_mean_deg":r["image_median_deg"],
            "image_median_deg":r["image_median_deg"],
            "S1_five_climate_means":r["S1_five_climate_means"]}
           for r in rows]
    results={}
    for adjusted in (False, True):
        r=fit_on_tree(tree,clean,"image_mean_deg",adjusted)
        results["five_fixed_climate_covariates" if adjusted else "intercept_only"]=r
    return {"n_photo_taxa":NUCLEAR_PHOTO_N,
            "source":"exact S1 continuous photographic angle matched to genuine author NA nuclear tree",
            "angle_is_gravity_referenced":False,"models":results}


def run(source,joint):
    tree_states,traits_sha=download_author_states()
    tree_bytes,tree_sha=fetch_author(AUTHOR_NUCLEAR_TREE_PATH,AUTHOR_NUCLEAR_TREE_GIT_BLOB)
    original_nuclear=read_author_nuclear_tree(tree_bytes,tree_states)
    trees=restore_trees()
    return {
        "schema":"azami_v4_T2_ALL_AVAILABLE_S1_BM_OU_VS_WN_v1",
        "status":"EMPIRICAL_ALL_AVAILABLE_MEASURED_S1_IMAGE_ANGLE_PHYLOGENETIC_SCREEN",
        "S1_all_matched":analyze_s1(source,trees),
        "published_independent_NA63_niche":analyze_independent_na63(),
        "independent_nuclear_S1_NA42_exact":analyze_nuclear_s1_overlap(joint,original_nuclear),
        "model":{
           "WN":"phylogenetically independent Gaussian tip residuals",
           "BM":"Brownian diffusion on source branch lengths",
           "OU1":"single optimum Ornstein-Uhlenbeck on source branch lengths",
           "comparison":"same tip cohort and same regressors; Gaussian ML and AICc",
           "model_nugget":"fixed 0.05 primary, 0.25 first-tree sensitivity; not inferred image error",
           "not_multiregime_selection_test":True,
        },
        "scientific_limits":{
           "S1_244_independent_nuclear_tree_complete":False,
           "placement_tree_does_not_prove_nuclear_ancestry":True,
           "photo_angle_is_not_gravitational_head_nodding":True,
           "S2_two_species_not_a_phylogenetic_trajectory":True,
           "no_evidence_of_adaptive_radiation_from_BM_OU_winner_alone":True,
           "no_independent_gain_counts_from_continuous_photo_change":True,
           "branch_length_and_lineage_uncertainty_unmodeled":True,
           "environment_tip_centroids_not_ancestral_rain":True
        }
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True)
    p.add_argument("--joint",required=True)
    p.add_argument("--output",required=True)
    args=p.parse_args()
    result=run(args.source,args.joint)
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,
                              allow_nan=False)+"\n",encoding="utf-8")
    print("T2_ALL_AVAILABLE_BM_OU_SUMMARY="+json.dumps({
      "s1":result["S1_all_matched"]["primary_52_placement_summary"],
      "NA63":{k:v["preferred"] for k,v in result["published_independent_NA63_niche"]["results"].items()},
      "NA42":{k:v["preferred"] for k,v in result["independent_nuclear_S1_NA42_exact"]["models"].items()}
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

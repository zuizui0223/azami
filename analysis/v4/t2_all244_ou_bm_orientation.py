"""Whole-S1-coverage continuous image-angle WN/BM/OU1 sensitivity.

Starts with every one of 244 original native source concepts (170066 source photos).
Fits all tree-covered taxa, with 191 uncollided concepts primary, n>=10
sensitivity, and a genuinely published 42-tip NA nuclear benchmark.
Uses EAzami's original finite-root OU covariance, not 52 independent histories.
Neither an OU preference nor alpha establishes adaptive radiation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

from analysis.v4.joint_s1_t2_phylogeny import (
    ENDPOINT, EXPOSURES, MODEL_PATH, COVERAGE_PATH, MAPPING_PATH,
    actual_native_roster, source_concept_index,
)
from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_history import (
    verified_species, audited_pruned_tree, shared_path_matrix, SHA,
)
from analysis.v4.t2_public_niche_bm_ou import _covariance
from analysis.v4.t2_source_blocked_five_climate import author_nuclear_tree

SOURCE_244_N=244
SOURCE_IMAGES_N=170066
TREE_194_N=194
CLEAN_191_N=191
TEN_143_N=143
NUCLEAR_NA42_N=42
COLLIDED_TIPS={"Cirsium_parryi","Cirsium_setidens","Cirsium_yezoense"}
S2_FOCAL={"Cirsium_arvense","Cirsium_vulgare"}
MATRIX_CASES=(
    ("source194","mean",True,0.05),
    ("clean191","mean",False,0.05),
    ("clean191","mean",True,0.05),
    ("clean191","median",True,0.05),
    ("clean191","mean",True,0.25),
    ("n10_143","mean",True,0.05),
)
ALPHA_GRID=np.r_[0., np.geomspace(1e-4,250,24)]


def source_roster(table: Path):
    """Never leak original observation/photographer/coordinate level data."""
    if hashlib.sha256(table.read_bytes()).hexdigest()!=SHA:
        raise ValueError("Native frozen photo source digest mismatch")
    models=json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    coverage=json.loads(COVERAGE_PATH.read_text(encoding="utf-8"))[ENDPOINT]
    source=actual_native_roster(pd.read_csv(table,low_memory=False),models,coverage)
    mapped=source_concept_index(json.loads(MAPPING_PATH.read_text(encoding="utf-8")))
    primary=verified_species(table)
    if (len(source)!=SOURCE_244_N or
        sum(row["n"] for row in source.values())!=SOURCE_IMAGES_N or
        len(primary)!=TREE_194_N or len({x["accepted_key"] for x in primary})!=TREE_194_N):
        raise ValueError("244 native and 194 placed source cohort changed")
    if {r["tip"] for r in primary if r["collision"]}!=COLLIDED_TIPS:
        raise ValueError("Source-name collision set drifted")
    mapped_rows={str(r["accepted_key"]):r for r in primary}
    roster=[]
    for key,v in sorted(source.items()):
        original_names=sorted({x.replace(" ","_") for x in v["source_names"]})
        m=mapped.get(key,{})
        strict_clean=len(original_names)==1 and not v["concept_collision_in_source"]
        placed=mapped_rows.get(key)
        if placed and placed["tip"] not in set(m.get("tree_matches",set())):
            raise ValueError("Wrong tree tip assigned to WCVP source concept")
        status=("collision_excluded" if not strict_clean
                else "tree_covered" if placed
                else "no_verified_placement")
        roster.append({
            "accepted_key":key,
            "source_taxon_names":"|".join(original_names),
            "n_source_orientation_photos":v["n"],
            "image_angle_mean_deg":v["mean_image_angle_deg"],
            "image_angle_median_deg":v["median_image_angle_deg"],
            "source_name_collision":not strict_clean,
            "tree_tip":placed["tip"] if placed else "",
            "tree194_covered":placed is not None,
            "strict_phylogenetic_eligible":bool(placed and strict_clean),
            "phylogenetic_coverage_status":status,
            "s2_focal":bool(S2_FOCAL.intersection(original_names)),
            **{e:float(v["mean_frozen_five_exposures"][e]) for e in EXPOSURES},
        })
    frame=pd.DataFrame(roster)
    if (len(frame)!=244 or int(frame.tree194_covered.sum())!=194 or
        int(frame.strict_phylogenetic_eligible.sum())!=191):
        raise ValueError("244/194/191 full-taxon ledger invariant failed")
    s2=frame.loc[frame.s2_focal].copy()
    if len(s2)!=2 or set(s2.source_taxon_names)!=S2_FOCAL:
        raise ValueError("S2 original source species do not both map uniquely into S1")
    clean=[r for r in primary if not r["collision"]]
    ten=[r for r in clean if r["n_eligible_S1_images"]>=10]
    if len(clean)!=191 or len(ten)!=143:
        raise ValueError("Original original-photo precision cohort drifted")
    return frame,{"source194":primary,"clean191":clean,"n10_143":ten}


def regression_design(rows,adjust_climate):
    y=np.array([float(r["mean_S1_image_axis_deg"]) for r in rows])
    if not np.isfinite(y).all():
        raise ValueError("Nonfinite source angles")
    X=np.ones((len(rows),1))
    if adjust_climate:
        original=np.array([[float(r["S1_five_climate_means"][e]) for e in EXPOSURES]
                           for r in rows],dtype=float)
        sd=original.std(axis=0)
        if not np.isfinite(original).all() or np.any(sd<=1e-10):
            raise ValueError("Bad original five-climate design")
        X=np.column_stack([np.ones(len(rows)),(original-original.mean(axis=0))/sd])
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise ValueError("Rank deficient original five-climate design")
    return X


def _design_and_values(rows,phenotype,adjust_climate):
    if phenotype not in ("mean","median"):
        raise ValueError("Unknown measured phenotype")
    key="mean_S1_image_axis_deg" if phenotype=="mean" else "median_S1_image_axis_deg"
    y=np.array([float(r[key]) for r in rows],float)
    scale=float(y.std(ddof=1))
    if not np.isfinite(scale) or scale<=1e-9:
        raise ValueError("Measured photo angle is constant")
    X=regression_design(rows,adjust_climate)
    return (y-float(y.mean()))/scale,X,scale


def _gls(y,X,V):
    """Gaussian ML covariance-profile at fixed phylogenetic kernel."""
    n=len(y)
    try:
        ch=cho_factor(V,lower=True,check_finite=False)
        viX=cho_solve(ch,X,check_finite=False)
        viY=cho_solve(ch,y,check_finite=False)
        info=X.T@viX
        beta=np.linalg.solve(info,X.T@viY)
        r=y-X@beta
        sigma2=float(r@cho_solve(ch,r,check_finite=False)/n)
        if sigma2<=1e-14 or not np.isfinite(sigma2):
            return None
        ll=-.5*(n*(np.log(2*np.pi)+1+np.log(sigma2))
                +2*np.log(np.diag(ch[0])).sum())
        beta_cov=sigma2*np.linalg.inv(info)
        return dict(loglik=float(ll),beta=beta,cov=beta_cov,sigma2=sigma2)
    except (ValueError,np.linalg.LinAlgError):
        return None


def _aicc(loglik,n,k):
    if n<=k+1:
        raise ValueError("Not enough species for AICc")
    return float(-2*loglik+2*k+2*k*(k+1)/(n-k-1))


def evaluate_models(y,X,C,phenotype_scale,nugget,*,alpha_grid=ALPHA_GRID):
    """WN/BM/OU1 with fixed-root EAzami kernel and identical taxon+X cohort."""
    n=len(y)
    if C.shape!=(n,n) or len(y)!=len(X) or X.ndim!=2 or not (0<=nugget<=1):
        raise ValueError("Incomplete matched species-level BM/OU inputs")
    C=(C+C.T)/2
    depth=np.diag(C).copy()
    if np.any(depth<=0) or np.min(np.linalg.eigvalsh(C))<=0:
        raise ValueError("No positive rooted BM covariance")
    I=np.eye(n)
    wn=_gls(y,X,I)
    bm=_gls(y,X,_covariance(depth,C,0.,nugget))
    if wn is None or bm is None:
        raise ValueError("Unidentified source WN/BM model")
    def eval_ou(alpha):
        return _gls(y,X,_covariance(depth,C,alpha,nugget))
    hits=[(float(alpha),val) for alpha in alpha_grid
          if (val:=eval_ou(float(alpha))) is not None]
    if not hits or hits[0][0]!=0:
        raise ValueError("OU model lost its nested alpha=0 BM limit")
    ix=max(range(len(hits)),key=lambda i:hits[i][1]["loglik"])
    opt_alpha,best=hits[ix]
    if 0<ix<len(hits)-1:
        lo,hi=hits[ix-1][0],hits[ix+1][0]
        if lo>0:
            opt=minimize_scalar(lambda la: (
                -fit["loglik"] if (fit:=eval_ou(float(np.exp(la)))) else 1e30),
                bounds=(float(np.log(lo)),float(np.log(hi))),method="bounded",
                options={"xatol":1e-4})
            if opt.success:
                fit=eval_ou(float(np.exp(opt.x)))
                if fit is not None and fit["loglik"]>best["loglik"]:
                    opt_alpha,best=float(np.exp(opt.x)),fit
    if best["loglik"]+1e-7<bm["loglik"]:
        raise ValueError("OU alpha=0 is nested BM but worse optimization")
    k=X.shape[1]+1
    aiccs={"WN":_aicc(wn["loglik"],n,k),
           "BM":_aicc(bm["loglik"],n,k),
           "OU1":_aicc(best["loglik"],n,k+1)}
    selected=min(aiccs,key=aiccs.get)
    # Biological degrees per +1 SD of annual precipitation on the original
    # S1 taxon mean; this remains a spatial observational coefficient.
    coef=None
    if X.shape[1]==6:
        original={"WN":wn,"BM":bm,"OU1":best}[selected]
        j=3 # [intercept, BIO1, BIO4, BIO12, BIO15, RSDS]
        coef=float(original["beta"][j]*phenotype_scale)
    return {
        "n_species":n,"n_covariates":X.shape[1]-1,
        "nugget_relative_to_median_root_depth":float(nugget),
        "WN_AICc":aiccs["WN"],"BM_AICc":aiccs["BM"],
        "OU1_AICc":aiccs["OU1"],
        "best_model_AICc":selected,
        "OU1_minus_BM_AICc":float(aiccs["OU1"]-aiccs["BM"]),
        "OU1_alpha_times_median_root_to_tip_branch_length":opt_alpha,
        "OU1_alpha_grid_edge":bool(opt_alpha>=float(alpha_grid[-1])-.01),
        "selected_model_BIO12_deg_per_sample_SD":coef,
        "OU1_is_single_optimum_no_optimum_shift_test":True,
        "branch_units_are_phylogenetic_placement_not_calendar_years":True,
    }


def white_noise_full_244(frame):
    """Actual full-cohort fit without pretending 50 unmatched taxa are in a tree."""
    frame=frame.sort_values("accepted_key")
    arr=frame[list(EXPOSURES)].to_numpy(float)
    z=(arr-arr.mean(axis=0))/arr.std(axis=0)
    X=np.column_stack((np.ones(len(frame)),z))
    out={}
    for v in ("image_angle_mean_deg","image_angle_median_deg"):
        raw=frame[v].to_numpy(float)
        scale=float(raw.std(ddof=1))
        yy=(raw-raw.mean())/scale
        fit=_gls(yy,X,np.eye(len(frame)))
        if fit is None:
            raise ValueError("Cannot fit 244-taxon WN source observational baseline")
        out[v]={
            "n_accepted_keys":len(frame),
            "full_five_climate_WN_AICc":_aicc(fit["loglik"],len(frame),X.shape[1]+1),
            "BIO12_image_angle_degrees_per_sample_SD":float(fit["beta"][3]*scale),
            "not_phylogenetic_or_spatially_corrected":True,
            "not_original_observation_level_S1_reestimate":True,
        }
    return out


def analyse(table:Path, tree_limit=52):
    frame,groups=source_roster(table)
    if not (1<=tree_limit<=52):
        raise ValueError("Specify between 1 and 52 original placement scenarios")
    source194=groups["source194"]
    lookup={r["tip"]:i for i,r in enumerate(source194)}
    cohort_views={}
    for label,rows in groups.items():
        names=[r["tip"] for r in rows]
        cohort_views[label]=([lookup[n] for n in names],rows)
    status=frame.groupby("phylogenetic_coverage_status").size().to_dict()
    result={
       "schema":"azami_V4_ALL244_S1_S2_IMAGE_ANGLE_PHYLOGENETIC_OU_BM_v1",
       "source_original_SHA256":SHA,
       "all_source_eligible":244,"original_source_images":170066,
       "original_tree_covered":194,"main_collision_free":191,
       "high_depth_n10":143,
       "all244_coverage_by_reason":{k:int(v) for k,v in status.items()},
       "s2_exact_source_species":[
           {"accepted_key":str(x["accepted_key"]),"taxon":x["source_taxon_names"],
            "in_tree194":bool(x["tree194_covered"]),
            "in_primary_clean191":bool(x["strict_phylogenetic_eligible"])}
           for x in frame.to_dict("records") if x["s2_focal"]],
       "all244_equal_taxon_WN_five_climate_baselines":white_noise_full_244(frame),
       "primary_phenotype":"image_vertical_angle_in_degrees_not_gravity_referenced_posture",
       "auxiliary_phenotype":"taxon_median_same_image_vertical_angle",
       "five_fixed_climate_covariates":list(EXPOSURES),
       "model_candidates":["WN","BM","OU1"],
       "model_input_phylogenies":"52 GBOTB/VPhyloMaker2 placement sensitivity, no independent global nuclear species phylogeny",
       "model_design":"ML source species means; optional simultaneous five-climate design; EAzami finite-root OU covariance; fixed tip nugget .05/.25, AICc; sensitivity across same 52 non-independent graft hypotheses",
       "claims_not_supported":{
           "same_full_244_in_phylogenetic_comparison":True,
           "OU1_support_proves_adaptive_radiation":True,
           "continuous_image_angle_is_botanical_downward_state":True,
           "grafted_52_placement_scenarios_are_independent_phylogenies":True,
           "chronogram_diversification_rate_estimated":True,
           "rain_caused_orientation_transitions":True,
       },
    }
    model_rows=[]
    for idx,(label,tree) in enumerate(restore_trees()[:tree_limit]):
        names=[r["tip"] for r in source194]
        pruned,_=audited_pruned_tree(tree,names)
        C,_,_,_=shared_path_matrix(pruned,names)
        for cohort,phenotype,adjust,nugget in MATRIX_CASES:
            subix,rows=cohort_views[cohort]
            subC=C[np.ix_(subix,subix)]
            yy,X,scale=_design_and_values(rows,phenotype,adjust)
            model_rows.append({
                "tree_scenario_id":idx,"tree_scenario":label,
                "cohort":cohort,"endpoint":phenotype,
                "adjust_fixed_five_climate":bool(adjust),
                **evaluate_models(yy,X,subC,scale,nugget),
            })
        if idx==0 or (idx+1)%10==0:
            print("ALL244_OU_BM_SCENARIOS_DONE="+str(idx+1),flush=True)
    summary={}
    for config,block in pd.DataFrame(model_rows).groupby(
         ["cohort","endpoint","adjust_fixed_five_climate",
          "nugget_relative_to_median_root_depth"],dropna=False):
        name="|".join(map(str,config))
        summary[name]={
           "n_nonindependent_placement_scenarios":len(block),
           "best_model_counts":{
                v:int((block.best_model_AICc==v).sum()) for v in ("WN","BM","OU1")},
           "OU1_minus_BM_AICc_min":float(block.OU1_minus_BM_AICc.min()),
           "OU1_minus_BM_AICc_median":float(block.OU1_minus_BM_AICc.median()),
           "OU1_minus_BM_AICc_max":float(block.OU1_minus_BM_AICc.max()),
           "OU1_alpha_at_grid_upper_count":int(block.OU1_alpha_grid_edge.sum()),
           "not_independent_tree_replicates":True,
        }
    result["comparison_by_case"]=summary
    result["n_placement_scenarios_used"]=tree_limit
    result["n_total_model_cases"]=len(model_rows)
    return frame,pd.DataFrame(model_rows),result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native-table",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--tree-limit",type=int,default=52)
    args=p.parse_args()
    frame,fits,summary=analyse(args.native_table,args.tree_limit)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    frame.to_csv(args.out_dir/"S1_ALL244_TAXON_COVERAGE.csv",index=False)
    fits.to_csv(args.out_dir/"S1_191_PHYLOGENETIC_WN_BM_OU1_ALL_PLACEMENTS.csv",index=False)
    (args.out_dir/"S1_ALL244_NUCLEAR_COVERAGE_AND_BM_OU_SUMMARY.json").write_text(
       json.dumps(summary,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n",
       encoding="utf-8")
    print("S1_ALL244_COMPLETE_WN_BM_OU1_SUMMARY="+json.dumps({
        "all_244":summary["all_source_eligible"],
        "tree_194":summary["original_tree_covered"],
        "clean_191":summary["main_collision_free"],
        "n_cases":summary["n_total_model_cases"],
        "core_case":summary["comparison_by_case"].get("clean191|mean|True|0.05",{})
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

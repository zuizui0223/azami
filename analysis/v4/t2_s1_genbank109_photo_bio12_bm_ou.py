"""BM, OU1, stationary OU, and WN on OWN GenBank ITS 109-tip S1 panel.

This analysis exclusively uses the genuine S1 photos' continuous image-axis
orientation (species mean and median) and original CHELSA BIO12. Binary U/D
obtained from flora/taxonomy or the retired classifier is NOT used.

109-tip ITS maximum-likelihood tree is built from GenBank records selected
from ALL 251 original S1 named units. Single multicopy ribosomal ITS is NOT a
multilocus species tree; the display root is unknown; ITS branches are in
substitutions per site, not years. OU alpha cannot be read as Myr^-1, and an
OU optimum is not proof of selection, especially when alpha is on a bound.

A fixed 5% *hypothetical* tip nugget is a transparent numerical/measurement
sensitivity, not a measured calibration; 0% is tested independently.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path

from Bio import Phylo
import numpy as np
from scipy import linalg, optimize, stats

from analysis.v4.t2_s1_genbank_its_original_photo_angle_rain import (
    load_evidence, source_design, RESPONSES, EXPOSURES
)

ROOTS=("midpoint_root_unverified","original_Newick_display_root_unverified")
NUGGETS=(0.05,0.0)
ALPHA_BOUNDS=(0.001,300.0)
EXPECTED_TIPS=109


def original_covariances(tree):
    """Source-only phylogenetic shared-path and patristic distances."""
    tips=tree.get_terminals()
    if len(tips)!=EXPECTED_TIPS or len({t.name for t in tips})!=EXPECTED_TIPS:
        raise ValueError("Original GenBank source tree tip scope changed")
    d=tree.depths()
    lengths=np.asarray([float(d[t]) for t in tips])
    if not np.all(np.isfinite(lengths)) or np.any(lengths<=0):
        raise ValueError("Invalid authentic ITS branch lengths")
    scale=float(np.median(lengths))
    n=len(tips);C=np.empty((n,n))
    for i,a in enumerate(tips):
        C[i,i]=lengths[i]/scale
        for j in range(i):
            val=float(d[tree.common_ancestor(a,tips[j])])/scale
            C[i,j]=C[j,i]=val
    dist=np.maximum(C.diagonal()[:,None]+C.diagonal()[None,:]-2*C,0.)
    if np.linalg.eigvalsh(C).min() <= 0:
        raise ValueError("Singular molecular-tree shared-ancestry matrix")
    return [t.name for t in tips],C,dist,scale


def kernel(C,D,model,alpha=None,tip_nugget=0.05):
    """Phylogenetic residual covariance in median-tip-variance units."""
    if not 0<=tip_nugget<=0.15:
        raise ValueError("Measurement nugget outside declared sensitivity")
    n=len(C)
    if model=="WN":return np.eye(n)
    if model=="BM":K=C.copy()
    elif model=="OU1":
        if alpha is None or not np.isfinite(alpha) or alpha<=0:
            raise ValueError("OU alpha must be finite positive")
        K=np.exp(-alpha*D)*(-np.expm1(-2*alpha*C))/(2*alpha)
    elif model=="OU_stationary":
        if alpha is None or not np.isfinite(alpha) or alpha<=0:
            raise ValueError("OU stationary alpha must be finite positive")
        K=np.exp(-alpha*D)
    else:raise ValueError("Unrecognized Gaussian evolutionary model")
    K=(K+K.T)*.5
    median=float(np.median(np.diag(K)))
    if median<=0:raise ValueError("Evolutionary covariance lacks variance")
    return K/median + np.eye(n)*tip_nugget


def gls(y,X,K,include_se=False):
    """Profiled ML Gaussian likelihood; all models use exact SAME tip set."""
    y=np.asarray(y,dtype=float);X=np.asarray(X,dtype=float)
    n=len(y);p=X.shape[1]
    if K.shape!=(n,n) or X.shape[0]!=n or n<=p+2:
        raise ValueError("Invalid same-species phylogeny/trait model dimensions")
    L=linalg.cholesky(K,lower=True,check_finite=True)
    A=linalg.solve_triangular(L,X,lower=True,check_finite=True)
    b=linalg.solve_triangular(L,y,lower=True,check_finite=True)
    if np.linalg.matrix_rank(A)<p:
        raise ValueError("Original climate axes collinear")
    beta=np.linalg.lstsq(A,b,rcond=None)[0]
    rss=float(np.sum((b-A@beta)**2))
    if rss<=1e-14:raise ValueError("No species phenotypic residual variance")
    logL=-0.5*(n*(np.log(2*np.pi)+1+np.log(rss/n))
               +2*sum(np.log(np.diag(L))))
    result={"logLik":float(logL),"n":n,"coef":beta.tolist()}
    if include_se:
        c=np.linalg.inv(A.T@A)*(rss/(n-p))
        se=np.sqrt(np.diag(c));q=float(stats.t.ppf(.975,n-p))
        result["std_error"]=se.tolist()
        result["ci95"]=[[float(x-q*s),float(x+q*s)] for x,s in zip(beta,se)]
        result["two_sided_p"]=[
            float(2*stats.t.sf(abs(x/s),n-p)) for x,s in zip(beta,se)]
    return result


def aicc(logL,n,parameters):
    if n<=parameters+1:raise ValueError("AICc sample size insufficient")
    return float(2*parameters-2*logL+
                 2*parameters*(parameters+1)/(n-parameters-1))


def fit_candidates(y,X,C,D,nugget):
    """Real original tree Gaussian WN/BM/one-optimum OU with global alpha scan."""
    results={}
    cache={}
    for model in ("WN","BM"):
        ll=gls(y,X,kernel(C,D,model,tip_nugget=nugget))
        k=X.shape[1]+1
        results[model]={"AICc":aicc(ll["logLik"],len(y),k),
                        "logLik":ll["logLik"],"n_parameters":k,
                        "mean_intercept_or_regression":ll["coef"]}
    for model in ("OU1","OU_stationary"):
        memo={}
        def fit_at_logalpha(loga):
            if loga not in memo:
                try:
                    a=float(np.exp(loga))
                    memo[loga]=gls(y,X,kernel(C,D,model,a,nugget))
                except (linalg.LinAlgError,ValueError):
                    memo[loga]={"logLik":float("-inf")}
            return memo[loga]["logLik"]
        grid=np.linspace(np.log(ALPHA_BOUNDS[0]),np.log(ALPHA_BOUNDS[1]),45)
        ygrid=np.array([fit_at_logalpha(float(z)) for z in grid])
        if not np.any(np.isfinite(ygrid)):
            raise ValueError("No finite OU likelihood over original ITS branches")
        best_indices=np.argsort(-ygrid)[:4]
        candidates=[float(grid[0]),float(grid[-1])]
        for ix in best_indices:
            lo=float(grid[max(0,ix-1)])
            hi=float(grid[min(len(grid)-1,ix+1)])
            if hi<=lo:continue
            result=optimize.minimize_scalar(
                lambda x:-fit_at_logalpha(float(x)),
                bounds=(lo,hi),method="bounded",options={"xatol":0.0005})
            if result.success:candidates.append(float(result.x))
        chosen=max(candidates,key=fit_at_logalpha)
        alpha=float(np.exp(chosen))
        ll=float(fit_at_logalpha(chosen))
        k=X.shape[1]+2
        results[model]={
            "AICc":aicc(ll,len(y),k),"logLik":ll,"n_parameters":k,
            "mean_intercept_or_regression":memo[chosen]["coef"],
            "alpha_per_median_ITS_root_tip_distance":alpha,
            "alpha_at_or_near_search_boundary":bool(
                alpha < ALPHA_BOUNDS[0]*1.5 or alpha>ALPHA_BOUNDS[1]*.93),
            "single_optimum_not_multi_regime":True,
            "stationary_root_covariance":model=="OU_stationary",
        }
    winner=min(results,key=lambda key:results[key]["AICc"])
    minimum=results[winner]["AICc"]
    for val in results.values():val["delta_AICc"]=float(val["AICc"]-minimum)
    return {"best_AICc":winner,"models":results,
            "same_nuclear_genbank_109_tips_in_all_models":True,
            "OU1_vs_BM_deltaAICc":float(results["OU1"]["AICc"]-results["BM"]["AICc"]),
            "tip_measurement_nugget_fraction_assumed":nugget}


def regress_BIO12_photo(y,photo_X,C,D,nugget):
    fit=fit_candidates(y,photo_X,C,D,nugget)
    slopes={}
    for model,raw in fit["models"].items():
        alpha=raw.get("alpha_per_median_ITS_root_tip_distance")
        info=gls(y,photo_X,kernel(C,D,model,alpha,nugget),include_se=True)
        # fixed S1 design: intercept; zBIO1; zBIO4; BIO12 centered/100mm;
        # zBIO15; zRSDS. This is NEVER a binary flora U/D predictor.
        at=3
        slopes[model]={
            "BIO12_image_angle_degrees_per_100mm":float(info["coef"][at]),
            "BIO12_ci95_degrees_per_100mm":info["ci95"][at],
            "BIO12_two_sided_p_conditional_on_fixed_tree_model":info["two_sided_p"][at],
        }
    return {"model_comparison":fit,"nuclear_five_S1_climate_axis_slopes":slopes,
            "S1_photo_angle_not_true_gravitational_head_posture":True,
            "alpha_and_tree_estimation_uncertainty_not_in_coefficient_ci":True}


def bm_whitened_rain_photo_correlation(y1,y2,C,D,nugget):
    L=linalg.cholesky(kernel(C,D,"BM",tip_nugget=nugget),lower=True)
    intercept=linalg.solve_triangular(L,np.ones(len(y1)),lower=True)
    def resid(y):
        yy=linalg.solve_triangular(L,np.asarray(y),lower=True)
        return yy-intercept*(intercept@yy)/(intercept@intercept)
    a=resid(y1);b=resid(y2)
    corr=float((a@b)/np.sqrt((a@a)*(b@b)))
    return {"r_Brownian_whitened_logBIO12_vs_S1_photo_mean_angle":corr,
            "causal_historical_correlated_evolution_not_tested":True}


def execute(s1,accessions,tree_file,source_receipt,sequence_qc):
    master,meta,tree,source,quality=load_evidence(
        s1,accessions,tree_file,source_receipt,sequence_qc)
    named={row["original_source_name"]:row for row in master}
    result={}
    for mode in ROOTS:
        phy=deepcopy(tree)
        if mode=="midpoint_root_unverified":phy.root_at_midpoint()
        tips,C,D,scale=original_covariances(phy)
        X,_=source_design(named,tips)
        axes={}
        angle_mean=np.array([float(named[t]["mean_image_angle_deg"]) for t in tips])
        angle_median=np.array([float(named[t]["median_image_angle_deg"]) for t in tips])
        lograin=np.array([np.log(float(named[t]["chelsa_bio12"])) for t in tips])
        for nugget in (NUGGETS if mode=="midpoint_root_unverified" else NUGGETS[:1]):
            per_trait={}
            for key,values in (("log_BIO12",lograin),
                               ("S1_photo_mean_angle",angle_mean),
                               ("S1_photo_median_angle",angle_median)):
                per_trait[key]=fit_candidates(values,np.ones((len(tips),1)),C,D,nugget)
            regression=regress_BIO12_photo(angle_mean,X,C,D,nugget)
            regress_median=regress_BIO12_photo(angle_median,X,C,D,nugget)
            result[mode+"__nugget"+str(nugget)]={
                "n_actual_original_ITS_taxa":len(tips),
                "n_S1_original_photo_measurements_on_this_tree":sum(
                    int(named[t]["original_photo_count"]) for t in tips),
                "normalized_median_root_to_tip_ITS_distance":scale,
                "root_is_arbitrary_unverified":True,
                "photo_measurement_nugget_is_assumed_not_measured":nugget,
                "trait_BM_OU_WN":per_trait,
                "photo_mean_slope_five_original_climate_axes":regression,
                "photo_median_slope_five_original_climate_axes":regress_median,
                "two_trait_BM_whitened_correlation":bm_whitened_rain_photo_correlation(
                    lograin,angle_mean,C,D,nugget),
            }
    return {
        "schema":"azami_T2_S1_109_actual_GenBank_ITS_BM_OU1_OUstationary_WN_photo_BIO12_v1",
        "status":"REAL_109_S1_SOURCE_GENBANK_NUCLEAR_ITS_PHYLOGENETIC_NICHE_AND_PHOTO_ANGLE",
        "original_S1_accept_keys":244,
        "original_S1_named_units":251,
        "original_S1_eligible_photographs":170066,
        "actual_genbank_nuclear_ITS_tip_names":109,
        "S1_original_photo_continuous_orientation_only":True,
        "retired_flora_U_D_label_not_used_for_phenotype":True,
        "S1_five_frozen_CHELSA_covariates":list(EXPOSURES),
        "models_are_phenotypic_sensitivities_not_phylo_adaptation_proof":True,
        "OU_single_optimum_parameter_is_in_ITS_substitutions_not_calendar_time":True,
        "ITS_nuclear_multicopy_not_orthology_validated_species_tree":True,
        "fully_independent_transition_and_paleoclimate_tests_completed":False,
        "results":result,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s1-251",required=True,type=Path)
    p.add_argument("--accessions",required=True,type=Path)
    p.add_argument("--its-tree",required=True,type=Path)
    p.add_argument("--source-receipt",required=True,type=Path)
    p.add_argument("--alignment-qc",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    x=p.parse_args()
    report=execute(x.s1_251,x.accessions,x.its_tree,
                   x.source_receipt,x.alignment_qc)
    x.output.parent.mkdir(parents=True,exist_ok=True)
    x.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,
                                   sort_keys=True,allow_nan=False)+"\n",
                        encoding="utf-8")
    for key,case in report["results"].items():
        print("S1_109_REAL_ITS_BM_OU_NICHE_PHOTO="+json.dumps({
            "case":key,
            "model_winners":{t:v["best_AICc"]
                             for t,v in case["trait_BM_OU_WN"].items()},
            "five_climate_angle_beta":{k:v["BIO12_image_angle_degrees_per_100mm"]
                 for k,v in case["photo_mean_slope_five_original_climate_axes"]
                      ["nuclear_five_S1_climate_axis_slopes"].items()},
            "rain_photo_BM_whitened_r":
                case["two_trait_BM_whitened_correlation"]
                    ["r_Brownian_whitened_logBIO12_vs_S1_photo_mean_angle"]
        },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

"""BM / single-optimum OU / independent white noise: BIO12 and original S1 head PHOTO angle.

The SAME original S1 251-source-name photographic measurements (170066 images)
supply log(BIO12) climate niche and unbinned continuous camera-frame head
angle (mean, median). The ONLY tree is the new 109-tip real GenBank nuclear ITS
MAFFT/IQ-TREE ML gene tree (Actions 39914024548), not NA63/Japan38 nor the
52 synthetic generic graft placements. All five S1 original CHELSA covariates
are retained in the angle~BIO12 complementary regression.

OU fixed-root: one global optimum and fixed-at-optimum root, C_ij =
exp(-alpha*(t_i+t_j-2*t_MRCA))*(1-exp(-2alpha*t_MRCA))/(2alpha).
OU stationary-root: one optimum and stationary common root, C_ij =
exp(-alpha*patristic_ij), invariant to unrooted ITS rerooting.
These root conditions are genuinely DIFFERENT statistical hypotheses.
BM C_ij=t_MRCA; WN IID C=I. Absolute ITS branch time cannot be inferred.
A fixed nonestimated 5% tip measurement noise is explicitly stress-tested
against zero noise; it is NOT independently calibrated camera error.

The original camera angle is NOT a gravity-confirmed upright/nodding category.
ITS is multi-copy, not a validated single-copy nuclear species tree.
OU likelihood preference is NOT proof of stabilizing selection or rainfall-
caused head orientation evolution. Photographic species mean BIO12 is a present
sampling niche, not ancestral or full species geographic climatic tolerance.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg, optimize, stats
from analysis.v4 import t2_s1_genbank_its_original_photo_angle_rain as source

TRAITS=("log_BIO12","photo_mean_angle_deg","photo_median_angle_deg")
MODES=("WN","BM","OU_stationary","OU_fixedroot")
N_S1=251
N_TREE=109
SEED=20261023
NUGGET=0.05


def prepare_data(tree, registry, root):
    phy=deepcopy(tree)
    if root=="midpoint":
        phy.root_at_midpoint()
    elif root!="source_Newick":
        raise ValueError("Only original unrooted Newick display and midpoint roots admissible")
    tips=phy.get_terminals()
    names=[x.name for x in tips]
    depths=phy.depths()
    age=np.array([float(depths[x]) for x in tips])
    if len(tips)!=N_TREE or np.any(age<=0):
        raise ValueError("Original real GenBank ITS taxon ages missing")
    med=np.median(age)
    M=np.empty((len(names),len(names)),float)
    for i,a in enumerate(tips):
        M[i,i]=age[i]/med
        for j in range(i):
            M[i,j]=M[j,i]=float(depths[phy.common_ancestor(a,tips[j])])/med
    age=age/med
    dist=np.maximum(age[:,None]+age[None,:]-2*M,0)
    vals={}
    for label in names:
        record=registry[label]
        if record["botanical_gravity_posture_inferred_from_photo_angle"]!="False":
            raise ValueError("Cannot substitute botany or inferred gravity U/D")
    vals["log_BIO12"]=np.log([float(registry[n]["chelsa_bio12"]) for n in names])
    vals["photo_mean_angle_deg"]=np.array(
        [float(registry[n]["mean_image_angle_deg"]) for n in names])
    vals["photo_median_angle_deg"]=np.array(
        [float(registry[n]["median_image_angle_deg"]) for n in names])
    vals["original_photo_count"]=np.array(
        [int(registry[n]["original_photo_count"]) for n in names])
    vals["BIO12_mm"]=np.array([float(registry[n]["chelsa_bio12"]) for n in names])
    for var in source.EXPOSURES:
        vals[var]=np.array([float(registry[n][var]) for n in names])
    if set(names)-set(registry) or len(set(names))!=N_TREE:
        raise ValueError("No exact original S1 taxon to actual GenBank tip match")
    return names,age,M,dist,vals,float(med)


def kernel(age,M,dist,model,alpha=None,nugget=NUGGET):
    if not 0<=nugget<=0.15:
        raise ValueError("Only zero to 15% technical nugget sensitivity")
    n=len(age)
    if model=="WN":return np.eye(n)
    if model=="BM":
        raw=M
    else:
        if alpha is None or not np.isfinite(alpha) or alpha<=0:
            raise ValueError("OU alpha must be positive")
        if model=="OU_stationary":
            raw=np.exp(-alpha*dist)
        elif model=="OU_fixedroot":
            raw=np.exp(-alpha*dist)*(-np.expm1(-2*alpha*M))/(2*alpha)
        else:raise ValueError("Only WN, BM, stationary OU and fixed-root OU")
    scale=float(np.median(np.diag(raw)))
    if not np.isfinite(scale) or scale<=0:
        raise ValueError("Invalid measured covariance from real molecular branches")
    return raw/scale+nugget*np.eye(n)


def ML(y,X,cov,SE=False):
    n=len(y);p=X.shape[1]
    if cov.shape!=(n,n) or X.shape[0]!=n or n-p<5 or np.linalg.matrix_rank(X)<p:
        raise ValueError("Unidentifiable same-tip phylogenetic regression")
    chol=linalg.cholesky(cov,lower=True,check_finite=True)
    x=linalg.solve_triangular(chol,X,lower=True)
    z=linalg.solve_triangular(chol,y,lower=True)
    beta=np.linalg.lstsq(x,z,rcond=None)[0]
    dev=z-x@beta
    rss=float(dev@dev)
    if rss<=1e-15:raise ValueError("No ecological trait or photo measurement residual variation")
    LL=float(-0.5*(n*(1+np.log(2*np.pi)+np.log(rss/n))
                  +2*np.log(np.diag(chol)).sum()))
    out={"logLik":LL,"mean_regression_beta":beta.tolist()}
    if SE:
        se=np.sqrt(np.diag(np.linalg.inv(x.T@x))*rss/(n-p))
        crit=stats.t.ppf(0.975,n-p)
        out["se"]=se.tolist()
        out["95CI"]=np.stack((beta-crit*se,beta+crit*se),axis=1).tolist()
        out["two_sided_t_p"]=(2*stats.t.sf(np.abs(beta/se),n-p)).tolist()
    return out


def AICc(ll,n,k):
    if n<=k+1:raise ValueError("Too few independent S1 taxon tips for AICc")
    return float(2*k-2*ll+2*k*(k+1)/(n-k-1))


def opt_ou(y,X,age,M,dist,model,nugget):
    # Search both boundary and interior optima; report range hits honestly.
    cache={}
    def score(logalpha):
        if logalpha not in cache:
            a=float(np.exp(logalpha))
            try:fit=ML(y,X,kernel(age,M,dist,model,a,nugget))
            except (ValueError,linalg.LinAlgError):fit={"logLik":-np.inf}
            cache[logalpha]=fit
        return cache[logalpha]["logLik"]
    grid=np.linspace(np.log(.001),np.log(300),53)
    values=np.array([score(float(v)) for v in grid])
    if not np.isfinite(values).any():
        raise ValueError("No identifiable OU source likelihood")
    options=[float(grid[0]),float(grid[-1])]
    for ix in np.argsort(-values)[:4]:
        l=float(grid[max(ix-1,0)]);u=float(grid[min(ix+1,len(grid)-1)])
        if u>l:
            ans=optimize.minimize_scalar(lambda x:-score(float(x)),bounds=(l,u),
                                         method="bounded",options={"xatol":.0005})
            if ans.success:options.append(float(ans.x))
    best=max(options,key=score)
    return {"logLik":score(best),"alpha":float(np.exp(best)),
            "alpha_bound_warning":bool(np.exp(best)<.0015 or np.exp(best)>280),
            "mean_regression_beta":cache[best]["mean_regression_beta"]}


def compare(y,X,age,M,dist,nugget=NUGGET):
    n=len(y);models={}
    for kind in MODES:
        if kind in ("OU_stationary","OU_fixedroot"):
            fit=opt_ou(y,X,age,M,dist,kind,nugget)
        else:fit=ML(y,X,kernel(age,M,dist,kind,nugget=nugget))
        params=X.shape[1]+1+(kind.startswith("OU"))
        models[kind]={"AICc":AICc(fit["logLik"],n,params),
                      "logLik":fit["logLik"],"parameters":params,
                      "alpha":fit.get("alpha"),"alpha_boundary":fit.get("alpha_bound_warning",False)}
    best=min(v["AICc"] for v in models.values())
    winner=min(models,key=lambda z:models[z]["AICc"])
    for m in models.values():m["delta_AICc"]=float(m["AICc"]-best)
    return {"best_AICc_model":winner,"models":models,
            "n_same_original_S1_photo_taxa":n,
            "one_optimum_OU_not_multiregime_ancestry":True,
            "photo_phenotype_not_biological_gravity_state":True}


def lambda_correlation(y,age,M,dist,nugget):
    """Supplemental Pagel correlation strength (different from BM/OU process)."""
    cov=kernel(age,M,dist,"BM",nugget=nugget)
    sd=np.sqrt(np.diag(cov))
    R=cov/(sd[:,None]*sd[None,:])
    np.fill_diagonal(R,1)
    w,V=np.linalg.eigh(R)
    ones=V.T@np.ones(len(y))
    proj=V.T@y
    n=len(y)
    def LL(lam):
        eigen=(1-lam)+lam*w
        if np.min(eigen)<=0:return -np.inf
        inv=1/eigen
        mu=(inv*ones*proj).sum()/(inv*ones*ones).sum()
        rss=float((inv*(proj-mu*ones)**2).sum())
        return float(-.5*(n*(1+np.log(2*np.pi)+np.log(rss/n))
                            +np.log(eigen).sum()))
    opt=optimize.minimize_scalar(lambda z:-LL(float(z)),bounds=(0,1),method="bounded",
                                  options={"xatol":.0005})
    candidates=[0.,1.,float(opt.x)]
    lam=max(candidates,key=LL)
    return {"Pagel_lambda_ML":float(lam),
            "WN_minus_lambda_model_AICc":float(AICc(LL(0),n,2)-AICc(LL(lam),n,3)),
            "likelihood_ratio_over_identity":float(max(0,2*(LL(lam)-LL(0)))),
            "branch_depths_standardized_to_tip_correlation":True,
            "not_BM_or_OU_evolutionary_process":True}


def null_lambda_permutation(y,age,M,dist,nugget,n_perm=999,seed=SEED):
    """Source-tip-label permutation of Pagel λ LRT, conditional on ONE ITS tree."""
    C=kernel(age,M,dist,"BM",nugget=nugget)
    sd=np.sqrt(np.diag(C))
    R=C/(sd[:,None]*sd[None,:])
    np.fill_diagonal(R,1)
    val,U=np.linalg.eigh(R)
    q=U.T@np.ones(len(y))
    n=len(y)
    def statistic(observations):
        transformed=U.T@observations
        def LL(lam):
            d=(1-lam)+lam*val
            if np.min(d)<=0:return -np.inf
            w=1/d
            mu=(w*q*transformed).sum()/(w*q*q).sum()
            ss=float(np.sum(w*(transformed-mu*q)**2))
            return -.5*(n*(1+np.log(2*np.pi)+np.log(ss/n))+np.log(d).sum())
        origin=LL(0)
        opt=optimize.minimize_scalar(lambda x:-LL(float(x)),bounds=(0,1),
                                    method="bounded",options={"xatol":.0005})
        maximum=max(origin,LL(1),-opt.fun)
        return float(max(0,2*(maximum-origin)))
    observed=statistic(y)
    random=np.random.default_rng(seed)
    count=0
    for _ in range(n_perm):
        if statistic(random.permutation(y))>=observed-1e-9:count+=1
    return {"observed_correlation_likelihood_ratio":observed,
            "n_original_taxon_label_permutations":n_perm,
            "permutations_as_or_more_phylogenetic_than_observed":count,
            "permutation_p":float((1+count)/(1+n_perm)),
            "test_is_conditional_on_fixed_ITS_tree_and_original_S1_tip_set":True,
            "not_a_geographic_or_ancestral_rain_selection_control":True}


def climate_regression(vals,age,M,dist,nugget):
    rain=vals["BIO12_mm"]
    angle=vals["photo_mean_angle_deg"]
    five=np.column_stack([
        np.ones(len(rain)),(rain-rain.mean())/100,
        *[(vals[x]-np.mean(vals[x]))/np.std(vals[x])
          for x in ("chelsa_bio01","chelsa_bio04","chelsa_bio15","chelsa_rsds_mean")]
    ])
    modes=compare(angle,five,age,M,dist,nugget)
    out={}
    for kind in MODES:
        a=modes["models"][kind]["alpha"]
        m=ML(angle,five,kernel(age,M,dist,kind,a,nugget),SE=True)
        out[kind]={"BIO12_angle_deg_per_100mm":m["mean_regression_beta"][1],
                   "BIO12_95CI":m["95CI"][1],
                   "BIO12_two_sided_p_conditional_on_covariance":m["two_sided_t_p"][1]}
    return {"covariate_labels":["Intercept","BIO12/100mm","BIO1","BIO4","BIO15","RSDS"],
            "same_109_taxa_original_five_S1_climate_axes":True,
            "angle_versus_current_precipitation_slopes":out,
            "best_residual_covariance_AICc":modes["best_AICc_model"],
            "not_original_observation_level_within_between_mixed_model":True}


def run(s1,accessions,tree_file,discovery,alignment,*,n_perm=999):
    master,_,tree,_,_=source.load_evidence(
        s1,accessions,tree_file,discovery,alignment)
    if len(master)!=N_S1:raise ValueError("Original S1 photo source has unexpected size")
    master_by_original_source={row["original_source_name"]:row for row in master}
    if len(master_by_original_source)!=N_S1:
        raise ValueError("Original S1 species-name entries not unique")
    results={}
    for rooting in ("midpoint","source_Newick"):
        names,age,M,D,vals,med=prepare_data(tree,master_by_original_source,rooting)
        noises=(("tip_nugget5percent",.05),("tip_nugget_zero",0.0)) if rooting=="midpoint" else (("tip_nugget5percent",.05),)
        for label,noise in noises:
            outcome={}
            for phen in TRAITS:
                y=vals[phen]
                X=np.ones((len(y),1))
                comp=compare(y,X,age,M,D,noise)
                comp["Pagel_lambda_correlated_tip_diagnostic"]=lambda_correlation(y,age,M,D,noise)
                if rooting=="midpoint" and noise==.05:
                    comp["Pagel_lambda_permutation"]=null_lambda_permutation(
                        y,age,M,D,noise,n_perm=n_perm,
                        seed=SEED+TRAITS.index(phen))
                outcome[phen]=comp
            results[rooting+"__"+label]={
                "root_scenario":rooting,
                "root_not_verified_by_outgroup":True,
                "original_ITS_median_root_tip_branch_distance_substitutions_per_site":med,
                "n_original_s1_genbank_its_photo_taxa":len(names),
                "n_original_photo_angle_measurements_in_matched_ITS_taxa":int(vals["original_photo_count"].sum()),
                "tip_noise_assumption":noise,
                "evolutionary_model_selection_per_trait":outcome,
                "photo_angle_vs_annual_rain_5covariates":climate_regression(vals,age,M,D,noise),
            }
    return {
        "schema":"azami_T2_true_original109_genbank_ITS_BM_OU_own_S1_photo_niche_v1",
        "source":[
            {"original_S1_251_photo_units_SHA256":source.S1_ALL251_CSV_SHA256,
             "original_109_GenBank_ML_tree_SHA256":source.ITS_ML_TREE_SHA256,
             "original_GenBank_ML_run":39914024548}],
        "original_S1_251_photo_name_units":251,
        "original_S1_WCVP_accepted_keys":244,
        "new_true_ITS_sequence_tree_species":109,
        "source_original_photo_angle_observation_universe":170066,
        "traits_exclusively":["log(BIO12)", "S1 original photo mean angle","S1 original photo median angle"],
        "root_and_tip_noise_model_sensitivity_reported":True,
        "species_evolutionary_tree_verified_from_single_copy_multilocus_data":False,
        "owner_photo_measurement_is_image_vertical_not_gravity":True,
        "outside_flora_botanical_U_D_not_used":True,
        "claimed_stabilizing_selection_or_evolution_of_adaptive_nodding":False,
        "results":results}


def render_AICc_figure(output,data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    group=data["results"]["midpoint__tip_nugget5percent"]["evolutionary_model_selection_per_trait"]
    taxa=["log BIO12 niche","Photo angle mean","Photo angle median"]
    names=["log_BIO12","photo_mean_angle_deg","photo_median_angle_deg"]
    fig,ax=plt.subplots(figsize=(10,5))
    x=np.arange(3)
    for i,model in enumerate(MODES):
        values=[group[k]["models"][model]["delta_AICc"] for k in names]
        ax.bar(x+(i-1.5)*.19,values,width=.18,label=model)
    ax.set_xticks(x,taxa)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("ΔAICc (relative to best model of this trait)")
    ax.legend(ncol=2,frameon=False)
    ax.set_title("S1 original 109 GenBank ITS tree: niche vs photographed head orientation")
    fig.text(.10,.018,"OU stationary-root, fixed-root OU, BM and WN; assumed 5% tip noise, ITS gene tree only.",fontsize=9)
    fig.tight_layout(rect=[0,.045,1,1])
    fig.savefig(output,dpi=160)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--S1-source",type=Path,required=True)
    ap.add_argument("--accessions",type=Path,required=True)
    ap.add_argument("--ITS-tree",type=Path,required=True)
    ap.add_argument("--discovery",type=Path,required=True)
    ap.add_argument("--alignment",type=Path,required=True)
    ap.add_argument("--out-dir",type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=999)
    a=ap.parse_args()
    if a.permutations<99 or a.permutations>9999:
        raise ValueError("No underpowered, unbounded or nonreproducible permutation request")
    results=run(a.S1_source,a.accessions,a.ITS_tree,a.discovery,a.alignment,
                n_perm=a.permutations)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    (a.out_dir/"S1_ORIGINAL_109_ITS_BM_OU_NICHE_PHOTO_RESULT.json").write_text(
        json.dumps(results,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    render_AICc_figure(a.out_dir/"S1_109_ITS_BM_OU_NICHE_VS_PHOTO_AICc.png",results)
    for label,scenario in results["results"].items():
        print("S1_GENBANK_ITS_BM_OU_OWN_PHOTO_ANGLE="+json.dumps({
            "root_and_noise":label,
            "best_models":{key:value["best_AICc_model"]
                for key,value in scenario["evolutionary_model_selection_per_trait"].items()},
            "BIO12_phylogenetic_lambda":scenario["evolutionary_model_selection_per_trait"][
                "log_BIO12"]["Pagel_lambda_correlated_tip_diagnostic"]["Pagel_lambda_ML"],
            "photo_angle_vs_BIO12":{name:val["BIO12_angle_deg_per_100mm"] for name,val in
                scenario["photo_angle_vs_annual_rain_5covariates"][
                "angle_versus_current_precipitation_slopes"].items()},
        },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":main()

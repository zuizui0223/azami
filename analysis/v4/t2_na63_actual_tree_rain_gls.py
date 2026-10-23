"""Direct phylogeny-first NA63 capitulum U/D vs BIO12 on published nuclear tree.

No photographs, camera angle, source-photo thresholds, 52 grafted trees,
invented global phylogeny or posthoc species selection enter this analysis.

One pinned Siniscalchi et al. (2023) published 63-tip nuclear chronology,
63 author-scored HeadPosition states and 63 author-derived BIO12/BIO1/Elevation
records are cross-checked through their immutable source Git blob IDs.

Fit trait-climate association on SAME 63 taxa under tip-independent errors,
strict Brownian shared ancestry and estimated Pagel lambda covariance. Also
show original true phylogenetic clades, sister-clade precipitation comparisons,
and conditional Mk internal-state uncertainty (not identified true gains).
Neither tip-only PGLS nor Mk ancestral probabilities can demonstrate
precipitation-caused adaptive nodding or historical precipitation at origin.
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
from collections import Counter

import numpy as np
from scipy import linalg, optimize, stats

from analysis.v4.t2_pgls_placement_screen import (
    SOURCE_PATH, SOURCE_GIT_BLOB, AUTHOR_NUCLEAR_TREE_PATH,
    AUTHOR_NUCLEAR_TREE_GIT_BLOB, read_author_states, read_author_nuclear_tree,
    parsimony, git_blob_sha, SOURCE_REPO, SOURCE_COMMIT, SOURCE_DOI
)
from analysis.v4.t2_phylogenetic_climate_null import (
    ENV_PATH, ENV_GIT_BLOB, fetch_author, read_published_climate,
    fitted_rate_given_root
)
from analysis.v4.t2_na63_sister_clade_rain import sister_clade_contrasts
from analysis.v4.rain_history import pruning_likelihood, transition_matrix

NA63 = 63
N_NODDING = 5
MIN_PAR = 4
ORDERED_MODELS = ("no_phylogeny", "Brownian_complete_phylogeny", "estimated_lambda")
FLAG_STOPS = {
    "individual_origin_events_identified": False,
    "ancestral_climate_observed": False,
    "precipitation_dependent_U_to_D_rate_fitted": False,
    "full_251_source_taxa_in_one_validated_molecular_phylogeny": False,
    "photos_or_exif_angles_used": False,
    "tree_posteriors_consider_sequence_topology_uncertainty": False,
}


def brownian_shared_ancestry(tree, ordered_tips):
    """Common-root branch-length covariance from the ACTUAL original phylogeny."""
    n = len(ordered_tips)
    if len(set(ordered_tips)) != n or set(ordered_tips) != {
        t.name for t in tree.get_terminals()
    }:
        raise ValueError("Actual nuclear root/tree/taxa not an identical source cohort")
    tips = {t.name:t for t in tree.get_terminals()}
    depths = tree.depths()
    mat = np.empty((n,n), dtype=float)
    for i,a in enumerate(ordered_tips):
        va = float(depths[tips[a]])
        if not np.isfinite(va) or va <= 0:
            raise ValueError("Missing terminal nuclear branch scale")
        mat[i,i] = va
        for j in range(i):
            b = ordered_tips[j]
            v = float(depths[tree.common_ancestor(tips[a], tips[b])])
            mat[i,j] = mat[j,i] = v
    if (np.any(np.diag(mat)<=0) or not np.allclose(mat,mat.T)
            or np.min(np.linalg.eigvalsh(mat)) <= 0):
        raise ValueError("Actual nuclear chronology is not a positive definite tip covariance")
    return mat


def gls(y, design, brownian, lam):
    """Gaussian phylogenetic GLS for PRESENT-DAY BIO12 (mm), not a transition model."""
    if not 0 <= lam <= 1 or not np.all(np.isfinite(y)):
        raise ValueError("Pagel lambda restricted to [0,1] on true source tree")
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Collinear biological trait, temperature, elevation model")
    n, k = design.shape
    if n != len(y) or brownian.shape != (n,n) or n-k < 5:
        raise ValueError("Insufficient original independent nuclear tip coverage")
    cov = lam * brownian + (1-lam)*np.diag(np.diag(brownian))
    ch = linalg.cholesky(cov, lower=True, check_finite=True)
    yw = linalg.solve_triangular(ch, y, lower=True, check_finite=True)
    xw = linalg.solve_triangular(ch, design, lower=True, check_finite=True)
    beta, *_ = np.linalg.lstsq(xw, yw, rcond=None)
    residual = yw - xw @ beta
    rss = float(residual @ residual)
    if rss <= 0:
        raise ValueError("No usable residual variation")
    cov_beta = np.linalg.inv(xw.T @ xw) * (rss/(n-k))
    se = np.sqrt(np.diag(cov_beta))
    tval = beta[1] / se[1]
    ll = -0.5 * (n*math.log(2*math.pi) + 2*np.log(np.diag(ch)).sum()
                 + n*math.log(rss/n) + n)
    tcrit = stats.t.ppf(.975, n-k)
    return {
        "pagel_lambda":float(lam),
        "beta_nodding_BIO12_mm":float(beta[1]),
        "beta_nodding_se_mm":float(se[1]),
        "beta_nodding_ci95_mm":[float(beta[1]-tcrit*se[1]),
                                float(beta[1]+tcrit*se[1])],
        "beta_nodding_one_sided_positive_t_p":float(stats.t.sf(tval,n-k)),
        "beta_nodding_two_sided_t_p":float(2*stats.t.sf(abs(tval),n-k)),
        "root_mean_square_error_cholesky":float(math.sqrt(rss/n)),
        "loglikelihood":float(ll),
        "n_original_nuclear_tips":n,
        "degrees_of_freedom":n-k,
        "linear_model_uses_phenotype_as_categorical_predictor":True,
        "one_lambda_uncertainty_not_propagated_in_Wald_se":True,
    }


def equal_variance_independent_OLS(y, design):
    """Truly unrelated species with IID residual covariance, not Pagel lambda=0.

    On an NON-ultrametric molecular tree, lambda=0 keeps unequal tip
    variances inherited from the branch diagonal and is NOT ordinary OLS.
    This explicit identity-covariance control isolates that difference.
    """
    y=np.asarray(y,dtype=float)
    design=np.asarray(design,dtype=float)
    model=gls(y,design,np.eye(len(y)),0.0)
    model["correlation_structure"]="identity_independent_equal_variance"
    model["no_nuclear_branch_length_information_in_covariance"]=True
    model["Pagel_lambda_0_on_nonultrametric_tree_is_not_IID"]=True
    return model


def fit_models(y, matrix, true_cov):
    models={}
    for label,lam in (("no_phylogeny",0.),("Brownian_complete_phylogeny",1.)):
        models[label]=gls(y,matrix,true_cov,lam)
    # Profile likelihood of the same biological estimand and tip cohort.
    objective=lambda lam: -gls(y,matrix,true_cov,float(lam))["loglikelihood"]
    optim=optimize.minimize_scalar(objective,bounds=(0,1),method="bounded",
                                  options={"xatol":1e-5})
    if not optim.success or not np.isfinite(optim.fun):
        raise ValueError("Phylogenetic covariance lambda ML failed")
    locations=[0.,float(optim.x),1.]
    best=min(locations,key=objective)
    models["estimated_lambda"]=gls(y,matrix,true_cov,best)
    models["estimated_lambda"]["lambda_model_AIC"]=float(
        2*(matrix.shape[1]+2)-2*models["estimated_lambda"]["loglikelihood"])
    models["no_phylogeny"]["lambda_model_AIC"]=float(
        2*(matrix.shape[1]+1)-2*models["no_phylogeny"]["loglikelihood"])
    models["Brownian_complete_phylogeny"]["lambda_model_AIC"]=float(
        2*(matrix.shape[1]+1)-2*models["Brownian_complete_phylogeny"]["loglikelihood"])
    return models


def internal_marginals(tree, observations, *,root_prior="stationary"):
    """Node/edge posterior under symmetric Mk; edge change is NOT event count."""
    if root_prior not in ("stationary","erect"):
        raise ValueError("Unsupported source phylogeny root sensitivity")
    rate, lnL, boundary = fitted_rate_given_root(tree, observations,root_prior)
    _, partial = pruning_likelihood(tree, observations,rate)
    prior=np.array([1.,0.] if root_prior=="erect" else [.5,.5])
    root_vec=prior * partial[tree.root]
    posterior={tree.root:root_vec/root_vec.sum()}
    edges=[]
    parent={c:p for p in tree.find_clades() for c in p.clades}
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:continue
        par=parent[clade]
        prob=transition_matrix(rate,clade.branch_length)
        joint=posterior[par][:,None]*prob*partial[clade][None,:]
        joint/= (prob @ partial[clade])[:,None]
        posterior[clade]=joint.sum(axis=0)
        if not np.isclose(posterior[clade].sum(),1.,atol=1e-8):
            raise ValueError("Internal ancestral-state posteriors not normalized")
        descend=sorted(x.name for x in clade.get_terminals())
        edges.append({
            "clade_tip_count":len(descend),
            "descendant_nodding_tip_count":sum(observations[t] for t in descend),
            "clade_member_names":descend if len(descend)<=6 else descend[:3],
            "original_branch_length":float(clade.branch_length),
            "parent_erect_child_nodding_endpoint_probability":float(joint[0,1]),
            "child_nodding_posterior":float(posterior[clade][1]),
            "historical_rain_at_node_observed":False,
            "independent_U_to_D_gain_count_on_branch_confirmed":False,
        })
    return {
        "root_prior":root_prior,"fitted_symmetric_Mk_rate":rate,
        "Mk_likelihood_at_rate":lnL,
        "rate_optimizer_boundary_warning":boundary,
        "root_nodding_posterior":float(posterior[tree.root][1]),
        "sum_expected_erect_to_nodding_endpoint_changes":
            float(sum(e["parent_erect_child_nodding_endpoint_probability"] for e in edges)),
        "top_edge_endpoint_U_to_D_probability_candidates":
            sorted(edges,key=lambda e:e["parent_erect_child_nodding_endpoint_probability"],reverse=True)[:12],
        "branch_node_posterior_does_not_count_multiple_unseen_flip_events":True,
    }


def analyze(tree,states,climate, *,original_shas=None):
    tip_names=[t.name for t in tree.get_terminals()]
    if (len(tip_names)!=NA63 or len(set(tip_names))!=NA63
        or set(tip_names)!=set(states) or set(tip_names)!=set(climate)
        or sum(states.values())!=N_NODDING
        or any(v not in (0,1) for v in states.values())):
        raise ValueError("Original 63 nuclear tree / 5 genuine botanical posture / climate set changed")
    if len(tree.root.clades)!=2:
        raise ValueError("Original author nuclear Newick has no source-bifurcated root")
    # The published Newick has an authored root split even when Bio.Phylo
    # conservatively marks its Newick rooted flag false. The root itself is a
    # working phylogenetic hypothesis and is tested in two Mk prior modes.
    true_cov=brownian_shared_ancestry(tree,tip_names)
    y=np.array([climate[t]["BIO12"] for t in tip_names],float)
    D=np.array([states[t] for t in tip_names],float)
    temp=np.array([climate[t]["BIO1"] for t in tip_names],float)
    elev=np.array([climate[t]["Elevation"] for t in tip_names],float)
    # Exact same 63 taxa in unadjusted and adjusted model. BIO1/elevation are
    # standardized only for numerical stability; D remains 0 or 1 in mm.
    scaled=np.column_stack((temp,elev))
    sd=scaled.std(axis=0)
    if np.any(sd<=0):
        raise ValueError("No variation in temperature/elevation controls")
    standardized=(scaled-scaled.mean(axis=0))/sd
    basic=np.column_stack((np.ones(len(D)),D))
    full=np.column_stack((basic,standardized))
    subset = sister_clade_contrasts(tree,states,climate)
    if subset["n_maximal_pure_nodding_tip_clades"]!=MIN_PAR or subset["n_nodding_clades_wetter_than_immediate_erect_sister"]!=2:
        raise ValueError("Published NA63 true clade comparisons drifted")
    desc={
        "n_tips":NA63,"D":N_NODDING,"U":NA63-N_NODDING,
        "n_terminal_pure_nodding_clades":subset["n_maximal_pure_nodding_tip_clades"],
        "mean_BIO12_D_mm":float(y[D==1].mean()),
        "mean_BIO12_U_mm":float(y[D==0].mean()),
        "mean_BIO12_D_minus_U_mm":float(y[D==1].mean()-y[D==0].mean()),
        "rain_mm_min_max":[float(y.min()),float(y.max())],
        "phylogenetic_tip_covariance_eigen_min":float(np.linalg.eigvalsh(true_cov).min()),
    }
    meta=original_shas or {}
    return {
        "schema":"azami_v4_T2_NA63_actual_author_phylo_U_D_BIO12_covariance_Mk_20261023_v1",
        "status":"ORIGINAL63_TRUE_NUCLEAR_PHYLOGENY_DIRECT_RAIN_TRAIT_ESTIMATION_NOT_PHOTOS",
        "source":{"article_DOI":SOURCE_DOI,"author_repo":SOURCE_REPO,
                  "author_git_commit":SOURCE_COMMIT,
                  "author_tree_blob":AUTHOR_NUCLEAR_TREE_GIT_BLOB,
                  "author_original_trait_blob":SOURCE_GIT_BLOB,
                  "author_original_climate_blob":ENV_GIT_BLOB,**meta},
        "original_scoped_taxa":desc,
        "unadjusted_BIO12_as_response_mm":fit_models(y,basic,true_cov),
        "BIO1_and_elevation_adjusted_BIO12_as_response_mm":fit_models(y,full,true_cov),
        "IID_species_effect_same_taxa_as_PGLS":{
            "BIO12_only":equal_variance_independent_OLS(y,basic),
            "BIO1_and_elevation_adjusted":equal_variance_independent_OLS(y,full),
        },
        "author_phylogenetic_sister_clades":subset["comparisons"],
        "source_root_Fitch_Sankoff_change_bounds":{
            **parsimony(tree,states),
            "confirmed_independent_U_to_D_gains":False},
        "symmetric_Mk_ancestral_edge_state_posteriors":{
            r:internal_marginals(tree,states,root_prior=r)
            for r in ("erect","stationary")},
        "interpretation":{
            "positive_author_descendant_tip_rain_not_equal_adaptive_event":True,
            "same_species_and_same_climate_naive_vs_true_phylogeny":True,
            "original_nuclear_branch_lengths_used_in_gls":True,
            "GLS_Gaussian_Brownian_Covariance_assumed":True,
            "phylogenetic_effect_estimated_not_years_or_causal_rate":True,
            "root_state_sensitive":True,
            "only_five_author_nodding_tips":True,
            "north_america_only_not_global_251":True,
            "taxon_species_climate_spatial_bias_not_fixed":True,
            **FLAG_STOPS,
        }
    }


def render_figure(tree,states,climate,out_path):
    """Actual published chronological nuclear phylogram plus original BIO12 tips."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    depth=tree.depths()
    terminals=tree.get_terminals()
    n=len(terminals)
    ypos={t:float(i) for i,t in enumerate(terminals)}
    for node in tree.find_clades(order="postorder"):
        if not node.is_terminal():
            ypos[node]=sum(ypos[c] for c in node.clades)/len(node.clades)
    maxdepth=max(float(depth[t]) for t in terminals)
    maxrain=max(float(climate[t.name]["BIO12"]) for t in terminals)
    bar_x=maxdepth+4.0
    bar_width=3.1
    fig,ax=plt.subplots(figsize=(16,20),dpi=175)
    for node in tree.find_clades(order="preorder"):
        x=float(depth[node])
        if not node.is_terminal():
            yy=[ypos[c] for c in node.clades]
            ax.plot([x,x],[min(yy),max(yy)],c="#87929b",linewidth=.65,zorder=1)
        for child in node.clades:
            ax.plot([x,float(depth[child])],[ypos[child],ypos[child]],
                    c="#87929b",linewidth=.65,zorder=1)
    for node in terminals:
        row=ypos[node]
        rain=float(climate[node.name]["BIO12"])
        down=states[node.name]==1
        label=node.name.replace("Cirsium_","C. ").replace("_"," ")
        ax.scatter([float(depth[node])],[row],marker="o",s=17 if down else 8,
                   color=("#bf5b40" if down else "#78939f"),zorder=4)
        ax.text(maxdepth+.20,row,label,fontweight="bold" if down else "normal",
                fontsize=7.3,va="center",color="#9b3d2a" if down else "#293540")
        ax.barh(row,bar_width*rain/maxrain,left=bar_x,height=.60,
                color=("#d47858" if down else "#b8c7cf"),edgecolor="none",zorder=3)
        ax.text(bar_x+bar_width+.12,row,f"{rain:,.0f}",fontsize=6.7,va="center")
    ax.axvline(bar_x-.16,color="#dfdfdf",lw=.7)
    ax.set_xlim(-.20,bar_x+bar_width+1.20)
    ax.set_ylim(n+.5,-2.2)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_frame_on(False)
    ax.text(0,-1.5,"True published 63-tip nuclear phylogeny (Siniscalchi et al. 2023)",
            fontsize=15,fontweight="bold")
    ax.text(0,-.50,"Author botanical head positions; tree branch lengths from source chronology",
            fontsize=9)
    ax.text(bar_x,-.55,"BIO12 • mean annual rain (mm)",fontsize=10,fontweight="bold")
    ax.legend(handles=[
        Line2D([0],[0],marker="o",linestyle="none",color="#bf5b40",label="Nodding (D), 5 tips"),
        Line2D([0],[0],marker="o",linestyle="none",color="#78939f",label="Erect (U), 58 tips")
    ],loc="lower left",fontsize=9,frameon=False,bbox_to_anchor=(0,0))
    fig.text(.1,.018,"Photographic EXIF angle not used. Tree-adjacent tip rain is modern niche, not ancestral rainfall.",
             fontsize=9)
    out_path=Path(out_path)
    out_path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out_path,dpi=175,bbox_inches="tight")
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--no-figure",action="store_true")
    a=p.parse_args()
    trait_data, trait_sha=fetch_author(SOURCE_PATH,SOURCE_GIT_BLOB)
    climate_data, climate_sha=fetch_author(ENV_PATH,ENV_GIT_BLOB)
    tree_data, tree_sha=fetch_author(AUTHOR_NUCLEAR_TREE_PATH,
                                    AUTHOR_NUCLEAR_TREE_GIT_BLOB)
    states=read_author_states(trait_data)
    climate=read_published_climate(climate_data)
    tree=read_author_nuclear_tree(tree_data,states)
    actual=analyze(tree,states,climate,original_shas={
        "traits_sha256":trait_sha,"climate_sha256":climate_sha,
        "nuclear_Newick_sha256":tree_sha})
    a.out_dir.mkdir(parents=True,exist_ok=True)
    out=a.out_dir/"T2_NA63_REAL_NUCLEAR_PHYLOGENETIC_RAIN_ASSOCIATION.json"
    out.write_text(json.dumps(actual,indent=2,ensure_ascii=False,
                              sort_keys=True,allow_nan=False)+"\n",
                   encoding="utf-8")
    with (a.out_dir/"T2_NA63_REAL_TREE_TIP_HEAD_POSITION_BIO12_BIO1_ELEVATION.csv").open(
            "w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=[
            "original_tree_tip","botanical_nodding_D1_U0","modern_BIO12_mm",
            "original_BIO1_units","elevation_m","source_tree_phylogeny_is_real"])
        w.writeheader()
        for t in tree.get_terminals():
            w.writerow({"original_tree_tip":t.name,
                        "botanical_nodding_D1_U0":states[t.name],
                        "modern_BIO12_mm":climate[t.name]["BIO12"],
                        "original_BIO1_units":climate[t.name]["BIO1"],
                        "elevation_m":climate[t.name]["Elevation"],
                        "source_tree_phylogeny_is_real":True})
    if not a.no_figure:
        render_figure(tree,states,climate,
            a.out_dir/"T2_NA63_REAL_NUCLEAR_TREE_POSTURE_RAIN.png")
    print("T2_NA63_REAL_TREE_GLS_RAIN="+json.dumps({
        "n":actual["original_scoped_taxa"]["n_tips"],
        "phylogenetic_BM":actual["BIO1_and_elevation_adjusted_BIO12_as_response_mm"]["Brownian_complete_phylogeny"],
        "estimated_lambda":actual["BIO1_and_elevation_adjusted_BIO12_as_response_mm"]["estimated_lambda"],
        "unrelated":actual["BIO1_and_elevation_adjusted_BIO12_as_response_mm"]["no_phylogeny"],
        "sister_clades":len(actual["author_phylogenetic_sister_clades"]),
        "mk_erect_root":actual["symmetric_Mk_ancestral_edge_state_posteriors"]["erect"]["root_nodding_posterior"]
    },ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

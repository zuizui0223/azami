"""T2 observational rain-vs-temperature test using real published head states.

The predeclared target is the signed standardized BIO12 contrast of nodding
versus erect tips after adjusting species' BIO1 and elevation. Two nulls:
fixed-prevalence label permutations and conditional-prevalence symmetric Mk
tip simulations on the 2023 nuclear tree or 52 historical PGLS placements.
This is NOT a likelihood test for a rain-dependent evolutionary rate.
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import urllib.request

from Bio import Phylo
import numpy as np
from scipy.optimize import minimize_scalar

from analysis.v4.rain_history import pruning_likelihood, rate_fit
from analysis.v4.t2_pgls_placement_screen import (
    SOURCE_REPO, SOURCE_COMMIT, SOURCE_DOI, SOURCE_GIT_BLOB,
    SOURCE_PATH, AUTHOR_NUCLEAR_TREE_PATH, AUTHOR_NUCLEAR_TREE_GIT_BLOB,
    git_blob_sha, read_author_states, read_author_nuclear_tree,
    check_placement_tree,
)

ENV_PATH = "environmental/cirsium_environmental_allcombined_DR.csv"
ENV_GIT_BLOB = "6ca878423881db13f2af4624380bdce342afdf1b"
PREDECLARED_TARGETS = ("BIO12", "BIO1")
OTHER_CONTROL = "Elevation"


def fetch_author(path: str, expected_blob: str) -> tuple[bytes, str]:
    url = f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{path}"
    with urllib.request.urlopen(url, timeout=90) as response:
        raw = response.read()
    if git_blob_sha(raw) != expected_blob:
        raise ValueError(f"{path}: author source Git blob mismatch")
    return raw, hashlib.sha256(raw).hexdigest()


def read_published_climate(raw: bytes) -> dict[str, dict[str, float]]:
    if git_blob_sha(raw) != ENV_GIT_BLOB:
        raise ValueError("Published environmental data Git blob mismatch")
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(rows) != 63:
        raise ValueError("Expected 63 original ecological rows")
    climate = {}
    for row in rows:
        name = row.get("taxa", "")
        if (not name.startswith("Cirsium_") or name in climate
                or row.get("Species", "") != name):
            raise ValueError("Duplicate/unaligned published taxon record")
        try:
            vals = {key: float(row[key]) for key in ("BIO12", "BIO1", "Elevation")}
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError(f"Incomplete published BIO12/BIO1/elevation: {name}") from exc
        if not all(np.isfinite(v) for v in vals.values()) or vals["BIO12"] <= 0:
            raise ValueError("Invalid published climatic values")
        climate[name] = vals
    return climate


def climate_contrast(names: list[str], states, climate):
    """Trait-blind linear residualisation on the identical complete tip cohort."""
    if (len(names) != len(set(names)) or not names
            or not set(names) <= states.keys() or not set(names) <= climate.keys()):
        raise ValueError("Invalid common original climate/tip cohort")
    y = np.array([states[t] for t in names], dtype=np.int8)
    k = int(y.sum())
    if k < 2 or len(y) - k < 2:
        raise ValueError("Both states need at least 2 observed taxa")
    values = np.array([[climate[t][var] for var in ("BIO12", "BIO1", "Elevation")]
                       for t in names], dtype=float)
    sd = np.std(values, axis=0)
    if np.any(sd <= 1e-12):
        raise ValueError("Environmental predictor has no variance")
    z = (values - values.mean(axis=0)) / sd
    if np.linalg.matrix_rank(np.column_stack((np.ones(len(y)),z))) < 4:
        raise ValueError("Climate/elevation design is collinear")
    adjusted = []
    for j in (0, 1):
        other = 1 - j
        cov = np.column_stack((np.ones(len(y)),z[:,other],z[:,2]))
        resid = z[:,j] - cov @ np.linalg.lstsq(cov, z[:,j], rcond=None)[0]
        residual_sd = np.std(resid)
        if residual_sd < 1e-10:
            raise ValueError("Residualized climate axis has no variation")
        adjusted.append(resid / residual_sd)
    r = np.column_stack(adjusted)
    observed = signed_contrast(y[None,:],r)[0]
    return y, r, {
        "n_tips": len(names), "n_nodding": k, "n_erect": len(y)-k,
        "measured_original_BIO12_units": "published source values, annual precipitation typically mm",
        "measured_original_BIO1_units": "published source values, standardised not interpreted as degrees C",
        "primary_estimand": "nodding-minus-erect contrast in standardized rainfall residual after BIO1 and elevation",
        "secondary_estimand": "nodding-minus-erect contrast in standardized temperature residual after BIO12 and elevation",
        "observed_adjusted_rain_contrast_sd_units": float(observed[0]),
        "observed_adjusted_temperature_contrast_sd_units": float(observed[1]),
        "raw_nodding_BIO12_mean": float(values[y==1,0].mean()),
        "raw_erect_BIO12_mean": float(values[y==0,0].mean()),
        "raw_nodding_BIO1_source_mean": float(values[y==1,1].mean()),
        "raw_erect_BIO1_source_mean": float(values[y==0,1].mean()),
        "nodding_taxa": [t for t,s in zip(names,y) if s==1],
    }


def signed_contrast(y_batch, environmental_resid):
    y = np.asarray(y_batch,dtype=float)
    if y.ndim != 2 or environmental_resid.shape!=(y.shape[1],2):
        raise ValueError("Expected [draws,tips] and [tips,two axes]")
    n=y.shape[1]
    count=y.sum(axis=1)
    if np.any(count < 1) or np.any(count >= n):
        raise ValueError("Invalid phenotype prevalence")
    # Environmental residuals have zero-mean for each complete tip cohort.
    return (y @ environmental_resid) * (n / (count * (n-count)))[:,None]


def simulate_mk_tip_matrix(tree, rate, n, rng, root_mode="stationary"):
    """Simulate independent CTMC process histories with one shared root per draw."""
    if rate <= 0 or not np.isfinite(rate):
        raise ValueError("Invalid Mk rate")
    if root_mode not in {"stationary","erect"}:
        raise ValueError("Unknown rooted null hypothesis")
    if n < 1:
        raise ValueError("At least one simulation required")
    nodes = list(tree.find_clades(order="preorder"))
    state = {tree.root: (np.zeros(n,dtype=np.int8) if root_mode=="erect"
                         else rng.integers(0,2,size=n,dtype=np.int8))}
    parent = {child: node for node in nodes for child in node.clades}
    for node in nodes[1:]:
        length = float(node.branch_length)
        if not np.isfinite(length) or length < 0:
            raise ValueError("Missing/nonfinite/negative evolutionary branches")
        p = -np.expm1(-2*rate*length) / 2.0
        flip = rng.random(n) < p
        state[node] = np.bitwise_xor(state[parent[node]],flip.astype(np.int8))
    return np.column_stack([state[t] for t in tree.get_terminals()])


def fitted_rate_given_root(tree, state_dict, mode):
    if mode == "stationary":
        rate, lnL, boundary = rate_fit(tree,state_dict)
        return rate, lnL, boundary
    if mode != "erect":
        raise ValueError("Only stationary or erect root allowed")
    lengths=[c.branch_length for c in tree.find_clades() if c is not tree.root]
    scale = float(np.median(lengths))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid branch scale")
    def obj(log_rate):
        q=np.exp(log_rate)/scale
        ll, partials=pruning_likelihood(tree,state_dict,q)
        root_vec=partials[tree.root]
        likelihood_ratio=float(root_vec[0] / (0.5*sum(root_vec)))
        return -(ll+np.log(likelihood_ratio)) if likelihood_ratio>0 else 1e100
    fit=minimize_scalar(obj,bounds=(-12,8),method="bounded",options={"xatol":1e-5})
    if not fit.success or not np.isfinite(fit.fun):
        raise ValueError("Fixed-erect-root fit failed")
    return float(np.exp(fit.x)/scale),float(-fit.fun),bool(min(abs(fit.x+12),abs(fit.x-8))<0.005)


def simulate_fixed_prevalence_shuffle(n,k,residual,proposals,seed):
    rng = np.random.default_rng(seed)
    results=[]
    completed=0
    while completed < proposals:
        batch=min(4096, proposals-completed)
        u=rng.random((batch,n))
        indices=np.argpartition(u,k-1,axis=1)[:,:k]
        state=np.zeros((batch,n),dtype=np.int8)
        state[np.arange(batch)[:,None],indices]=1
        results.append(signed_contrast(state,residual))
        completed+=batch
    return np.concatenate(results,axis=0)


def conditional_mk_null(tree,states,residual,*,root_mode="erect",seed=1,
                        target_draws=600,max_proposals=100000,batch=4096,
                        min_accepted=200):
    tips=[t.name for t in tree.get_terminals()]
    y=np.array([states[t] for t in tips],dtype=np.int8)
    k=int(y.sum())
    if k<2 or k>len(y)-2:
        raise ValueError("Uninformative binary cohort")
    rate, ll, boundary=fitted_rate_given_root(tree,states,root_mode)
    rng=np.random.default_rng(seed)
    results=[]
    accepted=0
    proposed=0
    while proposed<max_proposals and accepted<target_draws:
        n=min(batch,max_proposals-proposed)
        sims=simulate_mk_tip_matrix(tree,rate,n,rng,root_mode)
        proposed+=n
        matched=sims[sims.sum(axis=1)==k]
        if len(matched):
            results.append(signed_contrast(matched,residual))
            accepted+=len(matched)
    if accepted:
        null=np.concatenate(results,axis=0)[:target_draws]
    else:
        null=np.empty((0,2),dtype=float)
    return null, {
        "null_name":"conditional_binary_symmetric_Mk_taxon_prevalence",
        "null_root":root_mode, "fitted_Mk_rate_tree_units":rate,
        "fitted_Mk_loglikelihood":ll, "optimizer_near_rate_bound":boundary,
        "total_simulated_tip_histories":proposed,
        "accepted_occupancy_matched":accepted,
        "used_accepted_null_draws":len(null),
        "observed_tip_prevalence":k,
        "acceptance_fraction":accepted/proposed,
        "status":("calibrated_conditional_null" if len(null)>=min_accepted
                  else "insufficient_conditioned_simulations")
    }


def inference_from_null(observed,null,diagnostic,min_accepted=200):
    out={**diagnostic}
    if len(null)<min_accepted:
        return {**out,"observed_contrast":observed.tolist(),
                "raw_empirical_two_sided_p":None,
                "two_axis_maxT_adjusted_p":None}
    mean=null.mean(axis=0)
    sd=null.std(axis=0)
    if np.any(sd<1e-10):
        return {**out,"status":"degenerate_phylogenetic_null",
                "observed_contrast":observed.tolist(),
                "raw_empirical_two_sided_p":None,
                "two_axis_maxT_adjusted_p":None}
    observed_z=np.abs((observed-mean)/sd)
    null_z=np.abs((null-mean)/sd)
    p_raw=(1+np.sum(null_z>=observed_z,axis=0))/(1+len(null))
    null_max=null_z.max(axis=1)
    p_adjusted=(1+np.sum(null_max[:,None]>=observed_z,axis=0))/(1+len(null))
    return {**out,
       "observed_contrast":observed.tolist(),
       "null_mean":mean.tolist(),
       "null_sd":sd.tolist(),
       "null_q025_q50_q975":[np.quantile(null[:,j],[0.025,0.5,0.975]).tolist() for j in range(2)],
       "raw_empirical_two_sided_p":p_raw.tolist(),
       "two_axis_maxT_adjusted_p":p_adjusted.tolist(),
       "axis_order":["BIO12_adjusted_for_BIO1_elevation","BIO1_adjusted_for_BIO12_elevation"],
       "p_method":"two-sided null-centred Monte Carlo rank with plus-one correction; maxT controls two registered axes conditional on fitted null"}


def evaluate_tree(tree,states,climate,*,seed,root_modes=("erect","stationary"),
                  target_draws=600,max_proposals=100000,min_accepted=200):
    original={tip.name for tip in tree.get_terminals()}
    matched=sorted(original & set(states) & set(climate))
    if len(matched)<12 or sum(states[t] for t in matched)<2:
        return {"status":"insufficient_published_head_state_taxa","n_matched":len(matched)}
    phy=deepcopy(tree)
    for tip in list(phy.get_terminals()):
        if tip.name not in matched:
            phy.prune(tip)
    observed={t:states[t] for t in matched}
    # Phylo iteration order must be matched to the simulated tip-matrix order.
    tip_names=[t.name for t in phy.get_terminals()]
    y,resid,desc=climate_contrast(tip_names,states,climate)
    obs=signed_contrast(y[None,:],resid)[0]
    result={"status":"assessable_published_head_state_climate",
            "trait_and_environment":desc,
            "root_hypotheses":{}}
    shuffled=simulate_fixed_prevalence_shuffle(len(y),int(y.sum()),resid,
                 proposals=max(2000, target_draws * 4),seed=seed+10000)
    result["nonphylogenetic_fixed_prevalence_shuffle"]=inference_from_null(
       obs,shuffled,{"status":"calibrated_label_shuffling_not_phylogenetic",
                     "used_accepted_null_draws":len(shuffled),
                     "null_name":"nonphylogenetic_fixed_prevalence_tip_label_shuffling"},
       min_accepted=200)
    for i,root_mode in enumerate(root_modes):
        null, diag=conditional_mk_null(phy,observed,resid,root_mode=root_mode,
                                       seed=seed+i,target_draws=target_draws,
                                       max_proposals=max_proposals,
                                       min_accepted=min_accepted)
        result["root_hypotheses"][root_mode]=inference_from_null(obs,null,diag,min_accepted)
    return result


def execute(trees,states,climate,sha,*,draws=600,max_proposals=100000,seed=20261021):
    if len(trees)!=52 or len(states)!=63 or sum(states.values())!=5 or len(climate)!=63:
        raise ValueError("Frozen source-cohort contract changed")
    native=trees[0][1]
    ref_tips={t.name for t in native.get_terminals()}
    if any({t.name for t in t.get_terminals()} != ref_tips for _,t in trees):
        raise ValueError("PGLS trees must have identical taxon coverage")
    full=sha["nuclear_tree"]
    if full is None:
        raise ValueError("Missing independently published nuclear tree")
    results={}
    results["published_nuclear_63"]=evaluate_tree(full,states,climate,seed=seed,
                                                target_draws=draws,
                                                max_proposals=max_proposals)
    # Same scientific support for comparing nuclear and PGLS histories.
    shared=ref_tips & set(states) & set(climate)
    trimmed=deepcopy(full)
    for tip in list(trimmed.get_terminals()):
        if tip.name not in shared:
            trimmed.prune(tip)
    results["published_nuclear_shared_44"]=evaluate_tree(trimmed,states,climate,
                    seed=seed+100, target_draws=draws,max_proposals=max_proposals)
    for i,(name,tree) in enumerate(trees):
        check_placement_tree(tree)
        results[name]=evaluate_tree(tree,states,climate,seed=seed+200+2*i,
                     target_draws=draws,max_proposals=max_proposals)
    summary={
      "status":"EMPIRICAL_NA_AUTHOR_SOURCE_CLIMATE_ASSOCIATION_PHYLOGENETIC_NULL",
      "scientific_scope":"tip-level morphological precipitation association conditional on temperature/elevation and Mk null; no time-dependent rain-selection model",
      "author_source_doi":SOURCE_DOI,
      "source_provenance":{k:v for k,v in sha.items() if k!="nuclear_tree"},
      "registered_axis_order":["BIO12","BIO1"],"controlled_elevation":True,
      "n_placement_scenarios":52,
      "n_published_nuclear_taxa":results["published_nuclear_63"]["trait_and_environment"]["n_tips"],
      "n_shared_nuclear_and_pgls_taxa":results["published_nuclear_shared_44"]["trait_and_environment"]["n_tips"],
      "root_hypotheses":["erect","stationary"],
      "no_inference_of_rain_adaptation_or_transition_dates":True,
      "no_global_clade_generalization":True,
      "source_grouping":"author species climate mean, not flowering-time weather; taxon entries not independent environmental samples",
    }
    pvals={}
    for key,value in results.items():
        if "root_hypotheses" in value:
            for root,hyp in value["root_hypotheses"].items():
                pvals[(key,root)]=hyp.get("two_axis_maxT_adjusted_p")
    for mode in ("erect","stationary"):
        candidates=[v[0] for (tree,rt),v in pvals.items() if tree not in {
              "published_nuclear_63","published_nuclear_shared_44"} and rt==mode
             and v is not None]
        summary[mode+"_placement_calibrated_n"]=len(candidates)
        summary[mode+"_placement_BIO12_adjusted_p_range"]=(
            [float(min(candidates)),float(max(candidates))] if candidates else None)
    summary["interpretation_stop"]="no direct test of time-varying evolutionary rate; 5 pendulous states among 63; prior published diversification study also used rainfall; no external replication"
    return {"summary":summary,"trees":results}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",required=True)
    parser.add_argument("--draws",type=int,default=600)
    parser.add_argument("--max-proposals",type=int,default=100000)
    opts=parser.parse_args()
    if not 100<=opts.draws<=20000 or not 1000<=opts.max_proposals<=1000000:
        raise ValueError("Invalid simulation budget")
    from analysis.v4.native_phylogeny import restore_trees
    trait_raw,traitsha=fetch_author(SOURCE_PATH,SOURCE_GIT_BLOB)
    env_raw,envsha=fetch_author(ENV_PATH,ENV_GIT_BLOB)
    nuclear_raw,nuclearsha=fetch_author(AUTHOR_NUCLEAR_TREE_PATH,AUTHOR_NUCLEAR_TREE_GIT_BLOB)
    states=read_author_states(trait_raw)
    climate=read_published_climate(env_raw)
    nuclear=read_author_nuclear_tree(nuclear_raw,states)
    if set(states)!=set(climate):
        raise ValueError("Climate/taxon crosswalk mismatch")
    trees=restore_trees()
    receipt=execute(trees,states,climate,
          {"head_state_git_blob":SOURCE_GIT_BLOB,
           "head_state_sha256":traitsha,
           "environment_git_blob":ENV_GIT_BLOB,
           "environment_sha256":envsha,
           "nuclear_tree_git_blob":AUTHOR_NUCLEAR_TREE_GIT_BLOB,
           "nuclear_tree_sha256":nuclearsha,
           "nuclear_tree":nuclear},
          draws=opts.draws,max_proposals=opts.max_proposals)
    destination=Path(opts.out); destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("T2_PHYLOGENETIC_CLIMATE_NULL_SUMMARY="+json.dumps({
      **receipt["summary"],
      "nuclear63_results":receipt["trees"]["published_nuclear_63"],
      "nuclear44_results":receipt["trees"]["published_nuclear_shared_44"],
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    main()

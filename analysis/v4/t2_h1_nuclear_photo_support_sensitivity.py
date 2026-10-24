"""T2 H1: source-verified nuclear phylogeny under minimum photograph support.

Re-use EXACT original NA44/Japan12 independent author-scored botanical U/D,
the SAME photographic CHELSA BIO12/BIO1, and the SAME two authenticated
molecular trees as the original H1 execution (39903303479).

The >=5 photos cutoff is an expressly POST-HOC sensitivity (not a new
preregistered primary hypothesis or a way to chase P < 0.05).
The >=10 cutoff is descriptive only when Japan has fewer than 3 upright tips.
No fictitious 62-species phylogeny, no claim of independent adaptive origins.
"""
from __future__ import annotations
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.t2_h1_real_nuclear_phylo_same_s1 import (
    source_panels, make_panel, neutral_label_null, phylo_null
)

SEED = 20261023
PANEL_REQUIRE = {
    1: {"NA44": (44,3,41), "Japan12": (12,5,7)},
    5: {"NA44": (42,3,39), "Japan12": (6,3,3)},
    10: {"NA44": (41,3,38), "Japan12": (5,3,2)},
}
EA20_EXACT_NAMES = {
    "Cirsium tatakaense", "Cirsium morii", "Cirsium yezoense",
    "Cirsium suffultum", "Cirsium kamtschaticum",
    "Cirsium kawakamii", "Cirsium brevicaule", "Cirsium irumtiense",
}
MIN_TIPS_PER_STATE_FOR_MK = 3
SOURCE_ORIGINAL_H1_ACTION = 39903303479


def photo_scope(panel, nphotos):
    if nphotos not in PANEL_REQUIRE:
        raise ValueError("Only preregistered report scopes 1/5/10 photos are supported")
    name = panel["name"]
    if name not in PANEL_REQUIRE[nphotos]:
        raise ValueError("Unexpected original independent published nuclear study")
    original = panel["data"]
    thinned = original[original["original_photo_count"] >= nphotos].reset_index().copy()
    n = len(thinned)
    k = int((thinned["botanical_U_D_state"] == "D").sum())
    u = n-k
    if (n,k,u) != PANEL_REQUIRE[nphotos][name]:
        raise ValueError("Original named S1 data/EAzami/NA63 photo support changed")
    if n < 4 or not u or not k:
        raise ValueError("Both botanical postures must have observed tips")
    d = thinned[thinned["botanical_U_D_state"]=="D"]
    e = thinned[thinned["botanical_U_D_state"]=="U"]
    record={
        "exact_tree_panel": name, "n_original_taxa_with_at_least_nphotos": nphotos,
        "n_matched_source_taxa": n, "D": k, "U": u,
        "mean_BIO12_D_mm":float(d["chelsa_bio12"].mean()),
        "mean_BIO12_U_mm":float(e["chelsa_bio12"].mean()),
        "descriptive_D_minus_U_BIO12_mm":float(d["chelsa_bio12"].mean()-e["chelsa_bio12"].mean()),
        "source_names": sorted(thinned["original_source_name"].tolist()),
        "new_tip_phenotype_from_photo_angle": False,
        "taxon_concept_voucher_final_acceptance_unresolved": True,
        "photo_quality_and_sampling_not_random": True,
        "phylogenetic_Mk_test_status": (
            "eligible_posthoc_source_tip_panel" if min(k,u)>=MIN_TIPS_PER_STATE_FOR_MK
            else "BLOCKED_POSTHOC_TOO_FEW_INDEPENDENT_BOTANICAL_STATES"
        )
    }
    candidate=None
    if min(k,u)>=MIN_TIPS_PER_STATE_FOR_MK:
        candidate=make_panel(name+"_n"+str(nphotos),thinned,panel["tree"])
    return candidate,record


def summarize_original(df,base_panels,*,draws=30000):
    if draws < 2000 or draws > 100000:
        raise ValueError("No inadequately simulated or unlimited Mk null")
    if (len(df)!=251 or df["accepted_key"].nunique()!=244
        or df["original_photo_count"].sum()!=170066):
        raise ValueError("Original 251/244/170066 photograph source evidence changed")
    original_scored=df[df["botanical_U_D_state"].isin(["U","D"])]
    if (len(original_scored)!=66 or
        sum(original_scored["botanical_U_D_state"]=="D")!=12):
        raise ValueError("Original 66 U/D independently scored source changed")
    if len(base_panels)!=2 or set(p["name"] for p in base_panels)!={"NA44","Japan12"}:
        raise ValueError("Must use exactly NA44 and Japan12 genuine independently sourced trees")
    panel_names={p["name"]:set(p["data"]["original_source_name"]) for p in base_panels}
    if panel_names["NA44"]&panel_names["Japan12"]:
        raise ValueError("NA44 and Japan12 names overlap: cannot assume independent lineages")
    if len(panel_names["NA44"]|panel_names["Japan12"])!=56:
        raise ValueError("Source verified H1 56 real-molecular compatible photo taxa changed")
    ea20=df[df["original_source_name"].isin(EA20_EXACT_NAMES)]
    if len(ea20)!=8 or sum(ea20.botanical_U_D_state=="D")!=5:
        raise ValueError("Original EAzami20 exact eight source labels changed")
    east_union=panel_names["Japan12"]|EA20_EXACT_NAMES
    if len(east_union)!=18 or len((panel_names["NA44"]|east_union))!=62:
        raise ValueError("EA20 and Japan overlap count changed; never double-count")
    coverage={}
    all_panels={}
    for photos in (1,5,10):
        names=panel_names["NA44"]|east_union
        sample=df[(df["original_source_name"].isin(names))
                  & (df["original_photo_count"]>=photos)]
        d=sample[sample.botanical_U_D_state=="D"]
        u=sample[sample.botanical_U_D_state=="U"]
        coverage[str(photos)]={
            "S1_exact_named_taxa_in_true_nuclear_tree_sources":len(sample),
            "D":len(d), "U":len(u),
            "nuclear_tree_assembled_for_all_these_taxa":False,
            "union_across_three_source_panels_overlap_not_independent_replications":True,
            "descriptive_rain_D_minus_U_mm":float(d.chelsa_bio12.mean()-u.chelsa_bio12.mean()),
        }
        current={}
        for orig in base_panels:
            p, rec=photo_scope(orig,photos)
            current[orig["name"]]=(p,rec)
        all_panels[str(photos)]={k:record for k,(_p,record) in current.items()}
    # Compute same source tips (not different cohorts) under NA44/JP12
    # independently rooted molecular trees. Only n>=5 is newly fitted.
    true_panels={}
    observed={}
    exact={}
    nulls={}
    for i,orig in enumerate(base_panels):
        p, rec=photo_scope(orig,5)
        if p is None:
            raise ValueError("Five-photo source-state panel failed its original gate")
        true_panels[orig["name"]]=p
        observed[orig["name"]]=p["observed"]
        null=neutral_label_null(p)
        exact[orig["name"]]=null
        rec["observed_within_tree_standardized_BIO12_and_BIO1_residual"]=p["observed"].tolist()
        rec["nonphylo_frozen_prevalence_exact_reassignments"]=len(null)
        rec["nonphylo_H1_one_sided_p"]=(null>=p["observed"]-1e-12).mean(axis=0).tolist()
        rec["phylogeny_Mk_root_sensitivity"]={}
        for j,root in enumerate(("erect","stationary")):
            sims, meta=phylo_null(p,root,SEED+101*i+400*j+5,draws)
            nulls[(orig["name"],root)]=sims
            rec["phylogeny_Mk_root_sensitivity"][root]={
                "n_accepted_source_conditional_null_tip_states":len(sims),
                "fitted_symmetric_Mk_rate_in_own_tree_branch_units":meta["fitted_Mk_rate_tree_units"],
                "null_mean_of_standardized_contrasts":sims.mean(axis=0).tolist(),
                "null_q025_q975":np.quantile(sims,[.025,.975],axis=0).T.tolist(),
                "posthoc_H1_one_sided_p":((sims>=p["observed"]-1e-12).sum(axis=0)+1).__truediv__(draws+1).tolist(),
            }
        all_panels["5"][orig["name"]]=rec
    a,b=true_panels["NA44"],true_panels["Japan12"]
    w=a["weight"]+b["weight"]
    target=(a["weight"]*a["observed"]+b["weight"]*b["observed"])/w
    n0,n1=exact["NA44"],exact["Japan12"]
    exceed=np.zeros(2,dtype=np.int64)
    for start in range(0,len(n0),512):
        x=(a["weight"]*n0[start:start+512,None,:]+b["weight"]*n1[None,:,:])/w
        exceed+=(x>=target-1e-12).sum(axis=(0,1))
    denominator=len(n0)*len(n1)
    calibrated={}
    for first,second in product(("erect","stationary"),repeat=2):
        mix=(a["weight"]*nulls[("NA44",first)]
             +b["weight"]*nulls[("Japan12",second)])/w
        calibrated[first+"__"+second]={
            "n_conditioned_independent_source_tree_scenarios":len(mix),
            "null_mean":mix.mean(axis=0).tolist(),
            "H1_posthoc_one_sided_p":((mix>=target-1e-12).sum(axis=0)+1).__truediv__(draws+1).tolist(),
        }
    # EAzami20 has 8 real named taxa but overlaps Japan12 on 2.
    # Changes in sign on filtering to 5-photos are intentionally
    # DESCRIPTIVE ONLY; one upright 5-photo tip -> Mk cannot be fit.
    east_alternatives={}
    for photos in (1,5,10):
        f=ea20[ea20.original_photo_count>=photos]
        d=f[f.botanical_U_D_state=="D"]
        u=f[f.botanical_U_D_state=="U"]
        east_alternatives[str(photos)]={
            "exact_source_matched_EA20_tips":len(f), "D":len(d), "U":len(u),
            "rain_D_minus_U_mm":float(d.chelsa_bio12.mean()-u.chelsa_bio12.mean())
                                    if len(d) and len(u) else None,
            "valid_standalone_Mk_with_atleast3_per_state":min(len(d),len(u))>=3,
            "independent_eastasia20_phylogenetic_null_newly_refitted_here":False,
        }
    return {
        "schema":"azami_T2_H1_original56_nuclear_photo_min5_phylo_and_original62_overlap_v1",
        "primary_source_real_run":SOURCE_ORIGINAL_H1_ACTION,
        "status":"REAL_NUCLEAR_TWO_TREE_POSTHOC_PHOTO_DEPTH_SENSITIVITY",
        "source_original_251_photo_taxa_retained":251,
        "original_66_independent_botanical_scored":66,
        "original_66_descriptive_BIO12_D_minus_U_mm":
            float(original_scored.loc[original_scored.botanical_U_D_state=="D","chelsa_bio12"].mean()
                  -original_scored.loc[original_scored.botanical_U_D_state=="U","chelsa_bio12"].mean()),
        "two_original_independent_nuclear_trees":list(panel_names),
        "same_phylogeny_base56_then_photo_thinned_states":all_panels,
        "same_taxon_two_tree_min5_joint_H1":{
            "n_matched_original_taxa":48,"D":6,"U":42,
            "observed_standardized_rain_then_temperature_residual":target.tolist(),
            "n_exact_unrelated_tip_assignments":denominator,
            "unrelated_one_sided_exact_p":(exceed/denominator).tolist(),
            "within_each_own_nuclear_tree_conditional_Mk":calibrated,
            "no_fictitious_global_root_or_shared_branch_scaling":True,
        },
        "EA20_separate_unmerged_alternative_photo_depth":east_alternatives,
        "union_of_three_nuclear_source_name_panels_DESCRIPTIVE_ONLY":coverage,
        "missing_in_all_three_source_trees_original_scored_names":sorted(
            set(original_scored.original_source_name)-panel_names["NA44"]-east_union),
        "posthoc_photo_thresholds_not_prior_model_selection":True,
        "photo_centroid_not_full_native_range_or_ancestral_precipitation":True,
        "no_single_62tip_molecular_tree_or_independent_gain_count":True,
        "nuclear_single_copy_specimen_voucher_taxon_concept_review_incomplete":True,
        "source_66_nuclear_56_match_change_not_exclusively_phylogenetic_adjustment":True,
        "historical_rain_dependent_U_to_D_transitions_tested":False,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original251",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--draws",type=int,default=30000)
    args=ap.parse_args()
    base=source_panels(args.original251)
    original=pd.read_csv(args.original251)
    out=summarize_original(original,base,draws=args.draws)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    c=out["same_taxon_two_tree_min5_joint_H1"]
    print("T2_H1_SAME_NUCLEAR_PHOTO_DEPTH_SOURCE="+json.dumps({
        "n":c["n_matched_original_taxa"],
        "observed":c["observed_standardized_rain_then_temperature_residual"],
        "exact_p":c["unrelated_one_sided_exact_p"],
        "phylogenetic_p":{k:v["H1_posthoc_one_sided_p"] for k,v in
                         c["within_each_own_nuclear_tree_conditional_Mk"].items()},
        "EA20_photo_depth":out["EA20_separate_unmerged_alternative_photo_depth"],
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

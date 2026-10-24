"""Expand original S1 image-angle WN/BM/OU1 to all 202 source-name-matched tips.

Input is the source-verified 251 original taxon-name x accepted-key rows
from 170066 native photographs. Unlike legacy 194 accepted-key matches
(3 name collisions), each terminal is assigned ONLY its exact original
recorded species name. 202 matching original V.PhyloMaker/GBOTB placement
tips can be analyzed, versus prior 191 alias-clean accepted keys.
They are placement hypotheses, never independent Cirsium nuclear trees.

Use the frozen original five CHELSA covariates, the same EAzami finite-root
OU covariance, and WN independent-tips null, all on identical source-name
tip sets. Do not label EXIF vertical as gravitational nodding.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_history import (
    audited_pruned_tree, shared_path_matrix,
)
from analysis.v4.t2_all244_ou_bm_orientation import (
    _design_and_values, evaluate_models, EXPOSURES,
)

TOTAL_NATIVE=170066
SOURCE_UNITS=251
ACCEPTED_KEYS=244
MATCHED_TIPS=202
ORIGINAL_ALIAS_CLEAN=191
TREE_SCENARIOS=52
CASE_NAMES=(
    "matched202_all_images_mean_five_climates",
    "matched202_all_images_mean_intercept",
    "matched202_all_images_median_five_climates",
    "original191_noncollided_mean_five_climates",
    "matched202_n_ge_5_mean_five_climates",
    "matched202_n_ge_10_mean_five_climates",
)
SOURCE_COLS={"accepted_key","original_source_name","original_photo_count",
             "mean_image_angle_deg","median_image_angle_deg",
             "WCVP_key_collides_source_names",
             "exact_source_name_has_frozen_placement_tip",*EXPOSURES}


def read_source(path):
    df=pd.read_csv(path,dtype={"accepted_key":str},low_memory=False)
    if SOURCE_COLS-set(df):
        raise ValueError("Source-name photo/climate provenance columns missing")
    if (len(df)!=SOURCE_UNITS or
        df.accepted_key.nunique()!=ACCEPTED_KEYS or
        df.original_photo_count.sum()!=TOTAL_NATIVE or
        df.duplicated(["accepted_key","original_source_name"]).any()):
        raise ValueError("Original 251 source-name 244-key 170066-image corpus changed")
    matching=df.loc[df.exact_source_name_has_frozen_placement_tip].copy()
    clean=matching.loc[~matching.WCVP_key_collides_source_names].copy()
    if (len(matching)!=MATCHED_TIPS or len(clean)!=ORIGINAL_ALIAS_CLEAN or
        matching.original_source_name.nunique()!=MATCHED_TIPS):
        raise ValueError("202 exact original terminal / 191 alias-clean identity coverage changed")
    if matching.original_photo_count.min()<1 or \
            len(df.groupby("accepted_key").filter(lambda x:len(x)>1).accepted_key.unique())!=7:
        raise ValueError("Source-name separation or S1 photograph depth invalid")
    def to_rows(frame):
        out=[]
        for _,r in frame.sort_values("original_source_name").iterrows():
            out.append({
                "tip":str(r.original_source_name),
                "accepted_key":str(r.accepted_key),
                "n_eligible_S1_images":int(r.original_photo_count),
                "mean_S1_image_axis_deg":float(r.mean_image_angle_deg),
                "median_S1_image_axis_deg":float(r.median_image_angle_deg),
                "collision":bool(r.WCVP_key_collides_source_names),
                "S1_five_climate_means":{e:float(r[e]) for e in EXPOSURES},
            })
        return out
    allrows=to_rows(matching)
    return {
        "matched202":allrows,
        "clean191":to_rows(clean),
        "n_ge_5":[r for r in allrows if r["n_eligible_S1_images"]>=5],
        "n_ge_10":[r for r in allrows if r["n_eligible_S1_images"]>=10],
    }, hashlib.sha256(Path(path).read_bytes()).hexdigest()


def designs(panel):
    return (
        ("matched202_all_images_mean_five_climates","matched202","mean",True,.05),
        ("matched202_all_images_mean_intercept","matched202","mean",False,.05),
        ("matched202_all_images_median_five_climates","matched202","median",True,.05),
        ("original191_noncollided_mean_five_climates","clean191","mean",True,.05),
        ("matched202_n_ge_5_mean_five_climates","n_ge_5","mean",True,.05),
        ("matched202_n_ge_10_mean_five_climates","n_ge_10","mean",True,.05),
    )


def analyse(path,tree_limit=52):
    if not 1<=tree_limit<=TREE_SCENARIOS:
        raise ValueError("Unexpected original placement ensemble scope")
    panels,source_sha=read_source(path)
    full=panels["matched202"]
    tip_order=[r["tip"] for r in full]
    ix={tip:i for i,tip in enumerate(tip_order)}
    configs=[]
    for case,group,phenotype,adjust,nugget in designs(panels):
        rows=panels[group]
        if len(rows)<30:
            raise ValueError("Too few original source-named tips")
        seq=[ix[r["tip"]] for r in rows]
        y,x,scale=_design_and_values(rows,phenotype,adjust)
        configs.append((case,group,rows,seq,y,x,scale,nugget))
    trees=restore_trees()
    if len(trees)!=TREE_SCENARIOS:
        raise ValueError("Historical 52 placement sensitivity changed")
    results=[]
    for tree_i,(name,tree) in enumerate(trees[:tree_limit]):
        subtree,_=audited_pruned_tree(tree,tip_order)
        C=shared_path_matrix(subtree,tip_order)[0]
        for case,group,rows,seq,y,x,scale,nugget in configs:
            c=C[np.ix_(seq,seq)]
            ans=evaluate_models(y,x,c,scale,nugget)
            results.append({
                "case":case,"cohort":group,
                "tree_name":name,"tree_scenario_index":tree_i,
                "n_source_name_tips":len(rows),
                "n_collided_WCVP_taxon_keys_included":len({
                    r["accepted_key"] for r in rows if r["collision"]}),
                **ans,
            })
        if tree_i in (0,9,19,29,39,49,51):
            print("T2_S1_ALL202_NAME_RESOLVED_PLACEMENT_SCENARIOS="+str(tree_i+1),flush=True)
    bycase={}
    for case in CASE_NAMES:
        panel=[x for x in results if x["case"]==case]
        if len(panel)!=tree_limit:
            raise ValueError("Placement ensemble lost a condition")
        fits=Counter(r["best_model_AICc"] for r in panel)
        bycase[case]={
            "n_tips":panel[0]["n_source_name_tips"],
            "n_original_placement_scenarios":tree_limit,
            "wins":{m:fits[m] for m in ("WN","BM","OU1")},
            "OU_minus_BM_AICc_range":[float(min(z["OU1_minus_BM_AICc"] for z in panel)),
                                       float(max(z["OU1_minus_BM_AICc"] for z in panel))],
            "n_ou_alpha_upper_bound":sum(x["OU1_alpha_grid_edge"] for x in panel),
            "n_collided_source_keys_included":panel[0]["n_collided_WCVP_taxon_keys_included"],
        }
    return {
        "schema":"azami_v4_S1_original251_to_exact202_on52_placement_WN_BM_OU1_v1",
        "status":"ACTUAL_ORIGINAL_SOURCE_NAME_202_PHYLOGENETIC_PLACEMENT_SENSITIVITY",
        "source_sha256":source_sha,
        "original_photo_observations":TOTAL_NATIVE,
        "original_accepted_WCVP_keys":ACCEPTED_KEYS,
        "original_source_name_operational_units":SOURCE_UNITS,
        "exact_original_source_named_placement_tips":MATCHED_TIPS,
        "previous_clean_no_collided_keys_191":ORIGINAL_ALIAS_CLEAN,
        "source_name_tips_not_in_original_placement_tree":SOURCE_UNITS-MATCHED_TIPS,
        "n_nonindependent_placement_hypotheses":tree_limit,
        "original_fixed_five_CHELSA_covariates":list(EXPOSURES),
        "unit_is_photo_EXIF_vertical_not_botanical_gravity":True,
        "grafted_tree_not_independent_nuclear_species_phylogeny":True,
        "no_new_real_nuclear_tip_added_to_EAzami294":True,
        "no_independent_nodding_origin_count_estimated":True,
        "model_winner_summary":bycase,
        "each_source_original_placement_model":results,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-251",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--tree-limit",type=int,default=52)
    a=p.parse_args()
    result=analyse(a.source_251,a.tree_limit)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2,
             ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print("T2_S1_ALL202_SOURCE_NAMED_WN_BM_OU_SUMMARY="+json.dumps(
          result["model_winner_summary"],ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

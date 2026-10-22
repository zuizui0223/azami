"""Independent NA nuclear-tree BM/OU1 benchmark for EXACT same original S1 taxa.

Keep EAzami-style Gaussian phylogenetic models separate from the 52 synthetic
placement hypotheses. Input: original S1 source-authenticated aggregate with 42
author-published North American nuclear tips. Does NOT infer trait gain events.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

from Bio import Phylo
import numpy as np

from analysis.v4.s1_continuous_na_nuclear_signal import source_rows
from analysis.v4.s1_continuous_phylogenetic_history import shared_path_matrix
from analysis.v4.s1_continuous_phylogenetic_signal import verify_source
from analysis.v4.t2_source_blocked_five_climate import author_nuclear_tree
from analysis.v4.t2_all244_ou_bm_orientation import evaluate_models
from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES

NA_TIPS_EXPECTED=42


def prepare(source:Path,original_194:Path):
    """Bridge exact 42 concepts, but NEVER reuse a median as a mean."""
    bridge=source_rows(source)
    continuous=verify_source(original_194)
    verified194={row["tip"]:row for row in continuous}
    tree=author_nuclear_tree()
    if (len(bridge)!=NA_TIPS_EXPECTED or
        len({x["tip"] for x in bridge})!=NA_TIPS_EXPECTED):
        raise ValueError("Published nuclear 42-tip S1 overlap drifted")
    rows=[]
    for row in bridge:
        item=verified194.get(row["tip"])
        if item is None or item["source_concept_collision"]:
            raise ValueError("Independent 42-tip S1 original measured source concept missing")
        if (str(item["accepted_key"])!=str(row["accepted_key"])
            or int(item["n_eligible_images"])!=int(row["n_eligible_images"])
            or abs(float(item["image_median_deg"])-float(row["image_median_deg"]))>1e-6):
            raise ValueError("Exact source median/observation/concept mismatch")
        orig=item["S1_five_climate_means"]
        bridge_climate=row["S1_five_climate_means"]
        if set(orig)!=set(EXPOSURES) or set(bridge_climate)!=set(EXPOSURES) or any(
            abs(float(orig[e])-float(bridge_climate[e]))>1e-6 for e in EXPOSURES):
            raise ValueError("Independent 194-source five-climate crosswalk mismatch")
        rows.append({
            **row,
            "image_mean_deg":float(item["image_mean_deg"]),
            "image_median_deg":float(item["image_median_deg"]),
        })
    source_by={row["tip"]:row for row in rows}
    nuclear=copy.deepcopy(tree)
    for tip in list(nuclear.get_terminals()):
        if tip.name not in source_by:
            nuclear.prune(tip)
    if {x.name for x in nuclear.get_terminals()}!=set(source_by):
        raise ValueError("Author nuclear tree lacks exact S1 orientation concepts")
    tip_order=sorted(source_by)
    C,_,_,_=shared_path_matrix(nuclear,tip_order)
    rows=[source_by[n] for n in tip_order]
    return rows,C


def run(source:Path,original_194:Path):
    rows,C=prepare(source,original_194)
    ymean=np.array([float(r["image_mean_deg"]) for r in rows])
    ymedian=np.array([float(r["image_median_deg"]) for r in rows])
    assert 0<=ymean.min()<=ymean.max()<=180
    climate=np.array([[float(r["S1_five_climate_means"][e]) for e in EXPOSURES]
                      for r in rows],float)
    if not np.isfinite(climate).all() or np.any(climate.std(axis=0)<=1e-10):
        raise ValueError("No complete fixed 5-exposure original S1 NA42 climate matrix")
    X5=np.column_stack([np.ones(len(rows)),
                        (climate-climate.mean(axis=0))/climate.std(axis=0)])
    if np.linalg.matrix_rank(X5)!=6:
        raise ValueError("NA42 original five-climate fit not identifiable")
    cases=[
      ("mean",False,0.05,ymean,np.ones((42,1))),
      ("mean",True,0.05,ymean,X5),
      ("median",True,0.05,ymedian,X5),
      ("mean",True,0.25,ymean,X5),
    ]
    results=[]
    for endpoint,adjust,nugget,y,X in cases:
        scale=float(y.std(ddof=1))
        if scale<=0:
            raise ValueError("NA42 measured image angle lacks biological spread")
        fit=evaluate_models((y-y.mean())/scale,X,C,scale,nugget)
        results.append({
            "endpoint":"S1_original_taxon_"+endpoint+"_EXIF_vertical_image_angle",
            "simultaneous_five_climates":bool(adjust),
            "nugget":nugget,
            **fit,
        })
    return {
        "schema":"azami_V4_S1_NA42_PUBLISHED_NUCLEAR_TRUE_MEAN_AND_MEDIAN_WN_BM_OU1_20261022_v2",
        "source_S1_joint_run":39584744194,
        "source_S1_original_194_continuous_run":39594257539,
        "n_taxa_with_distinct_original_image_mean_median":int(np.sum(abs(ymean-ymedian)>1e-9)),
        "published_nuclear_tree_DOI":"10.1086/724310",
        "source_original_S1_taxa":244,
        "matched_independent_nuclear_taxa":42,
        "all_published_biological_nodding_scores_used_as_outcomes":False,
        "trait":"EXIF_image_vertical_relative_continuous_angle_NOT_biological_gravity_posture",
        "branch_length_unit":"published nuclear tree author ultrametric units, not independently calibrated paleotime",
        "same_five_annual_climate_predictors":list(EXPOSURES),
        "fits":results,
        "contrast_to_52_grafted_trees":"Same measured 42 tips; nuclear source validates tree-source sensitivity but does not generalize to remaining S1 taxa.",
        "claim_stop":{
            "not_global_Cirsium_nuclear_tree":True,
            "OU_alpha_not_fitness_selection_strength":True,
            "no_ancestral_nodding_gain_count":True,
            "no_adaptive_radiation_or_rain_protection_causality_proven":True
        }
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--joint",required=True,type=Path)
    p.add_argument("--original-194",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    result=run(a.joint,a.original_194)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,
                                sort_keys=True,allow_nan=False)+"\n")
    print("T2_NA42_TRUE_NUCLEAR_WN_BM_OU="+json.dumps([
        {"endpoint":r["endpoint"],
         "adjusted":r["simultaneous_five_climates"],
         "winner":r["best_model_AICc"],
         "deltaOU_BM":r["OU1_minus_BM_AICc"]}
        for r in result["fits"]],ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

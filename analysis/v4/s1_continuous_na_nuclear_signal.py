"""Control for tree source on EXACT SAME S1 continuous photo angles in NA42.

Does replacing V.PhyloMaker2 grafted placements with the original author's
2023 63-tip molecular NUCLEAR species tree recover stronger phylogenetic signal
for the exact same 42 S1-photographed North American Cirsium taxa?

No T2 photo-angle binarization, and no claim of global Cirsium nuclear history.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_signal import (
    LambdaProfile,model_design,covariance_for,
)
from analysis.v4.t2_source_blocked_five_climate import (
    authenticate_input, RECEIPT, NA,author_nuclear_tree,
)

NUCLEAR_SOURCE_DOI="10.1086/724310"
JOINT_SOURCE_RUN=39584744194


def source_rows(joint_path):
    raw=json.loads(Path(joint_path).read_text(encoding="utf-8"))
    receipt=json.loads(RECEIPT.read_text(encoding="utf-8"))
    all46=authenticate_input(raw,receipt)
    source=[r for r in all46 if r["source"]==NA]
    if len(source)!=42 or sum(r["state"]=="nodding" for r in source)!=3:
        raise ValueError("North-American 42 known botanical status cohort changed")
    selected=[]
    for r in source:
        tip=r["source_taxon"].replace(" ","_")
        yy=float(r["s1_image_orientation_median"])
        n=int(r["s1_observations"])
        climate=r["s1_five_climate_centroids"]
        if n<1 or not (0<=yy<=180) or set(climate)!=set(EXPOSURES):
            raise ValueError("Incomplete original S1 photo angle and climate source")
        selected.append({
           "tip":tip,"accepted_key":str(r["accepted_key"]),
           "n_eligible_images":n,
           "source_concept_collision":False,
           "image_mean_deg":yy,
           "image_median_deg":yy,
           "S1_five_climate_means":climate,
           "external_botanical_state_supplement_only":r["state"],
        })
    if len({x["tip"] for x in selected})!=42:
        raise ValueError("Source provenance duplicated 42 molecular tips")
    return selected


def one_tree(tree,label,rows,nugget=.05):
    ordered,X0=model_design(rows,False)
    _,X5=model_design(rows,True)
    C=covariance_for(tree,ordered)
    y=np.array([x["image_median_deg"] for x in ordered],dtype=float)
    return {"tree":label,"n_tips":len(rows),
        "unadjusted":LambdaProfile(C,X0,nugget).fit(y),
        "five_climate_adjusted":LambdaProfile(C,X5,nugget).fit(y)}


def analyze(rows,nuclear,trees,*,bootstrap_reps=399):
    if len(rows)!=42 or len(trees)!=52:
        raise ValueError("Only 42 exact S1-author-nuclear tips and 52 grafted sensitivity trees")
    ref=one_tree(nuclear,"independent_published_original_NA_nuclear",rows)
    ordered,X0=model_design(rows,False)
    _,X5=model_design(rows,True)
    y=np.array([x["image_median_deg"] for x in ordered],dtype=float)
    nuclear_C=covariance_for(nuclear,ordered)
    boots={
        "nuclear_unadjusted":LambdaProfile(nuclear_C,X0,.05).bootstrap_lambda0(
            y,reps=bootstrap_reps,seed=20261021),
        "nuclear_S1_five_climate_adjusted":LambdaProfile(nuclear_C,X5,.05).bootstrap_lambda0(
            y,reps=bootstrap_reps,seed=20261022),
    }
    # The generic Pagel helper's fixed-tree bootstrap wording describes
    # placement sensitivity by default. Override it for this source-authored
    # independently estimated North-American nuclear species tree.
    for record in boots.values():
        record["scope"]=(
            "One independently published author nuclear Cirsium species tree "
            "(Siniscalchi et al. 2023) pruned to exactly the same 42 S1 photo "
            "taxa; fixed topology, branch lengths, plug-in error variance and "
            "source photo-angle measurement; NOT 52 V.PhyloMaker2 placements, "
            "nor genus-wide ancestry, selection, or environmental causation."
        )
    # Same data, different tree source. Report all placement scenarios, not
    # 52 independent confirmations of any nuclear-tree result.
    fits=[one_tree(tree,name,rows) for name,tree in trees]
    def envelope(key,field):
        v=np.array([x[key][field] for x in fits],dtype=float)
        return {
            "min":float(v.min()),"median":float(np.median(v)),
            "max":float(v.max()),
            "n_positive_lambda_ge_0_1":int(np.sum(
               [x[key]["lambda_ML"]>=.1 for x in fits])),
            "n_deltaAIC_over_2":int(np.sum(
               [x[key]["deltaAIC_best_vs_lambda0"]>2 for x in fits])),
        }
    summary={
      "status":"REAL_S1_SAME_NA42_CONTINUOUS_ANGLE_INDEPENDENT_NUCLEAR_VS_52_GRAFTED_PHYLOGENETIC_SIGNAL",
      "source_genuine_nuclear_tree_DOI":NUCLEAR_SOURCE_DOI,
      "S1_native_photo_eligibility_total":170066,
      "S1_native_taxon_keys_total":244,
      "S1_available_GBOTB_matched_tip_count":194,
      "same_original_photo_phenotype_nuclear_vs_placements_n_taxa":42,
      "source_original_head_angle":"median_EXIF_image_vertical_degrees_NOT_gravity",
      "external_3_botanical_nodding_states_ignored_as_phenotype_response":True,
      "S1_five_original_covar":list(EXPOSURES),
      "fixed_tip_nugget_fraction":.05,
      "published_independent_NA_nuclear":ref,
      "matched_52_placement_phenotype_unadjusted_lambda":envelope("unadjusted","lambda_ML"),
      "matched_52_placement_phenotype_fiveclimate_lambda":envelope("five_climate_adjusted","lambda_ML"),
      "same_NA42_nuclear_bootstrap":boots,
      "NA42_is_not_global_Cirsium_independent_phylogeny":True,
      "52_placements_not_independent_evolutionary_replicates":True,
      "weak_signal_cannot_prove_no_genetic_structure":True,
      "nodding_gravity_independent_evolution_or_rain_selection_uninferred":True,
    }
    return {"summary":summary,"per_placement_results":fits,
            "source_frame":"S1 original measured image angle and same 5 climate covariates, with third-party botanical scores only as provenance"}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--joint",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--bootstrap-reps",type=int,default=399)
    opts=p.parse_args()
    if opts.bootstrap_reps<20 or opts.bootstrap_reps>5000:
        raise ValueError("Unsupported bootstrap depth")
    rows=source_rows(opts.joint)
    nuclear=author_nuclear_tree()
    placements=restore_trees()
    result=analyze(rows,nuclear,placements,bootstrap_reps=opts.bootstrap_reps)
    dest=Path(opts.out);dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False,
                               allow_nan=False)+"\n",encoding="utf-8")
    print("S1_NA42_SAME_PHENOTYPE_NUCLEAR_VS_GRAFTED_PHYLOGENETIC_SIGNAL="+
          json.dumps(result["summary"],sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

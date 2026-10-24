"""S1 image-angle phylogenetic signal: true NA nuclear vs SAME-tip 52 placements.

A matched original-source comparison, not a replacement of 194-species global
S1 results. Uses the SHA-verified S1 194-species continuous taxon aggregate
from completed GitHub Actions 39594257539 (NO private photo observations),
and the published source-pinned Siniscalchi et al. 2023 63-tip nuclear tree.
All five S1 CHELSA exposures fixed. No image-angle->botanical U/D conversion.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_signal import (
    verify_source, model_design, LambdaProfile, covariance_for,
)
from analysis.v4.t2_source_blocked_five_climate import author_nuclear_tree
from analysis.v4.t2_pgls_placement_screen import (
    AUTHOR_NUCLEAR_TREE_GIT_BLOB, SOURCE_COMMIT,
)

SOURCE_RUN=39594257539
SOURCE_ARTIFACT_ID=12343166776
SOURCE_ARTIFACT_SHA="sha256:345617fd8fef6bdb41a06a1a03fa8a776c5ad48ba13c2aa971deb1e929cbf8cb"
EXPECTED_NA_MATCH=42
NUCLEAR_SOURCE="Siniscalchi_Ackerfield_Folk_2023_DOI_10.1086/724310"
NUGGETS=(.05,.25)


def match_original_S1_to_nuclear(s1, original_nuclear):
    original_tipset={t.name for t in original_nuclear.get_terminals()}
    s1_tipset={r["tip"] for r in s1}
    if len(original_tipset)!=63:
        raise ValueError("Pinned independent NA nuclear tree must have 63 unique tips")
    match=original_tipset & s1_tipset
    if len(match)!=EXPECTED_NA_MATCH:
        raise ValueError("Published original nuclear and real S1 matched-tip count changed")
    rows=[r for r in s1 if r["tip"] in match]
    if len(rows)!=42 or any(r["source_concept_collision"] for r in rows):
        raise ValueError("NA nuclear source comparison contains identity collision")
    return rows


def fit_one_tree(tree,rows,*,root_label,source_type,nuggets=NUGGETS,
                 bootstrap_reps=399,bootstrap_only=False,seed=20261021):
    ordered,X0=model_design(rows,False)
    _,X5=model_design(rows,True)
    C=covariance_for(tree,ordered)
    names=[r["tip"] for r in ordered]
    y={
       "species_mean":np.array([r["image_mean_deg"] for r in ordered],dtype=float),
       "species_median":np.array([r["image_median_deg"] for r in ordered],dtype=float),
    }
    report=[]
    bootstrap={}
    for nugget in nuggets:
        for design_label,X in (("no_climate",X0),("all_five_original_S1_climates",X5)):
            evaluator=LambdaProfile(C,X,nugget)
            for response,observed in y.items():
                fit=evaluator.fit(observed)
                report.append({
                   "tree_source":source_type,"scenario":root_label,
                   "n_taxa":len(rows),"nugget_fraction":float(nugget),
                   "model":design_label,"response":response,
                   "lambda_ML":fit["lambda_ML"],
                   "LRT_against_lambda0":fit["likelihood_ratio_against_lambda0"],
                   "delta_AIC_against_lambda0":fit["deltaAIC_best_vs_lambda0"],
                   "covariance_model":"shared_root_Brownian_Pagel_lambda_0_to_1",
                   "not_true_gravity_nodding":True,
                })
                if (bootstrap_only and abs(nugget-.05)<1e-12
                       and response=="species_mean"):
                    key=f"{design_label}__{response}__nugget{nugget}"
                    bootstrap[key]=evaluator.bootstrap_lambda0(
                        observed,reps=bootstrap_reps,
                        seed=seed+(0 if design_label=="no_climate" else 1))
                    bootstrap[key]["scope"]=(
                        "Original author-published North American nuclear tree "
                        "pruned to identical 42 S1-measured taxa; fixed branch "
                        "lengths, fixed nugget, normal residuals, and fitted "
                        "lambda; not a genus-wide origin or photo-camera "
                        "accuracy test")
    return report,bootstrap


def compare(s1, published_nuclear, trees, *, bootstrap_reps=399):
    na=match_original_S1_to_nuclear(s1,published_nuclear)
    if len(trees)!=52:
        raise ValueError("52 original placement scenarios required")
    all_tips={t.name for t in trees[0][1].get_terminals()}
    if len(all_tips)!=len(trees[0][1].get_terminals()):
        raise ValueError("Duplicate placement tree tip")
    for _,t in trees[1:]:
        if {tip.name for tip in t.get_terminals()}!=all_tips:
            raise ValueError("Historical placements must have same tips")
    if not {r["tip"] for r in na} <= all_tips:
        raise ValueError("Independent nuclear S1 cohort mismatched placement tree")

    cohorts={"all42":na}
    cohorts["at_least5"]=[r for r in na if int(r["n_eligible_images"])>=5]
    cohorts["at_least10"]=[r for r in na if int(r["n_eligible_images"])>=10]
    if any(len(z)<10 for z in cohorts.values()):
        # 5 predictors + intercept + source limitations: avoid tiny submodels
        cohorts={k:v for k,v in cohorts.items() if len(v)>=10}

    details=[]
    bootstrap={}
    for cohort_name,subset in cohorts.items():
        for index,(label,tree) in enumerate([("original_NA_nuclear63",published_nuclear)]+trees):
            kind="independent_published_NA_nuclear" if index==0 else "GBOTB_grafted_same42_sensitivity"
            result,calibrated=fit_one_tree(tree,subset,root_label=label,source_type=kind,
                 bootstrap_reps=bootstrap_reps,
                 bootstrap_only=(index==0 and cohort_name=="all42"),
                 seed=20261021+10*list(cohorts).index(cohort_name))
            for row in result:
                row["cohort"]=cohort_name
            details.extend(result)
            if calibrated:
                bootstrap[cohort_name]=calibrated

    def aggregate(records):
        lam=[r["lambda_ML"] for r in records]
        lrt=[r["LRT_against_lambda0"] for r in records]
        return {
          "lambda_ML_range":[float(min(lam)),float(max(lam))],
          "lambda_ML_median":float(np.median(lam)),
          "n_scenarios":len(records),
          "n_lambda_gt0_1":sum(r>0.1 for r in lam),
          "n_deltaAIC_gt2":sum(r>4. for r in lrt),
          "LRT_range":[float(min(lrt)),float(max(lrt))],
        }

    summaries={}
    for key,subset in cohorts.items():
        section={}
        for response in ("species_mean","species_median"):
            for design in ("no_climate","all_five_original_S1_climates"):
                for nugget in NUGGETS:
                    matching=[r for r in details if r["cohort"]==key
                        and r["response"]==response and r["model"]==design
                        and r["nugget_fraction"]==nugget]
                    assert len(matching)==53
                    nuclear=[r for r in matching if r["tree_source"]=="independent_published_NA_nuclear"]
                    graft=[r for r in matching if r["tree_source"]!="independent_published_NA_nuclear"]
                    suffix=f"{response}|{design}|nugget{nugget}"
                    section[suffix]={
                       "n_tip_cohort":len(subset),
                       "nuclear":nuclear[0],
                       "grafted_52":aggregate(graft),
                    }
        summaries[key]=section
    return {
       "summary":{
         "status":"VERIFIED_PUBLISHED_NA_NUCLEAR_VS_GRAFTED_CONTINUOUS_S1_SIGNAL",
         "global_original_S1_native_eligible_taxa":244,
         "global_S1_phylogeny_covered_continuous_taxa":194,
         "original_published_na_nuclear_tip_count":63,
         "same_NA_S1_image_measured_taxa":42,
         "original_source_run_id":SOURCE_RUN,
         "original_source_artifact_id":SOURCE_ARTIFACT_ID,
         "original_source_zip_sha256":SOURCE_ARTIFACT_SHA,
         "NA_molecular_source_commit":SOURCE_COMMIT,
         "NA_molecular_original_Git_blob":AUTHOR_NUCLEAR_TREE_GIT_BLOB,
         "NA_molecular_DOI":"10.1086/724310",
         "cohort_image_count_sizes":{k:len(v) for k,v in cohorts.items()},
         "same_five_S1_climate_covariates": ["chelsa_bio01","chelsa_bio04",
                  "chelsa_bio12","chelsa_bio15","chelsa_rsds_mean"],
         "n_historical_GBOTB_grafts":52,
         "same_tip_comparison":True,
         "no_frozen_observational_S1_refit":True,
         "photo_angle_is_not_biological_gravity":True,
         "no_true_nodding_origin_count":True,
         "nuclear_tree_NA_scope_only":True,
         "root_age_or_tree_uncertainty_not_integrated":True,
         "nuclear_v_grafted_panels":summaries,
         "nuclear_original_fixed_tree_lambda0_parametric_bootstrap":bootstrap,
       },
       "full_53_tree_fit_details":details
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",required=True)
    ap.add_argument("--out",required=True)
    ap.add_argument("--bootstrap-reps",type=int,default=399)
    args=ap.parse_args()
    if not 20<=args.bootstrap_reps<=5000:
        raise ValueError("Bootstrap replicate count out of range")
    s1=verify_source(args.source)
    nuclear=author_nuclear_tree()
    trees=restore_trees()
    result=compare(s1,nuclear,trees,bootstrap_reps=args.bootstrap_reps)
    result["original_S1_continuous_full_source_file_sha256"]=hashlib.sha256(
       Path(args.source).read_bytes()).hexdigest()
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",
                   encoding="utf-8")
    panels=result["summary"]["nuclear_v_grafted_panels"]
    print("S1_NUCLEAR_V_GRAFTED_SAME42_SIGNAL_SUMMARY="+json.dumps({
       "status":result["summary"]["status"],
       "cohort_sizes":result["summary"]["cohort_image_count_sizes"],
       "first_na_nuclear":{
         cohort: {key:value for key,value in entries.items()
          if key in ("species_mean|no_climate|nugget0.05",
                     "species_mean|all_five_original_S1_climates|nugget0.05")}
         for cohort,entries in panels.items()},
       "null_calibration":result["summary"]["nuclear_original_fixed_tree_lambda0_parametric_bootstrap"],
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

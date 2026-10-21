"""Does S1 annual-rain angle association survive free-lambda phylogeny control?

Continuous trait, same original 194 source-measured S1 image-angle taxa and
same fixed five climate exposures BIO1/BIO4/BIO12/BIO15/RSDS. Include strict
191 no-WCVP-merge cohort, and original >=5 (171) / >=10 (143) photographs
per taxon precision sensitivity. Compare lambda=0, lambda_MLE, lambda=1
using the exact same taxon-level species-mean GLS design on the 52 frozen
GBOTB grafted placement scenarios. This is a complementary taxon-equal
diagnostic, NOT a replacement for the original S1 mixed/observer model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.native_phylogeny import restore_trees
from analysis.v4.s1_continuous_phylogenetic_signal import (
    LambdaProfile, covariance_for, model_design, taxon_sets, verify_source,
)

SOURCE_RUN=39594257539
SOURCE_ARTIFACT=12343166776
SOURCE_DIGEST="sha256:345617fd8fef6bdb41a06a1a03fa8a776c5ad48ba13c2aa971deb1e929cbf8cb"
INC=np.array([1.,100.,100.,10.,1.],dtype=float)
FIXED_MODEL_CLASS=("independent_lambda0","Pagel_lambda_MLE","Brownian_lambda1")
NUGGETS=(.05,.25)


def estimate_coefficients(profile,y,lam,sd_native):
    """Conditional GLS t intervals with ML/forced lambda, no lambda CI integration."""
    if not 0<=lam<=1:
        raise ValueError("Lambda outside permitted [0,1]")
    if sd_native.shape!=(5,) or np.any(sd_native<=0):
        raise ValueError("Original S1 climate standardization invalid")
    obs=np.asarray(y,float)
    Uy=profile.rotation.T@(obs/profile.depth)
    v=(1+profile.nugget)+lam*(profile.eig-1)
    if np.min(v)<=0:
        raise ValueError("Invalid shared covariance eigenvalue")
    weighted_X=profile.UX/np.sqrt(v[:,None])
    weighted_y=Uy/np.sqrt(v)
    information=weighted_X.T@weighted_X
    if np.linalg.cond(information)>1e11:
        raise ValueError("Unidentified original five predictor covariance matrix")
    beta=np.linalg.solve(information,weighted_X.T@weighted_y)
    resid=weighted_y-weighted_X@beta
    df=profile.n-profile.p
    mse=float(resid@resid/df)
    beta_se=np.sqrt(np.diag(np.linalg.inv(information))*mse)
    crit=float(student_t.ppf(.975,df))
    scales=INC/sd_native
    estimated={e:{
       "effect_image_deg_per_original_increment":float(beta[i+1]*scales[i]),
       "conditional_t_95pct":[float((beta[i+1]-crit*beta_se[i+1])*scales[i]),
                              float((beta[i+1]+crit*beta_se[i+1])*scales[i])],
       "increment":float(INC[i]),
       "original_climate_standard_deviation":float(sd_native[i]),
       "interval_conditional_fixed_tree_lambda_and_nugget":True,
     } for i,e in enumerate(EXPOSURES)}
    return estimated


def model_one_tree(cohort,original,*,nuggets=NUGGETS):
    rows,X=model_design(cohort,True)
    names=[r["tip"] for r in rows]
    y={
       "species_mean":np.array([r["image_mean_deg"] for r in rows],float),
       "species_median":np.array([r["image_median_deg"] for r in rows],float),
    }
    climate=np.array([[float(r["S1_five_climate_means"][v]) for v in EXPOSURES]
                      for r in rows])
    sd=climate.std(axis=0)
    C=covariance_for(original,rows)
    records=[]
    for nugget in nuggets:
        profile=LambdaProfile(C,X,nugget)
        for response,angle in y.items():
            m=profile.fit(angle)
            for model,lam in zip(FIXED_MODEL_CLASS,(0.,m["lambda_ML"],1.)):
                effect=estimate_coefficients(profile,angle,lam,sd)
                records.append({
                  "cohort_n":len(rows),
                  "response":response,
                  "fit":model,
                  "lambda":float(lam),
                  "nugget_fraction":nugget,
                  "LRT_model_lambda_MLE_against_lambda0":m["likelihood_ratio_against_lambda0"],
                  "source_five_climate_variables":list(EXPOSURES),
                  "coefficients":effect,
                  "not_observation_level_S1":True,
                  "not_gravity_head_angle":True,
                })
    return records


def compare(rows,trees,*,nuggets=NUGGETS):
    cohorts=taxon_sets(rows)
    expected={"original194":194,"collision_free191":191,
              "collision_free_at_least5":171,"collision_free_at_least10":143}
    if {key:len(group) for key,group in cohorts.items()} != expected or len(trees)!=52:
        raise ValueError("Frozen 194/191/171/143 S1 cohort and 52 trees changed")
    fit_records=[]
    for label,group in cohorts.items():
        for name,tree in trees:
            model=model_one_tree(group,tree,nuggets=nuggets)
            for r in model:
                r["cohort"]=label
                r["tree_placement"]=name
                fit_records.append(r)
    def stats(cohort,trait,model,nugget,exposure):
        matches=[r for r in fit_records if r["cohort"]==cohort
                 and r["response"]==trait and r["fit"]==model
                 and r["nugget_fraction"]==nugget]
        if len(matches)!=52:
            raise ValueError("Original PGLS tree scenario output incomplete")
        original=[r["coefficients"][exposure] for r in matches]
        estimates=[r["effect_image_deg_per_original_increment"] for r in original]
        lows=[r["conditional_t_95pct"][0] for r in original]
        highs=[r["conditional_t_95pct"][1] for r in original]
        lam=[r["lambda"] for r in matches]
        return {
           "n_nonindependent_scenario_trees":52,
           "range_effect_image_deg_per_exposure_increment":[float(min(estimates)),float(max(estimates))],
           "range_conditional_lower95":[float(min(lows)),float(max(lows))],
           "range_conditional_upper95":[float(min(highs)),float(max(highs))],
           "all_52_effects_positive":all(v>0 for v in estimates),
           "n_pointwise_95pct_positive":sum(v>0 for v in lows),
           "n_pointwise_95pct_negative":sum(v<0 for v in highs),
           "n_pointwise_interval_crosses_zero":sum(a<=0<=b for a,b in zip(lows,highs)),
           "lambda_scenario_median":float(np.median(lam)),
           "lambda_scenario_range":[float(min(lam)),float(max(lam))],
           "not_52_independent_studies":True,
        }
    panels={}
    for cohort in expected:
        for trait in ("species_mean","species_median"):
            for model in FIXED_MODEL_CLASS:
                for nugget in nuggets:
                    for exposure in ("chelsa_bio12","chelsa_bio15"):
                        key=f"{cohort}|{trait}|{model}|nugget{nugget}|{exposure}"
                        panels[key]=stats(cohort,trait,model,nugget,exposure)
    summary={
      "status":"EMPIRICAL_S1_CONTINUOUS_FIVE_CLIMATE_PAGEL_LAMBDA_PROFILED_EFFECTS",
      "original_S1_image_angle_eligible_observations":170066,
      "original_S1_accepted_taxon_keys":244,
      "original_S1_covar_species_analysis_source":194,
      "S1_collision_free_main_source":191,
      "S1_more_precisely_photographed_taxon_sensitivity":{"n_at_least5":171,"n_at_least10":143},
      "52_PGLS_historical_placements_not_independent_nuclear_trees":True,
      "same_five_S1_CHELSA_exposures":list(EXPOSURES),
      "registered_precipitation_axes":["chelsa_bio12","chelsa_bio15"],
      "model_comparators":list(FIXED_MODEL_CLASS),
      "fixed_tip_nuggets":list(nuggets),
      "no_observational_S1_or_V3_reestimate":True,
      "photographic_angle_not_real_gravity":True,
      "no_nodding_origin_frequency_inferred":True,
      "interval_limits":"pointwise t-model interval conditional on one tree, lambda, nugget and species mean; not a genus-wide phylogenetic hypothesis probability",
      "registered_panels":panels,
    }
    return {"summary":summary,"scenario_fits":fit_records}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source",required=True)
    ap.add_argument("--out",required=True)
    args=ap.parse_args()
    s1=verify_source(args.source)
    trees=restore_trees()
    report=compare(s1,trees)
    report["upstream_source_run_id"]=SOURCE_RUN
    report["upstream_source_artifact_id"]=SOURCE_ARTIFACT
    report["upstream_source_artifact_SHA"]=SOURCE_DIGEST
    report["original_continuous_S1_species_csv_source_sha256"]=hashlib.sha256(
                Path(args.source).read_bytes()).hexdigest()
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True,
                               allow_nan=False)+"\n",encoding="utf-8")
    panels=report["summary"]["registered_panels"]
    print("S1_FREE_LAMBDA_FIVE_CLIMATE_SUMMARY="+json.dumps({
       "status":report["summary"]["status"],
       "panels":{k:v for k,v in panels.items()
          if k.startswith("collision_free191|species_mean")
          and "nugget0.05" in k}
          | {k:v for k,v in panels.items()
          if k.startswith("collision_free_at_least10|species_mean|Pagel_lambda_MLE")
          and "nugget0.05" in k},
    },sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

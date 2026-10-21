"""Species-level rain association influence audit for ORIGINAL S1 image angles.

This uses source-hashed original 170,066 native S1 images aggregated onto the
same 194 tree-covered taxa, and prioritises the 191 collision-free taxa.
The five S1 covariates BIO1/BIO4/BIO12/BIO15/RSDS remain jointly adjusted.
All effects are image-axis degrees per original CHELSA increment; photo-derived
angle NEVER becomes a botanical nodding/erect label.

The inferential target here is taxon-equal five-covariate GLS at lambda=0
(favoured by the prior genuine-S1 52-tree lambda audit). This is an auxiliary
species-level robustness diagnostic, NOT S1's original observation-level
between/within, spatial/photographer, or frozen V3 analysis. HC3, taxon
leave-one-out and taxon-pairs bootstrap do not solve geographic confounding.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from analysis.v4.joint_s1_t2_phylogeny import EXPOSURES
from analysis.v4.s1_continuous_phylogenetic_signal import taxon_sets, verify_source

SOURCE_RUN=39594257539
SOURCE_ARTIFACT=12343166776
SOURCE_DIGEST="sha256:345617fd8fef6bdb41a06a1a03fa8a776c5ad48ba13c2aa971deb1e929cbf8cb"
SCALES=np.array([1.,100.,100.,10.,1.],dtype=float)
SELECTED_AXES=("chelsa_bio12","chelsa_bio15")
CI_NOMINAL=0.95
EXPECTED={"original194":194,"collision_free191":191,
          "collision_free_at_least5":171,"collision_free_at_least10":143}


def prepare(rows,field):
    if field not in {"image_mean_deg","image_median_deg"}:
        raise ValueError("Must use measured continuous S1 image angle")
    ordered=sorted(rows,key=lambda r:r["tip"])
    names=[r["tip"] for r in ordered]
    values=np.array([[float(r["S1_five_climate_means"][c]) for c in EXPOSURES]
                     for r in ordered],dtype=float)
    y=np.array([float(r[field]) for r in ordered],dtype=float)
    n=len(y)
    if n<50 or len(set(names))!=n or np.any(~np.isfinite(values)) or np.any(~np.isfinite(y)):
        raise ValueError("Unverified/source-collided or nonfinite S1 source cohort")
    X=np.column_stack([np.ones(n),(values-values.mean(axis=0))/SCALES])
    if np.linalg.matrix_rank(X)!=6 or np.linalg.cond(X)>1e5:
        raise ValueError("Climate predictors rank deficient or collinear")
    return names,X,y,values


def ols_hc3(X,y):
    n,p=X.shape
    if n<=p+2 or np.linalg.matrix_rank(X)!=p:
        raise ValueError("Insufficient independent taxon-level design")
    XtXi=np.linalg.inv(X.T@X)
    beta=XtXi@(X.T@y)
    residual=y-X@beta
    leverages=np.einsum("ni,ij,nj->n",X,XtXi,X)
    if np.any(leverages>=1-1e-8):
        raise ValueError("Extreme single-taxon leverage makes leave-one-out invalid")
    score=X * (residual/(1-leverages))[:,None]
    vcov=XtXi @ (score.T@score) @ XtXi
    se=np.sqrt(np.diag(vcov))
    if np.any(~np.isfinite(se)):
        raise ValueError("HC3 source uncertainty unavailable")
    crit=float(student_t.ppf(.975,n-p))
    return {
      "beta":beta,"se_hc3":se,
      "CI_HC3":np.stack([beta-crit*se,beta+crit*se],axis=1),
      "leverage":leverages,"residual":residual,
      "rank":p,"df":n-p,
      "max_abs_observed_angle":float(np.max(np.abs(y))),
    }


def audit_cohort(rows,field,*,reps=1999,seed=20261021):
    names,X,y,raw=prepare(rows,field)
    fit=ols_hc3(X,y)
    n=len(y); rng=np.random.default_rng(seed)
    source_axes={}
    for exposure in SELECTED_AXES:
        j=list(EXPOSURES).index(exposure)+1
        beta=float(fit["beta"][j])
        s=fit["se_hc3"][j]
        lower,upper=fit["CI_HC3"][j].tolist()
        # Recompute the original 5-covariate regression with every taxon
        # omitted. This is an exact influence analysis, not a null permutation.
        loo=[]
        loo_lower=[]
        loo_taxa=[]
        for leave in range(n):
            include=np.arange(n)!=leave
            # Translation of predictor centering does not change slopes.
            subfit=ols_hc3(X[include],y[include])
            loo.append(float(subfit["beta"][j]))
            loo_lower.append(float(subfit["CI_HC3"][j,0]))
            loo_taxa.append(names[leave])
        most_influential=int(np.argmax(np.abs(np.array(loo)-beta)))
        # Species-level pairs bootstrap; no within-taxon pseudoreplication.
        boot=[]
        draws_failed=0
        for i in range(reps):
            selected=rng.integers(0,n,size=n)
            try:
                betas=np.linalg.lstsq(X[selected],y[selected],rcond=None)[0]
                if np.linalg.matrix_rank(X[selected]) !=6:
                    raise ValueError("Rank-deficient bootstrap resample")
                boot.append(float(betas[j]))
            except (np.linalg.LinAlgError,ValueError):
                draws_failed+=1
        if len(boot)<max(100,reps*.95):
            raise ValueError("Species-level bootstrap rank failures excessive")
        values=np.array(boot)
        source_axes[exposure]={
           "source_original_5_joint_predictors":list(EXPOSURES),
           "increment_original_units":float(SCALES[j-1]),
           "estimate_image_degrees_per_increment":beta,
           "HC3_conditional_t95_interval":[float(lower),float(upper)],
           "robust_se_HC3":float(s),
           "species_leave_one_out_n":n,
           "species_leave_one_out_effect_range":[float(min(loo)),float(max(loo))],
           "species_leave_one_out_effect_all_positive":bool(all(z>0 for z in loo)),
           "species_leave_one_out_pointwise_HC3_lower95_positive_count":int(sum(v>0 for v in loo_lower)),
           "largest_single_taxon_effect_on_slope":{
                "taxon":loo_taxa[most_influential],
                "estimate_if_excluded":float(loo[most_influential]),
                "absolute_effect_change":float(abs(loo[most_influential]-beta)),
           },
           "species_pairs_bootstrap_requested":reps,
           "species_pairs_bootstrap_completed":int(len(boot)),
           "species_pairs_bootstrap_rank_failures":draws_failed,
           "species_pairs_bootstrap_percentile95":np.quantile(values,[.025,.975]).tolist(),
           "species_pairs_bootstrap_positive_fraction":float(np.mean(values>0)),
           "not_spatially_independent_sites_or_causal_evolution":True,
        }
    return {
      "n_original_S1_taxon_keys":n,
      "n_original_S1_images":int(sum(r["n_eligible_images"] for r in rows)),
      "phenotype":field,
      "five_S1_CHELSA_joint_predictors":list(EXPOSURES),
      "taxon_equal_weighting":True,
      "lambda_fixed_at_zero_due_to_source_checked_52_tree_signal_diagnostic":True,
      "not_a_replacement_for_original_observation_level_S1":True,
      "max_hat_leverage":float(np.max(fit["leverage"])),
      "n_hat_leverage_greater_0_2":int(np.sum(fit["leverage"]>0.2)),
      "predictor_rank":fit["rank"],
      "rain_exposure_robustness":source_axes,
    }


def run(source,out,*,bootstrap=1999):
    rows=verify_source(source)
    cohorts=taxon_sets(rows)
    if {k:len(v) for k,v in cohorts.items()}!=EXPECTED:
        raise ValueError("S1 data source/subset count drift")
    if bootstrap<100 or bootstrap>20000:
        raise ValueError("Bootstrap budget unsupported")
    result={}
    for i,(cohort,items) in enumerate(cohorts.items()):
        result[cohort]={
            field:audit_cohort(items,field,reps=bootstrap,seed=20261021+i*2+offset)
            for offset,field in enumerate(("image_mean_deg","image_median_deg"))}
    receipt={
       "status":"EMPIRICAL_S1_ANNUAL_RAINFALL_CONTINUOUS_ANGLE_SINGLE_TAXON_INFLUENCE_AUDIT",
       "registered_climate_vars":list(EXPOSURES),
       "primary":"chelsa_bio12",
       "secondary":"chelsa_bio15",
       "strict_main":"collision_free191",
       "n_194_all_tip_original_S1_observations":169912,
       "original_full_native_S1_image_cohort":170066,
       "original_full_native_S1_accepted_keys":244,
       "source_original_S1_run_id":SOURCE_RUN,
       "source_original_S1_artifact_id":SOURCE_ARTIFACT,
       "source_original_S1_artifact_digest":SOURCE_DIGEST,
       "S1_source_file_sha256":hashlib.sha256(Path(source).read_bytes()).hexdigest(),
       "same_5_climates_and_original_image_degrees":True,
       "all_52_grafted_phylogenetic_scenarios_not_independent_replication":True,
       "λ0_justification":"S1 clean191 five-climate lambda MLE median 0/52, source-original NA independent nuclear42 lambda MLE 0; this is an independent taxon-level diagnostic not phylogenetic causal proof",
       "original_20261021_S1_between_within_model_not_reestimated":True,
       "not_a_history_of_nodding_gains":True,
       "not_evidence_of_precipitation_evolutionary_causality":True,
       "cohort_results":result,
    }
    path=Path(out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n",
                    encoding="utf-8")
    selected={}
    for key in ("collision_free191","collision_free_at_least10"):
        s=result[key]["image_mean_deg"]["rain_exposure_robustness"]
        selected[key]={e:s[e] for e in SELECTED_AXES}
    print("S1_ANNUAL_RAIN_LOO_HC3_SUMMARY="+
          json.dumps({"status":receipt["status"],"rain_audit":selected},
                     ensure_ascii=False,sort_keys=True),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--bootstrap",type=int,default=1999)
    args=p.parse_args()
    run(args.source,args.out,bootstrap=args.bootstrap)


if __name__=="__main__":
    main()

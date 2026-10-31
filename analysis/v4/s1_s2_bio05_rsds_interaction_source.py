"""S1/S2 source-verified elevation vs BIO05×RSDS interaction audit.

This is a *post-result exploratory* comparison. Do NOT select a model by which
correlations you would like to retain. Same original all-12-environment-complete
photographs/weights are used within every set of four models:
  M3: BIO18 + BIO05 + RSDS
  M4: M3 + elevation (the archived V4 baseline)
  M3I: M3 + BIO05×RSDS
  M4I: M4 + BIO05×RSDS

S1 moderation is explicitly BETWEEN taxon climate means; no within-taxon
interaction is inferred. S2 moderation is within the fitted species×range
stratum. Interaction tests are Holm-adjusted separately across all 8
phenotypic constructs (S1) and all 8×4 strata (S2).
Original frozen S1/S2 coefficients and Holm families are never rewritten.
No observation identifiers, taxon keys, coordinates or uploader IDs exported.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_integration import FEATURES
from analysis.v4.s1_s2_extreme_lego import ALL_ENV, read_tables, regression_design
from analysis.v4.s1_s2_extreme_range_models import _get_dated_balanced, _formula
from analysis.v4.s1_s2_vif5_si10_alltraits import (
    standardized_y, hue_joint, selected_source as old_source
)
from analysis.v4.s1_s2_climate_only_3_vs_5 import (
    source_checks, source_s1, climate_only_source_mask,
    fit_regression, coefficient, no_private_keys
)

CLIMATE=("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean")
TERRAIN=(*CLIMATE,"topo_elevation")
INTERACTION="BIO05_x_RSDS"
MODELS={
    "M3_climate":(CLIMATE,False),
    "M4_original":(TERRAIN,False),
    "M3_plus_interaction":(CLIMATE,True),
    "M4_plus_interaction":(TERRAIN,True),
}
EXPECTED_S1_ANGLE_N=168907
EXPECTED_S1_ANGLE_TAXA=178
EXPECTED_OLD_BIO18_BETA=0.04544684307686103


def append_weighted_interaction(X, weights, *, mode):
    """Append a centred BIO05×RSDS product, keeping both main effects.

    Both regression designs have column 2 = temperature, column 3 = RSDS.
    S1 uses BETWEEN-taxon means at these locations. S2 uses their standardised
    climate observations within one species-range stratum.
    """
    if mode not in ("between_taxa","within_species_range"):
        raise ValueError("UNDECLARED_INTERACTION_MODE")
    if X.ndim != 2 or X.shape[1] < 4 or len(X)!=len(weights):
        raise ValueError("UNEXPECTED_INTERACTION_DESIGN")
    w=np.asarray(weights,float)
    if (not np.isfinite(X).all() or not np.isfinite(w).all() or
        (w<=0).any()):
        raise ValueError("INVALID_INTERACTION_INPUT")
    product=X[:,2]*X[:,3]
    mean=np.average(product,weights=w)
    standard=np.sqrt(np.average((product-mean)**2,weights=w))
    if not np.isfinite(standard) or standard < 1.e-12:
        raise ValueError("ZERO_INTERACTION_SPREAD")
    z=(product-mean)/standard
    expanded=np.column_stack((X,z))
    if np.linalg.matrix_rank(expanded)<expanded.shape[1]:
        raise ValueError("INTERACTION_NOT_IDENTIFIABLE")
    return expanded,{
        "interaction_mode":mode,
        "interacting_main_effects":["chelsa_bio05","chelsa_rsds_mean"],
        "weighted_product_mean":float(mean),
        "weighted_product_sd":float(standard),
        "product_standardized":True,
        "all_main_effects_retained":True,
    }


def fit_one_cohort(data, outcomes, variables, interaction, *, scope):
    if scope=="S1_between_taxa":
        X,groups=regression_design(data,list(variables))
        n=data.groupby("accepted_key").size()
        weights=data.accepted_key.map(1/n).to_numpy(float)
        mode="between_taxa"
    elif scope=="S2_within_range":
        X,_,_,kept,removed=_formula(data,list(variables))
        groups=data["_cell"].to_numpy()
        n=data.groupby("_cell").size()
        weights=data._cell.map(1/n).to_numpy(float)
        mode="within_species_range"
    else:
        raise ValueError("INVALID_COHORT_SCOPE")
    dims=X.shape[1]
    diag={"interaction_added":bool(interaction)}
    if scope=="S2_within_range":
        diag["date_nuisance_retained"]=list(kept)
        diag["date_nuisance_omitted"]=list(removed)
    if interaction:
        X,diag_extra=append_weighted_interaction(X,weights,mode=mode)
        diag.update(diag_extra)
    fits={
        name:fit_regression(
            X,standardized_y(data[name].to_numpy(float)),weights,groups
        ) for name in outcomes
    }
    aic={name:float(fit.aic) for name,fit in fits.items()}
    sse={name:float(np.sum(weights*fit.resid**2)) for name,fit in fits.items()}
    slopes={}
    for name,fit in fits.items():
        if name.startswith("corolla_hue_"):
            continue
        scales=(
            [("between_equal_taxa",1),("within_taxon",len(variables)+1)]
            if scope=="S1_between_taxa" else [("within_species_range",1)]
        )
        slopes[name]={
            scope_name:{
                env:coefficient(fit,offset+i)
                for i,env in enumerate(variables)
            } for scope_name,offset in scales
        }
        if interaction:
            slopes[name]["interaction"]={INTERACTION:coefficient(fit,X.shape[1]-1)}
    if interaction:
        rows=[]
        if list(outcomes)==list(FEATURES["hue"]):
            hue=hue_joint(
                [fits[outcomes[0]],fits[outcomes[1]]],
                X,weights,groups,[INTERACTION],
                [("between_equal_taxa" if scope=="S1_between_taxa" else "within_species_range",X.shape[1]-1)]
            )
            if len(hue)!=1 or hue[0]["status"]!="JOINT_CIRCULAR_HUE_WALD_2DF":
                raise ValueError("CIRCULAR_HUE_INTERACTION_UNIDENTIFIABLE")
            rows=[{
                "trait":"circular_hue_joint","p_uncorrected":float(hue[0]["p_uncorrected"]),
                "chi_square_2df":float(hue[0]["wald_chi2"]),
                "beta_std":None,"ci95_pointwise":None,
                "scale":hue[0]["scale"],"kind":"joint_circular_hue_2df",
            }]
        else:
            rows=[{
                "trait":name,"p_uncorrected":float(fit.pvalues[-1]),
                "beta_std":float(fit.params[-1]),
                "ci95_pointwise":[float(v) for v in np.asarray(fit.conf_int())[-1]],
                "chi_square_2df":None,
                "scale":"between_equal_taxa" if scope=="S1_between_taxa" else "within_species_range",
                "kind":"scalar_cluster_robust",
            } for name,fit in fits.items()]
    else:
        rows=[]
    return {
        "n_photos":int(len(data)),
        "n_clusters":int(pd.Series(groups).nunique()),
        "predictors":list(variables),
        "n_design_columns_before_interaction":dims,
        "model":diag,"AIC":aic,"weighted_SSE":sse,
        "main_effect_slopes":slopes,
        "interaction_test_rows":rows,
    }


def summarize_nested(a,b):
    """Same-cohort in-sample fit comparison; no causal or CV claim."""
    keys=set(a["AIC"])&set(b["AIC"])
    if len(keys)!=len(a["AIC"]) or len(keys)!=len(b["AIC"]):
        raise ValueError("NESTED_OUTCOME_KEYS_CHANGED")
    if a["n_photos"]!=b["n_photos"] or a["n_clusters"]!=b["n_clusters"]:
        raise ValueError("SAMPLE_OR_CLUSTER_CHANGED_IN_COMPARISON")
    return {
        k:{
            "AIC_new_minus_reference":float(b["AIC"][k]-a["AIC"][k]),
            "weighted_residual_fraction_removed":float(
                1-b["weighted_SSE"][k]/a["weighted_SSE"][k]),
            "same_photos_weights_and_clusters":True,
            "AIC_in_sample_not_independent_validation":True,
        } for k in sorted(keys)
    }


def correct_interactions(rows,scope,model_name):
    if not rows or any(not np.isfinite(r["p_uncorrected"]) for r in rows):
        raise ValueError("INVALID_INTERACTION_TEST_FAMILY")
    p=np.array([r["p_uncorrected"] for r in rows])
    holm=multipletests(p,method="holm")[1]
    for r,v in zip(rows,holm):
        r["holm_p"]=float(v)
        r["holm_005"]=bool(v<.05)
    return {
        "scope":scope,"model":model_name,
        "family":"one predetermined BIO05×RSDS interaction per trait and stratum",
        "n_tests":int(len(rows)),
        "n_holm":int(np.sum(holm<.05)),
        "rows":rows,
        "post_result_exploratory":True,
        "not_original_october_Holm_family":True,
    }


def run_s1(native):
    # a single old ALL12-complete cohort for ALL models; new 3-variable
    # unrelated-photo set is separately recorded only for sample audit.
    cohorts={}
    interactions={"M3_plus_interaction":[],"M4_plus_interaction":[]}
    models={}
    for block,columns in FEATURES.items():
        original,info=old_source(native,columns,min_taxa=20)
        if original is None:
            raise ValueError("S1_ORIGINAL_EIGHT_TRAIT_QC_MISSING")
        group={
            "n_photos":int(len(original)),
            "n_taxa":int(original.accepted_key.nunique()),
        }
        if block=="angle":
            if group != {"n_photos":EXPECTED_S1_ANGLE_N,"n_taxa":EXPECTED_S1_ANGLE_TAXA}:
                raise ValueError("S1_ORIGINAL_ANGLE_COHORT_DRIFT")
        cohorts[block]=group
        fresh, fresh_status=source_s1(native,columns)
        if fresh is None or len(fresh)<len(original):
            raise ValueError("S1_FRESH_CLIMATE_PHOTOS_INSUFFICIENT")
        cohorts[block]["fresh6_photos"]=int(len(fresh))
        cohorts[block]["photos_added_without_elevation_gate"]=int(len(fresh)-len(original))
        fresh3=fit_one_cohort(fresh,columns,CLIMATE,False,scope="S1_between_taxa")
        specout={
            name:fit_one_cohort(original,columns,vars_,interactive,scope="S1_between_taxa")
            for name,(vars_,interactive) in MODELS.items()
        }
        if block=="angle":
            old=specout["M4_original"]["main_effect_slopes"][columns[0]]["between_equal_taxa"]["chelsa_bio18"]["beta_std"]
            if not np.isclose(old,EXPECTED_OLD_BIO18_BETA,atol=3e-7,rtol=0):
                raise ValueError("FROZEN_BIO18_BETA_NOT_REPRODUCED")
        for name in interactions:
            for row in specout[name]["interaction_test_rows"]:
                interactions[name].append({"block":block,**row,"n_photos":len(original)})
        specout["fresh6_climate3_unmatched_sample_diagnostic"]=fresh3
        specout["contrasts"]={
            "add_elevation_to_three":summarize_nested(specout["M3_climate"],specout["M4_original"]),
            "add_interaction_to_three":summarize_nested(specout["M3_climate"],specout["M3_plus_interaction"]),
            "add_interaction_to_original_four":summarize_nested(specout["M4_original"],specout["M4_plus_interaction"]),
            "three_interaction_vs_original_four":summarize_nested(specout["M4_original"],specout["M3_plus_interaction"]),
        }
        models[block]=specout
    tests={name:correct_interactions(rows,"S1_between_taxa",name)
           for name,rows in interactions.items()}
    if any(r["n_tests"]!=8 for r in tests.values()):
        raise ValueError("S1_EIGHT_CONSTRUCTS_NOT_COVERED")
    return {"cohorts":cohorts,"models":models,"interaction_families":tests,
            "all_original_12_env_rows_fixed":True,
            "interactions_between_taxon_only":True}


def run_s2(focal):
    cohorts={}
    tests={"M3_plus_interaction":[],"M4_plus_interaction":[]}
    models={}
    for (species,state),group in focal.groupby(["taxon_name","range_state"],sort=True):
        original=group.loc[np.isfinite(group[ALL_ENV].to_numpy(float)).all(axis=1)].copy()
        label=species+" / "+state
        if original.empty:
            raise ValueError("S2_ALL12_ORIGINAL_COHORT_MISSING")
        blocks={}
        for block,columns in FEATURES.items():
            if not set(columns)<=set(original):
                raise ValueError("S2_PHOTOGRAPHIC_TRAIT_MISSING")
            mask=np.isfinite(original[columns].to_numpy(float)).all(axis=1)
            sample=original.loc[mask].copy()
            if "orientation_image_vertical_angle" in columns:
                sample=sample.loc[sample.orientation_image_vertical_angle.between(0,180)].copy()
            sample=_get_dated_balanced(sample)
            if len(sample)<100 or sample._cell.nunique()<15 or sample.user_id.nunique()<10:
                raise ValueError("S2_ORIGINAL_COHORT_INSUFFICIENT")
            fits={
                name:fit_one_cohort(sample,columns,variables,interactive,
                                    scope="S2_within_range")
                for name,(variables,interactive) in MODELS.items()
            }
            for name in tests:
                for row in fits[name]["interaction_test_rows"]:
                    tests[name].append({
                        "species":species,"range":state,"block":block,
                        "n_photos":len(sample),"n_cells":int(sample._cell.nunique()),
                        **row
                    })
            fits["contrasts"]={
                "add_elevation_to_three":summarize_nested(fits["M3_climate"],fits["M4_original"]),
                "add_interaction_to_three":summarize_nested(fits["M3_climate"],fits["M3_plus_interaction"]),
                "add_interaction_to_original_four":summarize_nested(fits["M4_original"],fits["M4_plus_interaction"]),
                "three_interaction_vs_original_four":summarize_nested(fits["M4_original"],fits["M3_plus_interaction"]),
            }
            # month/year dummies are also kept on the SAME photos. The time
            # terms selected by rank should be documented for all fits.
            fresh_basis=group.loc[climate_only_source_mask(group)].copy()
            fresh_trait=fresh_basis.loc[
                np.isfinite(fresh_basis[columns].to_numpy(float)).all(axis=1)].copy()
            if "orientation_image_vertical_angle" in columns:
                fresh_trait=fresh_trait.loc[
                    fresh_trait.orientation_image_vertical_angle.between(0,180)].copy()
            fresh_sample=_get_dated_balanced(fresh_trait)
            if len(fresh_sample)<len(sample):
                raise ValueError("S2_FRESH_PHOTO_ELIGIBILITY_SMALLER")
            fits["fresh6_climate3_unmatched_sample_diagnostic"]=fit_one_cohort(
                fresh_sample,columns,CLIMATE,False,scope="S2_within_range")
            fits["fresh_sample_diagnostic"]={
                "fresh6_photos":int(len(fresh_sample)),
                "old12_photos":int(len(sample)),
                "additional_photos_without_elevation_gate":int(len(fresh_sample)-len(sample)),
                "comparisons_with_fresh6_not_same_photos":True,
                "AIC_between_old12_and_fresh6_not_comparable":True,
            }
            blocks[block]=fits
        cohorts[label]={"wild_source_observations":int(len(group)),
                        "all12_complete_before_photo_QC":int(len(original))}
        models[label]=blocks
    if len(cohorts)!=4:
        raise ValueError("FOUR_S2_STRATA_MISSING")
    families={name:correct_interactions(rows,"S2_four_species_range_strata",name)
              for name,rows in tests.items()}
    if any(f["n_tests"]!=32 for f in families.values()):
        raise ValueError("S2_32_INTERACTION_TESTS_NOT_COVERED")
    return {"cohorts":cohorts,"models":models,"interaction_families":families,
            "same_original_all12_climate_completeness_for_all_four_models":True}


def run(native, focal, extra, measurements, output):
    s1,s2=read_tables(native,focal,extra,measurements)
    original=source_checks(s1,s2)
    result={
        "schema":"azami_v4_s1_s2_BIO05_RSDS_interaction_vs_elevation_source_audit_v1",
        "status":"EXPLORATORY_SOURCE_VERIFIED_IF_COMPLETED",
        "design":{
            "model_set":{name:{"predictors":list(vars_),"has_BIO05_RSDS_interaction":flag}
                        for name,(vars_,flag) in MODELS.items()},
            "old12_source_complete_cohorts_held_constant":True,
            "main_four_covariates_historical":"BIO18/BIO05/RSDS/elevation",
            "S1_interaction_scale":"between-taxon only",
            "S2_interaction_scale":"within species-range stratum",
            "study_3_3_original_BIO12_BIO01_RSDS_untouched":True,
            "no_significance_driven_selection_or_discarding_unwanted_results":True,
            "old_frozen_Holm_80_S1_and_128_S2_untouched":True,
        },
        "source_checks":original,
        "S1":run_s1(s1),
        "S2":run_s2(s2),
    }
    no_private_keys(result)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,allow_nan=False,indent=2)+"\n",
                      encoding="utf-8")
    return result


def main():
    p=argparse.ArgumentParser()
    for flag in ("native","focal","extra","measurements","output"):
        p.add_argument("--"+flag,required=True,type=Path)
    args=p.parse_args()
    result=run(args.native,args.focal,args.extra,args.measurements,args.output)
    print(json.dumps({
        "status":result["status"],
        "S1_interaction_tests":{
            k:{"tests":v["n_tests"],"holm_positive":v["n_holm"]}
            for k,v in result["S1"]["interaction_families"].items()},
        "S2_interaction_tests":{
            k:{"tests":v["n_tests"],"holm_positive":v["n_holm"]}
            for k,v in result["S2"]["interaction_families"].items()},
        "public_output":"ANONYMISED_AGGREGATES_ONLY",
    }))


if __name__=="__main__":
    main()

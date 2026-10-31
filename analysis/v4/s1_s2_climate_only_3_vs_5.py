"""Source-gated climate-only 3/5-predictor reanalysis of V4 S1 and S2.

New results 3.1/3.2: use a common climate-only completeness gate that does
NOT condition on elevation; three predeclared climate specifications share
each outcome's photographs and original sampling weights. Fully redo BH/Holm
for the eight actual photographic constructs, with circular hue as one 2df
biological endpoint. Preserve the original four-predictor results untouched.

Source photographs, uploader IDs, source taxon keys and precise coordinates
are private to the authorised GitHub Actions runner; export aggregates only.
"""
from __future__ import annotations

import argparse
import json
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from analysis.v4.native_integration import FEATURES
from analysis.v4.s1_s2_extreme_lego import (
    ALL_ENV, read_tables, regression_design, require, diagnostic_matrix,
)
from analysis.v4.s1_s2_vif5_si10_alltraits import (
    standardized_y, hue_joint,
    selected_source as frozen_s1_sample,
)
from analysis.v4.s1_s2_extreme_range_models import _get_dated_balanced, _formula
from analysis.v4.s2_original_colour_spatial_block_sensitivity import (
    EXPECTED as ORIGINAL_S2_COLOUR,
    original_sample, _original_fit,
)

MODELS = OrderedDict((
    ("annual_temperature_three", (
        "chelsa_bio18", "chelsa_bio01", "chelsa_rsds_mean")),
    ("warmest_temperature_three", (
        "chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean")),
    ("seasonal_five", (
        "chelsa_bio18", "chelsa_bio01", "chelsa_rsds_mean",
        "chelsa_bio04", "chelsa_bio15")),
))
PRIMARY = "annual_temperature_three"
HOT = "warmest_temperature_three"
FIVE = "seasonal_five"
CLIMATE_REQUIRED = tuple(dict.fromkeys(
    term for vars_ in MODELS.values() for term in vars_))
ANGLE = "orientation_image_vertical_angle"
ORIGINAL_S1_N = 168907
ORIGINAL_S1_TAXA = 178
ORIGINAL_S1_BETA = 0.04544684307686103
MIN_TAXON_N = 5
MIN_TAXA_S1 = 20
MIN_PHOTOS_S2 = 100
MIN_CELLS_S2 = 15
MIN_OBSERVERS_S2 = 10
FORBIDDEN_OUTPUT = {
    "obs_id", "photo_id", "user_id", "photographer_id",
    "accepted_key", "taxon_name", "latitude", "longitude",
    "individual", "sample_id",
}


def source_checks(native: pd.DataFrame, focal: pd.DataFrame):
    """Reproduce old source results in private memory; report only pass flags."""
    frozen, _ = frozen_s1_sample(native, [ANGLE], min_taxa=20)
    if frozen is None or len(frozen) != ORIGINAL_S1_N or frozen.accepted_key.nunique() != ORIGINAL_S1_TAXA:
        raise ValueError("FROZEN_S1_ANGLE_ELIGIBILITY_FAILED")
    X, taxon = regression_design(frozen, list((
        "chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean", "topo_elevation")))
    group_n = frozen.groupby("accepted_key").size()
    w = frozen.accepted_key.map(1/group_n).to_numpy(float)
    y = standardized_y(frozen[ANGLE].to_numpy(float))
    baseline = sm.WLS(y, X, weights=w).fit(
        cov_type="cluster", cov_kwds={"groups": taxon}, use_t=False)
    if not np.isclose(baseline.params[1], ORIGINAL_S1_BETA, atol=3e-7, rtol=0):
        raise ValueError("FROZEN_S1_COEFFICIENT_FAILED")
    checked=[]
    for species, state, trait, env, required_n, beta in ORIGINAL_S2_COLOUR:
        sample, _ = original_sample(focal, species, state, trait)
        if len(sample) != required_n:
            raise ValueError("FROZEN_S2_COLOUR_COHORT_FAILED")
        *_, result = _original_fit(sample, trait, env)
        if not np.isclose(result["beta"], beta, atol=3e-7, rtol=0):
            raise ValueError("FROZEN_S2_COEFFICIENT_FAILED")
        checked.append({"species":species, "range":state, "n":required_n,
                        "outcome":trait, "old_beta":float(beta), "verified":True})
    return {
        "original_S1_168907_178_and_BIO18_beta_verified":True,
        "original_S2_six_colour_coefficients_verified":True,
        "old_S1_BIO18_standardized_beta":float(baseline.params[1]),
        "old_S2_provenance":checked,
        "original_four_environment_variables_only_used_to_verify_source":True,
    }


def climate_only_source_mask(frame:pd.DataFrame):
    return np.isfinite(frame[list(CLIMATE_REQUIRED)].to_numpy(float)).all(axis=1)


def source_s1(frame:pd.DataFrame, members:list[str]):
    """Unlike frozen selection, do NOT filter on absent elevation."""
    if not all(trait in frame and trait+"__eligible" in frame for trait in members):
        return None, {"status":"PHOTO_MEASURE_OR_QC_COLUMN_MISSING"}
    mask=climate_only_source_mask(frame)
    for t in members:
        mask &= np.isfinite(pd.to_numeric(frame[t],errors="coerce").to_numpy(float))
        mask &= frame[t+"__eligible"].fillna(False).astype(bool).to_numpy()
    if ANGLE in members:
        mask &= frame[ANGLE].between(0,180).fillna(False).to_numpy(bool)
    source=frame.loc[mask].copy()
    count=source.groupby("accepted_key").size()
    source=source.loc[source.accepted_key.isin(count[count>=MIN_TAXON_N].index)].copy()
    if source.accepted_key.nunique()<MIN_TAXA_S1:
        return None,{"status":"INSUFFICIENT_CLIMATE_ONLY_TAXON_SUPPORT"}
    return source,{"status":"CLIMATE_ONLY_SOURCE_LINKED",
                   "n":int(len(source)), "taxa":int(source.accepted_key.nunique())}


def fit_regression(X,y,weights,groups):
    fit=sm.WLS(y,X,weights=weights).fit(
        cov_type="cluster",cov_kwds={"groups":groups},use_t=False)
    if not (np.isfinite(fit.params).all() and np.isfinite(fit.bse).all()):
        raise ValueError("NONFINITE_CLIMATE_ONLY_REGRESSION")
    return fit


def coefficient(fit, j):
    return {
        "beta_std":float(fit.params[j]),
        "se_cluster":float(fit.bse[j]),
        "p_uncorrected":float(fit.pvalues[j]),
        "ci95_pointwise": [float(q) for q in np.asarray(fit.conf_int())[j]],
    }


def vif(frame,variables):
    d=diagnostic_matrix(frame[list(variables)],list(variables))
    return {
        "status":d["status"],
        "max_vif":None if d["status"]!="OK" else float(max(d["vif"].values())),
        "max_abs_pair_r":d.get("max_abs_pair_r"),
        "vif":d.get("vif"),
        "vif_gt5_flag":bool(d["status"]=="OK" and max(d["vif"].values())>5),
    }


def append_fitted_rows(rows, fit_map, variables, cohort_label, sample_n,
                       *, scale_info, taxon_or_cell_n, weighting):
    for name,fit in fit_map.items():
        if name in ("corolla_hue_sin","corolla_hue_cos"):
            continue
        for scale_name,base in scale_info:
            for k,env in enumerate(variables):
                rows.append({
                    "model":cohort_label, "trait":name, "scale":scale_name,
                    "environment":env, "n_photos":sample_n,
                    "groups":taxon_or_cell_n, "weighting":weighting,
                    "test_kind":"single_continuous_construct",
                    **coefficient(fit,base+k),
                })


def append_hue_rows(rows, fits, X, w, group, variables,
                    cohort_label, n, scales, cluster_n, weighting):
    joint=hue_joint(
        [fits["corolla_hue_sin"],fits["corolla_hue_cos"]],
        X,w,group,list(variables),scales,
    )
    if any(r["status"]!="JOINT_CIRCULAR_HUE_WALD_2DF" for r in joint):
        raise ValueError("CIRCULAR_HUE_UNIDENTIFIABLE_IN_CLIMATE_FAMILY")
    for r in joint:
        rows.append({
            "model":cohort_label, "trait":"circular_hue_joint",
            "scale":r["scale"], "environment":r["environment"],
            "n_photos":n,"groups":cluster_n,"weighting":weighting,
            "test_kind":"joint_sine_cosine_2df",
            "beta_std":None,"se_cluster":None,
            "p_uncorrected":float(r["p_uncorrected"]),
            "ci95_pointwise":None,
            "joint_chisq_2df":float(r["wald_chi2"]),
        })


def run_s1(native):
    cohorts={}
    tests={k:[] for k in MODELS}
    for block, columns in FEATURES.items():
        part,eligibility=source_s1(native,columns)
        if part is None:
            raise ValueError("S1_CLIMATE_ONLY_TRAIT_INELIGIBLE_"+block)
        taxon=part.accepted_key.astype(str).to_numpy()
        counts=part.groupby("accepted_key").size()
        weight=part.accepted_key.map(1/counts).to_numpy(float)
        y={col:standardized_y(part[col].to_numpy(float)) for col in columns}
        cohort={"n":len(part),"taxa":int(part.accepted_key.nunique()),
                "variables_used_for_completeness":list(CLIMATE_REQUIRED),
                "models":{}}
        for spec,variables in MODELS.items():
            X, group=regression_design(part,list(variables))
            if not np.array_equal(group,taxon):
                raise ValueError("S1_GROUPS_CHANGED")
            fits={name:fit_regression(X,y[name],weight,group) for name in columns}
            means={
                name:{
                    "AIC":float(fit.aic),
                    "weighted_SSE":float(np.sum(weight*np.asarray(fit.resid)**2)),
                } for name,fit in fits.items()
            }
            cohort["models"][spec]={"predictors":list(variables),
                                    "VIF_on_identical_phenotype_cohort":vif(part,variables),
                                    "diagnostic_fit":means}
            if block=="hue":
                append_hue_rows(tests[spec],fits,X,weight,group,variables,
                    spec,len(part),
                    [("between_equal_taxa",1),("within_taxon",1+len(variables))],
                    int(part.accepted_key.nunique()),"taxon_equal_weight")
            else:
                append_fitted_rows(tests[spec],fits,variables,spec,len(part),
                    scale_info=[("between_equal_taxa",1),("within_taxon",1+len(variables))],
                    taxon_or_cell_n=int(part.accepted_key.nunique()),
                    weighting="taxon_equal_weight")
        for name in columns:
            three=cohort["models"][PRIMARY]["diagnostic_fit"][name]
            five=cohort["models"][FIVE]["diagnostic_fit"][name]
            cohort["models"][FIVE]["diagnostic_fit"][name]["AIC_five_minus_three_annual_temp"] = float(five["AIC"]-three["AIC"])
        cohorts[block]=cohort
    return {"status":"CLIMATE_ONLY_S1_SOURCE_FITTED",
            "traits":cohorts,
            "families":apply_multiplicity(tests,"S1"),
            "source_observation_total":len(native)}


def fit_s2_part(part,species,state,block,columns):
    sample=_get_dated_balanced(part)
    if (len(sample)<MIN_PHOTOS_S2
        or sample._cell.nunique()<MIN_CELLS_S2
        or sample.user_id.nunique()<MIN_OBSERVERS_S2):
        raise ValueError("S2_CLIMATE_ONLY_SPATIAL_SUPPORT_INSUFFICIENT")
    group=sample._cell.to_numpy()
    sizes=sample.groupby("_cell").size()
    weight=sample._cell.map(1/sizes).to_numpy(float)
    y={name:standardized_y(sample[name].to_numpy(float)) for name in columns}
    models={}
    tests={k:[] for k in MODELS}
    for key,variables in MODELS.items():
        X,_,_,retained,omitted=_formula(sample,list(variables))
        fits={name:fit_regression(X,y[name],weight,group) for name in columns}
        models[key]={
            "predictors":list(variables),
            "VIF_on_identical_phenotype_cohort":vif(sample,variables),
            "nuisance_month_year_retained":list(retained),
            "nuisance_month_year_omitted":list(omitted),
            "diagnostic_fit":{
                name:{"AIC":float(fit.aic),
                      "weighted_SSE":float(np.sum(weight*np.asarray(fit.resid)**2))}
                for name,fit in fits.items()
            },
        }
        if block=="hue":
            append_hue_rows(tests[key],fits,X,weight,group,variables,key,len(sample),
                [("within_species_range",1)],
                int(sample._cell.nunique()),"equal_5degree_cell_weight")
        else:
            append_fitted_rows(tests[key],fits,variables,key,len(sample),
                scale_info=[("within_species_range",1)],
                taxon_or_cell_n=int(sample._cell.nunique()),
                weighting="equal_5degree_cell_weight")
    comparable=(models[PRIMARY]["nuisance_month_year_retained"]==
                models[FIVE]["nuisance_month_year_retained"])
    for name in columns:
        if comparable:
            a=models[PRIMARY]["diagnostic_fit"][name]
            b=models[FIVE]["diagnostic_fit"][name]
            b["AIC_five_minus_three_annual_temp"]=float(b["AIC"]-a["AIC"])
        else:
            models[FIVE]["diagnostic_fit"][name]["AIC_five_minus_three_annual_temp"]=None
    return {
        "n":int(len(sample)), "cells":int(sample._cell.nunique()),
        "observers":int(sample.user_id.nunique()),
        "same_photographs_and_weights_three_vs_five":True,
        "same_date_nuisance_three_vs_five":bool(comparable),
        "models":models,
    },tests


def run_s2(focal):
    blocks={}
    table={k:[] for k in MODELS}
    for (species,state),g in focal.groupby(["taxon_name","range_state"],sort=True):
        complete=climate_only_source_mask(g)
        basis=g.loc[complete].copy()
        if basis.empty:
            raise ValueError("S2_CLIMATE_ONLY_ENVIRONMENT_NOT_AVAILABLE")
        group_label=species+"/"+state
        results={"original_wild_observations":int(len(g)),
                 "climate_six_complete_original_observations":int(len(basis)),
                 "blocks":{}}
        for block,columns in FEATURES.items():
            if not set(columns)<=set(basis):
                raise ValueError("S2_ORIGINAL_FEATURE_MISSING_"+block)
            p=basis.loc[np.isfinite(basis[columns].to_numpy(float)).all(axis=1)].copy()
            if ANGLE in columns:
                p=p.loc[p[ANGLE].between(0,180)].copy()
            entry,entries=fit_s2_part(p,species,state,block,columns)
            results["blocks"][block]=entry
            for k in MODELS:
                for row in entries[k]:
                    row["species"]=species
                    row["range"]=state
                    table[k].append(row)
        blocks[group_label]=results
    if len(blocks)!=4:
        raise ValueError("S2_FOUR_RANGE_STRATA_CHANGED")
    return {"status":"CLIMATE_ONLY_S2_SOURCE_FITTED","source_original_observations":len(focal),
            "strata":blocks, "families":apply_multiplicity(table,"S2")}


def apply_multiplicity(families,scope):
    out={}
    for name,rows in families.items():
        p=np.array([row["p_uncorrected"] for row in rows],float)
        if not np.isfinite(p).all() or len(p)==0:
            raise ValueError("NEW_CLIMATE_ONLY_FAMILY_HAS_NONFINITE_P")
        bh=multipletests(p,method="fdr_bh")[1]
        holm=multipletests(p,method="holm")[1]
        for row,q1,q2 in zip(rows,bh,holm):
            row["p_BH"]=float(q1)
            row["p_Holm"]=float(q2)
            row["Holm_lt_0_05"]=bool(q2<.05)
        out[name]={
            "status":"NEW_CLIMATE_ONLY_EXPLORATORY_MULTIPLICITY",
            "scope":scope,
            "n_tests":len(rows),
            "n_holm_005":sum(row["Holm_lt_0_05"] for row in rows),
            "n_bh_005":int((bh<.05).sum()),
            "n_nominal_005":int((p<.05).sum()),
            "results":rows,
            "not_original_20261024_Holm_family":True,
            "no_after_result_predictor_selection":True,
        }
    return out


def no_private_keys(value):
    def check(node):
        if isinstance(node,dict):
            for key,val in node.items():
                if key in FORBIDDEN_OUTPUT:
                    raise ValueError("PRIVATE_SOURCE_FIELD_IN_AGGREGATE:"+key)
                check(val)
        elif isinstance(node,list):
            for val in node: check(val)
    check(value)


def execute(native_path,focal_path,extra_path,measurements_path,output):
    native,focal=read_tables(native_path,focal_path,extra_path,measurements_path)
    require(set(CLIMATE_REQUIRED)<=set(native) and set(CLIMATE_REQUIRED)<=set(focal),
            "CLIMATE_ONLY_INPUT_NOT_COMPLETE")
    verification=source_checks(native,focal)
    out={
        "schema":"V4_S1_S2_CLIMATE_ONLY_3_VS_5_20261031",
        "status":"SOURCE_VERIFIED_NEW_POST_FREEZE_CLIMATE_ONLY_RESULTS",
        "design":{
            "models":{k:list(v) for k,v in MODELS.items()},
            "all_models_exclude_elevation":True,
            "climate_only_complete_row_eligibility_excludes_elevation":True,
            "same_phenotype_QC_cohort_per_model":True,
            "new_Holm_families_separate_per_source_and_model":True,
            "hue_as_joint_two_df_not_independent_sin_cos_tests":True,
            "no_spatial_causal_identification":True,
            "original_S1_S2_Holm_20261024_not_replaced":True,
            "H3c_3_and_5_axes_separately_source_verified":True,
        },
        "original_source_gates":verification,
        "S1":run_s1(native),
        "S2":run_s2(focal),
    }
    no_private_keys(out)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,ensure_ascii=False,allow_nan=False,indent=2)+"\n",encoding="utf-8")
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native",type=Path,required=True)
    p.add_argument("--focal",type=Path,required=True)
    p.add_argument("--environment-ext",type=Path,required=True)
    p.add_argument("--measurements",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=execute(a.native,a.focal,a.environment_ext,a.measurements,a.output)
    print(json.dumps({
        "status":r["status"],
        "S1":{k:{"n_tests":v["n_tests"],"n_holm":v["n_holm_005"]}
              for k,v in r["S1"]["families"].items()},
        "S2":{k:{"n_tests":v["n_tests"],"n_holm":v["n_holm_005"]}
              for k,v in r["S2"]["families"].items()},
        "output":"SOURCE_FREE_AGGREGATES_ONLY",
    }))


if __name__=="__main__":
    main()

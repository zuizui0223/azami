"""V4 S1/S2 source-verified 4-var annual/warm vs VIF5/VIF10: all phenotypes.

Compare two registered 4-variable rainfall models, the 5-variable VIFstep5
and 7-variable VIFstep10 source-blind sets. Raw p, BH and Holm in parallel.
Both 4-variable models are *not* outcome-blind VIFstep selections.
The original 12-variable source is held constant. All actual head traits,
including two LIMITED image-derived bract proxies, are fitted as available.
No per-observation identifiers, photographers or coordinates are released.

Observational associations do NOT establish adaptive change, genetic LEGO,
true phyllary spines, mucilage, or gravity-calibrated capitulum direction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

from analysis.v4.s1_s2_extreme_lego import (
    ALL_ENV, FEATURES, FEATURE_COLS, MODULES, PAIRS, read_tables,
    complete_lego_cohort, standardize_traits, split_covariances, rv_rows,
    lego_module_summary, residualize, regression_design,
    diagnostic_matrix, all_diagnostics, require,
    COHORT_ORIGINAL_S1, COHORT_ORIGINAL_S2,
)
from analysis.v4.s1_s2_extreme_range_models import (
    _get_dated_balanced, _formula, MIN_PHOTOS, MIN_INDEPENDENT_CELLS,
    MIN_PHOTOGRAPHERS,
)
from analysis.v4.native_integration import association
from analysis.v4.s1_s2_bract_lego import (
    COLUMNS as BRACT_COMMON, MIN_TAXA as BRACT_MIN_TAXA,
    MIN_N_PER_TAXON, valid_source as bract_source,
)

BR = ("bract_projection_p95", "bract_spread_fraction")
MAIN = ("chelsa_bio01", "chelsa_bio05", "chelsa_bio14",
        "chelsa_bio18", "topo_elevation")
SI = (*MAIN, "chelsa_rsds_mean", "chelsa_bio13")
FOUR_ANNUAL = ("chelsa_bio12", "chelsa_bio05",
               "chelsa_rsds_mean", "topo_elevation")
FOUR_WARM = ("chelsa_bio18", "chelsa_bio05",
             "chelsa_rsds_mean", "topo_elevation")
SPEC = OrderedDict((("four_bio12", FOUR_ANNUAL),
                    ("four_bio18", FOUR_WARM),
                    ("main_vif5", MAIN), ("si_vif10", SI)))
LIMITS = {"four_bio12": 5., "four_bio18": 5.,
          "main_vif5": 5., "si_vif10": 10.}
VIFSTEP_NAMES = {"main_vif5", "si_vif10"}
ENDPOINTS = list(dict.fromkeys(FEATURE_COLS + list(BR)))
SEED = 20261024


def vif_for(diag, names):
    if diag.get("status") != "OK" or not set(names) <= set(diag["variables"]):
        raise ValueError("VIF_SOURCE_NOT_IDENTIFIABLE")
    idx = [diag["variables"].index(name) for name in names]
    corr = np.asarray(diag["pearson"], float)[np.ix_(idx, idx)]
    if np.linalg.eigvalsh(corr).min() < 1.e-8:
        raise ValueError("VIF_SOURCE_SINGULAR")
    values = np.diag(np.linalg.inv(corr))
    return {key: float(value) for key, value in zip(names, values)}


def vifstep_all_contexts(diagnostics, cutoff):
    """Drop the largest VIF across seven source-only contexts, never by trait p."""
    selected = list(ALL_ENV)
    drops = []
    while len(selected) > 1:
        values = [(float(vif), key, context)
                  for context, diagnostic in diagnostics.items()
                  for key, vif in vif_for(diagnostic, selected).items()]
        value, column, context = max(values, key=lambda x: x[0])
        if value <= cutoff:
            return {"status": "VIFSTEP_CONVERGED", "cutoff": cutoff,
                    "retained": selected, "drops": drops, "max_vif": value,
                    "max_vif_context": context, "n_contexts": len(diagnostics)}
        drops.append({"removed": column, "vif_before_removal": value,
                      "source_context": context})
        selected.remove(column)
    raise ValueError("VIFSTEP_NO_MULTIVARIATE_SET")


def source_vif(native, focal):
    env_mask = np.isfinite(native[ALL_ENV].to_numpy(float)).all(axis=1)
    require(int(env_mask.sum()) > 300000, "UNEXPECTED_S1_ENV_COVERAGE")
    d = all_diagnostics(native, ALL_ENV)
    combined = {"S1/" + name: v for name, v in d.items()}
    populations = {}
    for (sp, state), g in focal.groupby(["taxon_name", "range_state"]):
        good = g.loc[np.isfinite(g[ALL_ENV].to_numpy(float)).all(axis=1)]
        require(len(good) >= 500, "S2_SOURCE_GEOGRAPHIC_COVERAGE_CHANGED")
        label = "S2/" + sp + "/" + state
        combined[label] = diagnostic_matrix(good[ALL_ENV], ALL_ENV)
        populations[label] = len(good)
    require(len(combined) == 7 and len(populations) == 4, "SOURCE_CONTEXTS_NOT_SEVEN")
    main = vifstep_all_contexts(combined, 5.)
    si = vifstep_all_contexts(combined, 10.)
    require(set(main["retained"]) == set(MAIN) and len(main["retained"]) == 5,
            "MAIN_5_SOURCE_VIFSTEP_SET_CHANGED")
    require(set(si["retained"]) == set(SI) and len(si["retained"]) == 7,
            "SI_7_SOURCE_VIFSTEP_SET_CHANGED")
    receipt = {"source_only_before_any_trait_screen": True,
               "original_S1_all12_environment_rows": int(env_mask.sum()),
               "original_S2_stratum_complete_rows": populations,
               "7_source_contexts": {name: int(diag["n"])
                                     for name, diag in combined.items()},
               "selection": {"main_vif5": main, "si_vif10": si},
               "fixed_model_set_diagnostics": {}}
    for name, variables in SPEC.items():
        vv = {key: vif_for(v, variables) for key, v in combined.items()}
        m = max((value, key, col) for key, row in vv.items()
                for col, value in row.items())
        if name in VIFSTEP_NAMES:
            require(m[0] <= LIMITS[name] + 1.e-6,
                    "VIFSTEP_SELECTED_SET_EXCEEDS_SOURCE_THRESHOLD")
        # The BIO12 four-variable comparison is deliberately retained even
        # when C. arvense/introduced source VIF is slightly above 5.
        # This is NOT a VIFstep-selected model; report that failure explicitly.
        receipt["fixed_model_set_diagnostics"][name] = {
            "columns": list(variables), "diagnostic_vif_cutoff": LIMITS[name],
            "selected_by_vifstep": name in VIFSTEP_NAMES,
            "source_vif_pass": bool(m[0] <= LIMITS[name]+1.e-6),
            "max_vif": m[0], "max_vif_context": m[1],
            "max_vif_variable": m[2],
            "per_context_vif": vv}
    return receipt


def selected_source(frame, traits, *, min_taxa=20):
    """Use QC flags in S1; do not impute absent traits or locations."""
    if not all(c in frame and c + "__eligible" in frame for c in traits):
        return None, {"status": "ORIGINAL_S1_TRAIT_OR_QC_UNAVAILABLE",
                      "required": list(traits)}
    mask = np.isfinite(frame[ALL_ENV].to_numpy(float)).all(axis=1)
    for col in traits:
        mask &= np.isfinite(pd.to_numeric(frame[col], errors="coerce").to_numpy(float))
        mask &= frame[col + "__eligible"].fillna(False).astype(bool).to_numpy()
    d = frame.loc[mask].copy()
    if "orientation_image_vertical_angle" in traits:
        d = d.loc[d.orientation_image_vertical_angle.between(0, 180)].copy()
    counts = d.groupby("accepted_key").size()
    d = d.loc[d.accepted_key.isin(counts[counts >= 5].index)].copy()
    coverage = {"n": len(d), "taxa": int(d.accepted_key.nunique()),
                "matched_original_12_climate": True, "min_taxon_n": 5}
    if coverage["taxa"] < min_taxa:
        return None, {"status": "INSUFFICIENT_S1_TRAIT_TAXA", **coverage}
    return d, {"status": "SOURCE_LINKED", **coverage}


def standardized_y(y):
    arr = np.asarray(y, float)
    sd = arr.std()
    if not np.isfinite(sd) or sd < 1.e-10:
        raise ValueError("SINGLE_TRAIT_ZERO_VARIANCE")
    return (arr-arr.mean())/sd


def _cluster_cross_cov(fit_a, fit_b, design, weight, groups):
    """Cross-sandwich from same-design clustered WLS for joint circular hue."""
    if len(groups) != len(design):
        raise ValueError("INVALID_CLUSTER_DESIGN")
    codes, distinct = pd.factorize(groups, sort=True)
    n = len(groups)
    G = len(distinct)
    p = design.shape[1]
    if G <= p + 2:
        raise ValueError("TOO_FEW_CLUSTERS_FOR_JOINT_HUE")
    gram_inv = np.linalg.inv(design.T @ (weight[:, None] * design))
    residual_a = np.asarray(fit_a.resid)
    residual_b = np.asarray(fit_b.resid)
    score_a = np.zeros((G, p))
    score_b = np.zeros((G, p))
    np.add.at(score_a, codes, (weight*residual_a)[:, None]*design)
    np.add.at(score_b, codes, (weight*residual_b)[:, None]*design)
    correction = (G/(G-1))*((n-1)/(n-p))
    return correction*gram_inv @ (score_a.T @ score_b) @ gram_inv


def hue_joint(fits, design, weights, groups, variables, scales):
    if len(fits) != 2:
        raise ValueError("HUE_MUST_HAVE_SIN_AND_COS")
    left, right = fits
    cross = _cluster_cross_cov(left, right, design, weights, groups)
    out = []
    for suffix, offset in scales:
        for i, v in enumerate(variables):
            col = offset + i
            covariance = np.array([[left.cov_params()[col, col], cross[col, col]],
                                   [cross[col, col], right.cov_params()[col, col]]])
            betas = np.array([left.params[col], right.params[col]])
            if not np.isfinite(covariance).all() or np.linalg.eigvalsh(covariance).min() <= 1.e-14:
                out.append({"environment": v, "scale": suffix,
                            "status": "JOINT_HUE_COVARIANCE_UNIDENTIFIABLE"})
                continue
            stat = float(betas @ np.linalg.solve(covariance, betas))
            out.append({"environment": v, "scale": suffix,
                        "status": "JOINT_CIRCULAR_HUE_WALD_2DF",
                        "wald_chi2": max(0., stat),
                        "p_uncorrected": float(chi2.sf(max(0., stat), 2)),
                        "n_coordinates": 2,
                        "not_two_independent_hue_traits": True})
    return out


def cohort_vif(frame, variables, spec):
    """Re-diagnose the matched trait subset without reselecting variables."""
    diagnostic = diagnostic_matrix(frame[list(variables)], list(variables))
    vif = max(diagnostic.get("vif", {}).values()) if diagnostic.get("status") == "OK" else None
    return {"max_vif_on_trait_cohort": vif,
            "source_blind_cutoff": LIMITS[spec],
            "trait_cohort_vif_pass": bool(vif is not None and vif <= LIMITS[spec] + 1.e-5),
            "not_used_to_select_or_drop_predictors": True}


def _coefficient_rows(fit, variables, scales):
    out = []
    for label, offset in scales:
        for i, variable in enumerate(variables):
            j = offset+i
            coef, se, p = float(fit.params[j]), float(fit.bse[j]), float(fit.pvalues[j])
            out.append({"environment": variable, "scale": label,
                        "beta_std": coef, "se_cluster": se,
                        "p_uncorrected": p,
                        "ci95_wald": [coef-1.95996398454*se,
                                      coef+1.95996398454*se],
                        "status": "ESTIMATED" if np.isfinite(p) else "NONFINITE_P"})
    return out


def estimate_s1_blocks(native):
    output = {}
    blocks = OrderedDict(FEATURES)
    blocks["bract_projection"] = [BR[0]]
    blocks["bract_spread"] = [BR[1]]
    for block, traits in blocks.items():
        part, coverage = selected_source(native, traits, min_taxa=20)
        if part is None:
            output[block] = coverage
            continue
        answer = {**coverage, "trait_coordinates": traits, "models": {}}
        groups = part.accepted_key.astype(str).to_numpy()
        n_taxa = part.accepted_key.nunique()
        counts = part.groupby("accepted_key").size()
        weights = part.accepted_key.map(1/counts).to_numpy(float)
        Y = [standardized_y(part[c].to_numpy(float)) for c in traits]
        for spec, variables in SPEC.items():
            status = cohort_vif(part, variables, spec)
            try:
                design, _ = regression_design(part, list(variables))
                fits = [sm.WLS(y, design, weights=weights).fit(
                    cov_type="cluster", cov_kwds={"groups": groups}, use_t=False)
                    for y in Y]
                components = {
                    col: _coefficient_rows(fit, variables,
                        [("between_equal_taxa", 1),
                         ("within_taxon", 1+len(variables))])
                    for col, fit in zip(traits, fits)}
                circular = hue_joint(fits, design, weights, groups, variables,
                    [("between_equal_taxa", 1),
                     ("within_taxon", 1+len(variables))]) if block == "hue" else []
                answer["models"][spec] = {**status, "status": "ESTIMATED",
                    "predictors": list(variables), "n": len(part),
                    "n_taxa": int(n_taxa), "component_effects": components,
                    "joint_circular_hue_tests": circular,
                    "equal_taxon_total_weight": True,
                    "accepted_key_cluster_robust": True}
            except (ValueError, np.linalg.LinAlgError) as exc:
                answer["models"][spec] = {**status,
                    "status": "NOT_IDENTIFIABLE", "reason": str(exc)[:180]}
        output[block] = answer
    return output


def s1_joint_28(native):
    # Restrict only multi-trait relationships to the SAME eligible photographs.
    before = complete_lego_cohort(native)
    cohort = before.loc[np.isfinite(before[ALL_ENV].to_numpy(float)).all(axis=1)].copy()
    cnt = cohort.groupby("accepted_key").size()
    cohort = cohort.loc[cohort.accepted_key.isin(cnt[cnt>=5].index)].copy()
    if len(cohort)<1000 or cohort.accepted_key.nunique()<20:
        raise ValueError("V4_S1_MULTI_TRAIT_COHORT_MISSING")
    v = standardize_traits(cohort)
    within, among, *_ = split_covariances(v, cohort.accepted_key)
    raw = rv_rows(within, among)
    output = {"n": len(cohort), "taxa": int(cohort.accepted_key.nunique()),
              "n_visible_constructs": 8, "n_pairs": len(raw),
              "raw_pairs": raw, "models": {}}
    for spec, variables in SPEC.items():
        try:
            res, *_ = residualize(cohort, v, list(variables))
            w, a, *_ = split_covariances(res, cohort.accepted_key)
            pairs = rv_rows(w, a)
            output["models"][spec] = {**cohort_vif(cohort, variables, spec),
                "status": "ESTIMATED", "predictors": list(variables),
                "n_pairs": len(pairs),
                "n_among_stronger_raw": int(sum(x["among_minus_within"]>0 for x in raw)),
                "n_among_stronger_adjusted": int(sum(x["among_minus_within"]>0 for x in pairs)),
                "module_association": lego_module_summary(pairs),
                "pairs": pairs, "not_evolutionary_transitions": True}
        except (ValueError, np.linalg.LinAlgError) as exc:
            output["models"][spec] = {"status":"NOT_IDENTIFIABLE",
                                      "reason": str(exc)[:180]}
    return output


def s1_bract_edges(native):
    # Bracts are limited original V4 PHOTO PROJECTION proxies, no defence trait.
    links = {}
    from analysis.v4.s1_s2_extreme_lego import FEATURES as CORE_FEATURES
    for core, cols in CORE_FEATURES.items():
        traits = list(dict.fromkeys(cols + list(BR)))
        cohort, cover = selected_source(native, traits, min_taxa=BRACT_MIN_TAXA)
        if cohort is None:
            links[core] = cover
            continue
        Y = cohort[traits].to_numpy(float)
        std = Y.std(axis=0)
        if min(std) <= 1.e-9:
            links[core] = {"status": "INVALID_TRAIT_SCALE", **cover}
            continue
        Y = (Y-Y.mean(axis=0))/std
        ncore = len(cols)
        left = list(range(ncore))
        right = list(range(ncore,ncore+len(BR)))
        def compare(dat):
            within, among, *_ = split_covariances(dat, cohort.accepted_key)
            return {"within_rv": float(association(within,left,right)[0]),
                    "among_rv": float(association(among,left,right)[0])}
        result = {**cover, "status": "OBSERVED_PHOTO_BRACT_BRIDGE",
                  "core_construct": core,
                  "bract_proxy_columns": list(BR),
                  "raw": compare(Y), "models": {}}
        for spec, variables in SPEC.items():
            try:
                adj,*_ = residualize(cohort, Y, list(variables))
                result["models"][spec] = {"status": "ESTIMATED",
                     "predictors": list(variables),
                     **cohort_vif(cohort, variables, spec),
                     **compare(adj)}
            except (ValueError, np.linalg.LinAlgError) as exc:
                result["models"][spec] = {"status": "NOT_IDENTIFIABLE",
                                          "reason": str(exc)[:180]}
        links[core] = result
    return {"status": "PHOTO_PROXY_ONLY_NO_SPINES_NO_MUCILAGE",
            "bridges_to_all_eight_original_constructs": links,
            "pair_specific_photo_cohorts_NOT_same_as_28core": True}


def s2_extra_bract(focal, original_measurement_path):
    """Only include original bracts when original S2 archived source has QC."""
    if original_measurement_path is None:
        return focal, {"status": "S2_BRACT_SOURCE_NOT_SUPPLIED"}
    m = pd.read_csv(original_measurement_path,
                    dtype={"obs_id":str}, low_memory=False)
    if not m.obs_id.is_unique:
        raise ValueError("S2_BRACT_MEASUREMENT_SOURCE_DUPLICATE")
    matched = {}
    for endpoint in BR:
        aliases = (endpoint, endpoint+"_median")
        eligible = [(s, s+"__eligible") for s in aliases
                    if s in m and s+"__eligible" in m]
        if not eligible:
            # Median without documented eligible indicator cannot be treated
            # as quality controlled. No false negative claim of biological absence.
            return focal, {"status": "S2_BRACT_NOT_VERIFIABLY_QC_MEASURED",
                           "required": list(BR),
                           "matched_columns": list(matched)}
        matched[endpoint] = eligible[0]
    columns = ["obs_id"] + [entry for a, b in matched.values() for entry in (a, b)]
    rename_map = {}
    for canonical, (measurement, flag) in matched.items():
        rename_map[measurement] = canonical
        rename_map[flag] = canonical + "__eligible"
    mm = m[columns].rename(columns=rename_map)
    # Source WCVP eligibility and dates remain the focal source, never inferred.
    out = focal.merge(mm, on="obs_id", how="left", validate="one_to_one")
    if len(out)!=len(focal):
        raise ValueError("S2_BRACT_JOIN_CHANGED_GRAIN")
    return out, {"status":"S2_BRACT_ORIGINAL_QC_AVAILABLE",
                 "source_columns": matched}


def s2_trait_blocks(focal, bract_receipt):
    blocks = OrderedDict(FEATURES)
    if bract_receipt["status"] == "S2_BRACT_ORIGINAL_QC_AVAILABLE":
        blocks["bract_projection"] = [BR[0]]
        blocks["bract_spread"] = [BR[1]]
    table = {}
    for (species,state), g in focal.groupby(["taxon_name","range_state"]):
        all12 = np.isfinite(g[ALL_ENV].to_numpy(float)).all(axis=1)
        basis = g.loc[all12].copy()
        answer = {"original_wild_observation_n": len(g),
                  "all12_environment_n": len(basis),
                  "source_measured_visible_traits": {}, "joint_28_pairs": {}}
        for block, traits in blocks.items():
            if not set(traits) <= set(basis):
                answer["source_measured_visible_traits"][block] = {
                    "status": "NO_SOURCE_LINKED_MEASUREMENT_COLUMNS"}
                continue
            part = basis.loc[np.isfinite(basis[traits].to_numpy(float)).all(axis=1)].copy()
            if "orientation_image_vertical_angle" in traits:
                part = part.loc[part.orientation_image_vertical_angle.between(0,180)].copy()
            if block.startswith("bract"):
                part = part.loc[part[traits[0]+"__eligible"].fillna(False).astype(bool)].copy()
            raw_n = len(part)
            cohort = _get_dated_balanced(part)
            coverage = {"source_trait_n": raw_n, "balanced_n": len(cohort),
                        "balanced_cells": int(cohort._cell.nunique()) if len(cohort) else 0,
                        "balanced_observers": int(cohort.user_id.nunique()) if len(cohort) else 0,
                        "one_observer_per_5deg_cell": True}
            if (len(cohort)<MIN_PHOTOS or cohort._cell.nunique()<MIN_INDEPENDENT_CELLS
                    or cohort.user_id.nunique()<MIN_PHOTOGRAPHERS):
                answer["source_measured_visible_traits"][block] = {
                    **coverage, "status":"INSUFFICIENT_BALANCED_S2_SUPPORT"}
                continue
            results = {**coverage,"status":"SOURCE_LINKED", "models": {}}
            groups = cohort._cell.to_numpy()
            weights = cohort._cell.map(1/cohort.groupby("_cell").size()).to_numpy(float)
            for spec, vars_tuple in SPEC.items():
                vif = cohort_vif(cohort, vars_tuple, spec)
                try:
                    design, mean, sd, retained, omitted = _formula(cohort,list(vars_tuple))
                    fits = [sm.WLS(standardized_y(cohort[col]),design,
                        weights=weights).fit(cov_type="cluster",
                        cov_kwds={"groups":groups},use_t=False) for col in traits]
                    tests = {col: _coefficient_rows(fit, vars_tuple,
                        [("within_species_range", 1)])
                        for col,fit in zip(traits,fits)}
                    circular = (hue_joint(fits,design,weights,groups,vars_tuple,
                        [("within_species_range",1)]) if block=="hue" else [])
                    results["models"][spec] = {**vif, "status":"ESTIMATED",
                        "predictors": list(vars_tuple),
                        "component_effects":tests,
                        "joint_circular_hue_tests":circular,
                        "nuisance_month_year_retained":retained,
                        "nuisance_month_year_omitted_rank_deficient":omitted,
                        "equal_5deg_cell_weights_cluster_robust":True}
                except (ValueError,np.linalg.LinAlgError) as exc:
                    results["models"][spec] = {**vif,
                        "status": "NOT_IDENTIFIABLE", "reason":str(exc)[:180]}
            answer["source_measured_visible_traits"][block] = results

        # S2 covariance of all eight observable constructs on exactly same
        # balanced records, distinct from S1 across-taxa/within-taxa contrast.
        core = basis.loc[np.isfinite(basis[FEATURE_COLS].to_numpy(float)).all(axis=1)].copy()
        core = core.loc[core.orientation_image_vertical_angle.between(0,180)].copy()
        core = _get_dated_balanced(core)
        if (len(core)<MIN_PHOTOS or core._cell.nunique()<MIN_INDEPENDENT_CELLS
                or core.user_id.nunique()<MIN_PHOTOGRAPHERS):
            answer["joint_28_pairs"] = {"status":"INSUFFICIENT_S2_JOINT_SUPPORT",
                                        "n":len(core)}
        else:
            response = core[FEATURE_COLS].to_numpy(float)
            std = response.std(axis=0)
            std_hue = float(np.sqrt((std[FEATURE_COLS.index("corolla_hue_sin")]**2+
                                      std[FEATURE_COLS.index("corolla_hue_cos")]**2)/2))
            std[FEATURE_COLS.index("corolla_hue_sin")] = std_hue
            std[FEATURE_COLS.index("corolla_hue_cos")] = std_hue
            if np.min(std)<1.e-9:
                raise ValueError("S2_JOINT_TRAIT_COORDINATE_SCALE")
            y=(response-response.mean(axis=0))/std
            ng=core.groupby("_cell").size()
            wt=core._cell.map(1/ng).to_numpy(float)
            def pairs_of(arr):
                mean=np.average(arr,axis=0,weights=wt)
                centered=arr-mean
                cov=(centered*wt[:,None]).T@centered/wt.sum()
                return rv_rows(cov,cov)
            raw_pairs=pairs_of(y)
            joint = {"status":"S2_WITHIN_RANGE_PHENOTYPE_COVARIATION",
                     "n":len(core),"cells":len(ng),"n_pairs":28,
                     "raw_pairs":raw_pairs, "models":{}}
            for spec,vars_tuple in SPEC.items():
                try:
                    X, *_ = _formula(core,list(vars_tuple))
                    beta = np.linalg.lstsq((X*np.sqrt(wt[:,None])),
                                           y*np.sqrt(wt[:,None]),rcond=None)[0]
                    residual = y-X@beta
                    pairs = pairs_of(residual)
                    joint["models"][spec] = {**cohort_vif(core,vars_tuple,spec),
                        "status":"ESTIMATED", "predictors":list(vars_tuple),
                        "pairs":pairs, "not_between_species_modularity":True}
                except (ValueError,np.linalg.LinAlgError) as exc:
                    joint["models"][spec]={"status":"NOT_IDENTIFIABLE",
                                           "reason":str(exc)[:180]}
            answer["joint_28_pairs"]=joint
        table[species+"/"+state]=answer
    return table


def family_corrections(s1, s2):
    """Four transparent families, excluding duplicate hue components as claims."""
    families = {}
    for layer, source in (("S1",s1),("S2",s2)):
        for spec in SPEC:
            rows=[]
            if layer=="S1":
                members=[(key, v) for key,v in source.items()]
            else:
                members=[(key+"/"+block, outcome) for key,v in source.items()
                         for block,outcome in v["source_measured_visible_traits"].items()]
            for label, entry in members:
                model=entry.get("models",{}).get(spec,{})
                if model.get("status")!="ESTIMATED":
                    continue
                for col, effects in model.get("component_effects",{}).items():
                    # Biological colour hue is one 2-dimensional trait.
                    # Retain individual sin/cos estimates in component_effects,
                    # but never double-count them alongside the joint 2df test.
                    if col in ("corolla_hue_sin","corolla_hue_cos"):
                        continue
                    for v in effects:
                        if v["status"]=="ESTIMATED" and np.isfinite(v["p_uncorrected"]):
                            rows.append({"group":label,"coordinate":col,
                                         "environment":v["environment"],"scale":v["scale"],
                                         "p_original":v["p_uncorrected"],
                                         "beta_std":v["beta_std"],
                                         "ci95_wald":v["ci95_wald"],
                                         "n":int(model.get("n",entry.get("balanced_n",0))),
                                         "trait_cohort_vif_pass":bool(
                                             model["trait_cohort_vif_pass"]),
                                         "test_kind":"coordinate_exploratory"})
                for v in model.get("joint_circular_hue_tests",[]):
                    if v["status"] == "JOINT_CIRCULAR_HUE_WALD_2DF":
                        rows.append({"group":label,"coordinate":"hue_joint_sin_cos",
                                     "environment":v["environment"],"scale":v["scale"],
                                     "p_original":v["p_uncorrected"],
                                     "beta_std":None,
                                     "ci95_wald":None,
                                     "n":int(model.get("n",entry.get("balanced_n",0))),
                                     "trait_cohort_vif_pass":bool(
                                         model["trait_cohort_vif_pass"]),
                                     "test_kind":"joint_circular_hue_2df"})
            if rows:
                adj_bh = multipletests([r["p_original"] for r in rows],
                                      method="fdr_bh")[1]
                adj_holm = multipletests([r["p_original"] for r in rows],
                                        method="holm")[1]
                for row,pbh,pholm in zip(rows,adj_bh,adj_holm):
                    row["p_BH"] = float(pbh)
                    row["p_Holm"] = float(pholm)
                    row["raw_lt_0_05"] = bool(row["p_original"]<.05)
                    row["BH_lt_0_05"] = bool(pbh<.05)
                    row["Holm_lt_0_05"] = bool(pholm<.05)
            families[layer+"/"+spec] = {
                "status":"COMPLETE_NEW_EXPLORATORY_HYPOTHESIS_FAMILY",
                "n_tests":len(rows),
                "n_unadjusted":sum(x["raw_lt_0_05"] for x in rows),
                "n_BH":sum(x["BH_lt_0_05"] for x in rows),
                "n_Holm":sum(x["Holm_lt_0_05"] for x in rows),
                "n_raw_with_model_VIF_pass":sum(
                    x["raw_lt_0_05"] and x["trait_cohort_vif_pass"] for x in rows),
                "n_BH_with_model_VIF_pass":sum(
                    x["BH_lt_0_05"] and x["trait_cohort_vif_pass"] for x in rows),
                "n_Holm_with_model_VIF_pass":sum(
                    x["Holm_lt_0_05"] and x["trait_cohort_vif_pass"] for x in rows),
                "rows":rows,
                "hue_interpretation": "Only joint sin/cos 2df hue tests count; individual sin/cos estimates preserved solely as diagnostics.",
                "family_scope": "one V4 layer and one model: all eligible actual biological constructs, all environmental terms, both S1 scales if applicable",
                "q_values_per_model_not_joint_correction_for_selecting_the_best_of_4_models": True}
    return families


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--native",type=Path,required=True)
    ap.add_argument("--focal",type=Path,required=True)
    ap.add_argument("--environment-ext",type=Path,required=True)
    ap.add_argument("--s2-measurements",type=Path)
    ap.add_argument("--out-public",type=Path,required=True)
    args=ap.parse_args()
    args.out_public.mkdir(parents=True, exist_ok=True)
    native,focal = read_tables(args.native,args.focal,args.environment_ext,
                               args.s2_measurements)
    if not set(ALL_ENV) <= set(focal):
        from analysis.v4.s1_s2_restore_original_five import recover_original_five
        focal,_=recover_original_five(focal,args.out_public.parent/"private"/"original_s2_five",
                                      Path("analysis/v4/environment_source_metadata.json"))
    require(len(native)==COHORT_ORIGINAL_S1 and len(focal)==COHORT_ORIGINAL_S2,
            "ORIGINAL_S1_S2_SOURCE_IDENTITY_CHANGED")
    selection=source_vif(native,focal)
    focal,bract_receipt=s2_extra_bract(focal,args.s2_measurements)
    s1=estimate_s1_blocks(native)
    joint=s1_joint_28(native)
    bract=s1_bract_edges(native)
    s2=s2_trait_blocks(focal,bract_receipt)
    families=family_corrections(s1,s2)
    output={"analysis_id":"azami_v4_s1_s2_four_bio12_bio18_vs_vifstep5_vifstep10_alltraits_20261024",
            "status":"SOURCE_VERIFIED_S1_S2_FOUR_MODELS_ALL_AVAILABLE_PHENOTYPES",
            "V4_NOT_V3":True,"selection":selection,
            "S1": {"original_native_n":len(native),"traits":s1,
                   "joint_8_constructs_28_pairs":joint,
                   "bract_to_all_8_constructs":bract},
            "S2": {"original_wild_source_n":len(focal),
                   "optional_bract_source":bract_receipt,
                   "strata":s2},
            "multiplicity":families,
            "analysis_plan":{
                "four_bio12": "BIO12 + BIO05 + RSDS + elevation; registered non-VIFstep sensitivity",
                "four_bio18": "BIO18 + BIO05 + RSDS + elevation; registered non-VIFstep sensitivity",
                "main_vif5": "source-only stepwise max-VIF cutoff 5 from 12 candidate predictors",
                "si_vif10": "source-only stepwise max-VIF cutoff 10 from same candidates",
                "cohorts": "all12-complete source with outcome-specific QC, identical per biological trait across four models",
                "inferential_contrasts": "unadjusted p, BH q and Holm p side by side for distinct biological endpoint tests; 95% pointwise Wald CI",
                "corrected_results_not_precondition_for_reporting": True,
            },
            "interpretation_boundaries":{
                "all_P_exploratory_and_spatial_phylogeny_unresolved":True,
                "image_orientation_not_gravity_calibrated":True,
                "bract_projection_and_spread_not_true_spines_or_mucilage":True,
                "same_photo_outline_colour_calculated_coordinates_not_independent":True,
                "S1_S2_native_observations_overlap":True,
                "joint_circular_hue_reported_separately":True,
                "no_heritable_genomic_LEGO_or_fitness_test":True},
            "no_public_observation_identifiers_coordinates_or_photographers":True}
    path=args.out_public/"V4_S1_S2_four_vs_vif5_vif10_alltraits.json"
    path.write_text(json.dumps(output,ensure_ascii=False,allow_nan=False,indent=2)+"\n",
                    encoding="utf-8")
    print(json.dumps({"status":output["status"],
       "four_bio12":list(FOUR_ANNUAL),"four_bio18":list(FOUR_WARM),
       "vif5":list(MAIN),"vif10":list(SI),
       "source_S1":len(native),"source_S2":len(focal),
       "S1_trait_blocks":{k:v["status"] for k,v in s1.items()},
       "S1_common":joint["n"],"S2_strata":list(s2),
       "hypothesis_family":{k:{"n":v["n_tests"],"raw":v["n_unadjusted"],
                                    "BH":v["n_BH"],"holm":v["n_Holm"]}
                            for k,v in families.items()}},
       ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

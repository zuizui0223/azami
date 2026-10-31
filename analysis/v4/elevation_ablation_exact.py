"""Matched-cohort elevation-ablation sensitivity for V4 S1/S2.

M3: CHELSA BIO18 + BIO05 + RSDS (climate only).
M4: identical predictors plus WorldClim elevation (original frozen model).
Never choose M3/M4 by significance, never overwrite the 2026-10-24 Holm
families, never export observation, taxon-key, camera-user or location records.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from analysis.v4.s1_s2_extreme_lego import (
    ALL_ENV, read_tables, regression_design, diagnostic_matrix, require
)
from analysis.v4.s1_s2_vif5_si10_alltraits import (
    FOUR_WARM, selected_source, standardized_y
)
from analysis.v4.s2_original_colour_spatial_block_sensitivity import (
    EXPECTED, original_sample, _original_fit
)

M3 = tuple(FOUR_WARM[:-1])
M4 = tuple(FOUR_WARM)
ANGLE = "orientation_image_vertical_angle"
EXPECTED_S1_N = 168907
EXPECTED_S1_TAXA = 178
EXPECTED_S1_BIO18_BETA = 0.04544684307686103
FORBIDDEN_OUTPUT_KEYS = {
    "obs_id", "photo_id", "user_id", "photographer_id",
    "latitude", "longitude", "accepted_key", "source_obs_id"
}


def wls_design_response(x: np.ndarray, y: np.ndarray, weight: np.ndarray, group: np.ndarray):
    model = sm.WLS(y, x, weights=weight).fit(
        cov_type="cluster", cov_kwds={"groups": group}, use_t=False
    )
    if not np.isfinite(model.params).all() or not np.isfinite(model.bse).all():
        raise ValueError("NONFINITE_SOURCE_MODEL")
    return model


def summary(fit, predictors, *, offset=1):
    return {
        name: {
            "beta_std": float(fit.params[i + offset]),
            "cluster_se": float(fit.bse[i + offset]),
            "raw_p": float(fit.pvalues[i + offset]),
            "pointwise_95ci": [float(x) for x in fit.conf_int()[i + offset]],
        }
        for i, name in enumerate(predictors)
    }


def diagnostics(frame, predictors):
    stat = diagnostic_matrix(frame[list(predictors)], list(predictors))
    if stat["status"] != "OK":
        raise ValueError("CLIMATE_COLLINEARITY_SINGULAR")
    return {
        "max_vif": float(max(stat["vif"].values())),
        "vif": {k: float(v) for k, v in stat["vif"].items()},
        "max_abs_pair_r": float(stat["max_abs_pair_r"]),
        "note": "Diagnostics on the same eligibility-restricted observations",
    }


def compare_s1(native):
    # The 12-covariate missingness gate remains IDENTICAL in both models.
    sample, coverage = selected_source(native, [ANGLE], min_taxa=20)
    if (sample is None or len(sample) != EXPECTED_S1_N or
            sample.accepted_key.nunique() != EXPECTED_S1_TAXA):
        raise ValueError("S1_FROZEN_COHORT_NOT_REPRODUCED")
    source_groups = sample.accepted_key.astype(str).to_numpy()
    count = sample.groupby("accepted_key").size()
    weight = sample.accepted_key.map(1 / count).to_numpy(float)
    y = standardized_y(sample[ANGLE].to_numpy(float))
    fits = {}
    for key, variables in (("M3_climate_only", M3), ("M4_original_elevation", M4)):
        x, group = regression_design(sample, list(variables))
        if not np.array_equal(group, source_groups):
            raise ValueError("S1_GROUPS_CHANGED_BETWEEN_MODELS")
        fit = wls_design_response(x, y, weight, group)
        fits[key] = fit
    original = float(fits["M4_original_elevation"].params[1])
    if not np.isclose(original, EXPECTED_S1_BIO18_BETA, atol=3e-7, rtol=0):
        raise ValueError("S1_FROZEN_BIO18_COEFFICIENT_NOT_REPRODUCED")

    result = {}
    for key, variables in (("M3_climate_only", M3), ("M4_original_elevation", M4)):
        fit = fits[key]
        result[key] = {
            "n": len(sample), "n_taxa": EXPECTED_S1_TAXA,
            "predictors": list(variables),
            "between_equal_taxa": summary(fit, variables, offset=1),
            "within_taxon": summary(fit, variables, offset=1 + len(variables)),
            "aic_same_sample": float(fit.aic),
            "weighted_sse": float(np.sum(weight * fit.resid**2)),
            "diagnostics": diagnostics(sample, variables),
        }
    m3, m4 = fits["M3_climate_only"], fits["M4_original_elevation"]
    return {
        "status": "EXACT_FROZEN_S1_REPRODUCED_NEW_EXPLORATORY_M3",
        "shared_cohort": {"n": len(sample), "taxa": EXPECTED_S1_TAXA},
        "models": result,
        "changes": {
            "BIO18_between_beta_M3_minus_M4": float(m3.params[1] - m4.params[1]),
            "BIO18_between_sign_same": bool(np.sign(m3.params[1]) == np.sign(m4.params[1])),
            "AIC_M3_minus_M4": float(m3.aic - m4.aic),
            "partial_weighted_sse_reduction_from_elevation": float(
                1 - (np.sum(weight * m4.resid**2) / np.sum(weight * m3.resid**2))
            ),
        },
        "frozen_S1_80_Holm_family_unchanged": True,
        "new_three_predictor_p_values_exploratory_uncorrected": True,
    }


def compare_s2(focal):
    from analysis.v4.s1_s2_extreme_range_models import _formula
    rows = []
    for species, state, phenotype, original_exposure, original_n, original_beta in EXPECTED:
        sample, cohort = original_sample(focal, species, state, phenotype)
        if len(sample) != original_n:
            raise ValueError("S2_FROZEN_PHOTO_COUNT_NOT_REPRODUCED")
        X4, y4, wt4, index4, original = _original_fit(sample, phenotype, original_exposure)
        if not np.isclose(original["beta"], original_beta, atol=3e-7, rtol=0):
            raise ValueError("S2_ORIGINAL_FROZEN_COEFFICIENT_NOT_REPRODUCED")
        X3, _, _, nuis3, omitted3 = _formula(sample, list(M3))
        X4bis, _, _, nuis4, omitted4 = _formula(sample, list(M4))
        if not np.allclose(X4, X4bis, rtol=0, atol=1e-12):
            raise ValueError("S2_BASELINE_DESIGN_CHANGED")
        if nuis3 != nuis4 or omitted3 != omitted4:
            raise ValueError("S2_NUISANCE_DESIGN_DIFFERS_BETWEEN_NESTED_MODELS")
        groups = sample._cell.to_numpy()
        model3 = wls_design_response(X3, y4, wt4, groups)
        model4 = wls_design_response(X4, y4, wt4, groups)
        m3_sse = float(np.sum(wt4 * model3.resid**2))
        m4_sse = float(np.sum(wt4 * model4.resid**2))
        if m4_sse > m3_sse + 1e-7:
            raise ValueError("S2_NESTED_ELEVATION_MODEL_NOT_NESTED")
        rows.append({
            "species": species,
            "range": state,
            "trait": phenotype,
            "n": original_n,
            "n_geographic_cells": int(sample._cell.nunique()),
            "frozen_original_exposure": original_exposure,
            "original_exposure_beta_reproduced": float(original["beta"]),
            "M3_climate_only": {
                "predictors": list(M3), "coefficients": summary(model3, M3),
                "aic_same_sample": float(model3.aic),
                "diagnostics": diagnostics(sample, M3),
            },
            "M4_original_elevation": {
                "predictors": list(M4), "coefficients": summary(model4, M4),
                "aic_same_sample": float(model4.aic),
                "diagnostics": diagnostics(sample, M4),
            },
            "AIC_M3_minus_M4": float(model3.aic - model4.aic),
            "partial_weighted_sse_reduction_from_elevation": float(1 - m4_sse / m3_sse),
            "same_observations_weights_and_dates_in_both_models": True,
            "no_native_vs_introduced_interaction_test": True,
            "new_three_predictor_p_values_exploratory_uncorrected": True,
        })
    return {
        "status": "ALL_SIX_FROZEN_S2_COEFFICIENTS_REPRODUCED_NEW_EXPLORATORY_M3",
        "strata": rows,
        "frozen_S2_128_Holm_family_unchanged": True,
        "no_claim_about_introduction_caused_evolution": True,
    }


def check_public_payload(result):
    def walk(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                if key in FORBIDDEN_OUTPUT_KEYS:
                    raise ValueError("PRIVATE_SOURCE_KEY_IN_PUBLIC_JSON:" + key)
                walk(value)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)
    walk(result)
    if result["S1"]["status"] != "EXACT_FROZEN_S1_REPRODUCED_NEW_EXPLORATORY_M3":
        raise ValueError("S1_FROZEN_GATE_DID_NOT_PASS")
    if len(result["S2"]["strata"]) != 6:
        raise ValueError("MISSING_ONE_OF_SIX_FROZEN_S2_COMPARISONS")


def run(native_path: Path, focal_path: Path, extra_path: Path,
        measurements_path: Path, output_path: Path):
    native, focal = read_tables(native_path, focal_path, extra_path, measurements_path)
    if not set(ALL_ENV) <= set(native) or not set(ALL_ENV) <= set(focal):
        raise ValueError("MISSING_ALL12_ENVIRONMENT_FOR_SAME_COHORT")
    result = {
        "schema": "V4_EXACT_MATCHED_ELEVATION_ABLATION_20261031",
        "status": "NEW_EXPLORATORY_SENSITIVITY_NOT_ORIGINAL_HOLM_INFERENCE",
        "original_20261024_frozen_estimates_unchanged": True,
        "comparison": {
            "M3": list(M3), "M4": list(M4),
            "cohort_restricted_to_original_all12_climate_complete": True,
            "no_new_tuning_by_trait_p_values": True,
            "no_elevation_exposure_coefficient_in_M3": True,
            "no_spatial_or_phylogenetic_causal_identification": True,
        },
        "S1": compare_s1(native),
        "S2": compare_s2(focal),
    }
    check_public_payload(result)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--focal", required=True, type=Path)
    parser.add_argument("--environment-ext", required=True, type=Path)
    parser.add_argument("--s2-measurements", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    r = run(args.native, args.focal, args.environment_ext,
            args.s2_measurements, args.output)
    print(json.dumps({
        "status": r["status"], "S1_n": r["S1"]["shared_cohort"]["n"],
        "S1_taxa": r["S1"]["shared_cohort"]["taxa"],
        "S2_rows": len(r["S2"]["strata"]),
        "export": "AGGREGATE_ONLY"
    }))


if __name__ == "__main__":
    main()

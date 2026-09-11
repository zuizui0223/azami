#!/usr/bin/env python3
"""Fail-closed validation for the GEB v2 full-27/full-environment atlas."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


EXPECTED_PREDICTORS = [
    "chelsa_bio01",
    "chelsa_bio04",
    "chelsa_bio12",
    "chelsa_bio15",
    "chelsa_rsds_mean",
    "chelsa_vpd_mean",
    "chelsa_sfcwind_mean",
    "chelsa_gsp",
    "chelsa_npp",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--traits-long", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--reference-environment", type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    return parser.parse_args()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty(len(p), dtype=float)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def check(checks: list[dict[str, object]], name: str, condition: bool, detail: object) -> None:
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


def demean_standardize(frame: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(frame[column], errors="coerce")
    centred = values - values.groupby(frame["taxon_name"]).transform("mean")
    return (centred / float(centred.std(ddof=0))).to_numpy(float)


def standardized_slope(y: np.ndarray, x: np.ndarray) -> float:
    yz = (y - np.mean(y)) / np.std(y, ddof=0)
    xz = (x - np.mean(x)) / np.std(x, ddof=0)
    return float(np.dot(xz, yz) / np.dot(xz, xz))


def main() -> int:
    args = parse_args()
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    inventory = pd.read_csv(args.results_dir / "v2_full27_endpoint_inventory.csv", low_memory=False)
    variance = pd.read_csv(args.results_dir / "v2_full27_variance_decomposition.csv", low_memory=False)
    geography = pd.read_csv(args.results_dir / "v2_full27_trait_geography.csv", low_memory=False)
    within = pd.read_csv(args.results_dir / "v2_full27_environment_within.csv", low_memory=False)
    among = pd.read_csv(args.results_dir / "v2_full27_environment_among.csv", low_memory=False)
    cross = pd.read_csv(args.results_dir / "v2_full27_environment_cross_scale.csv", low_memory=False)
    run_report = json.loads((args.results_dir / "v2_full27_environment_report.json").read_text(encoding="utf-8"))
    checks: list[dict[str, object]] = []

    endpoint_ids = set(contract["endpoint_id"])
    check(checks, "contract_has_27_endpoints", len(contract) == 27, len(contract))
    for label, frame in (("inventory", inventory), ("variance", variance), ("geography", geography)):
        check(
            checks,
            f"{label}_retains_exact_27",
            len(frame) == 27 and set(frame["endpoint_id"]) == endpoint_ids,
            {"rows": len(frame), "unique_endpoints": frame["endpoint_id"].nunique()},
        )
    check(
        checks,
        "measured_and_unexecuted_counts",
        inventory["measurement_status"].eq("measured").sum() == 22
        and inventory["measurement_status"].eq("unexecuted_no_measurement").sum() == 5,
        inventory["measurement_status"].value_counts().to_dict(),
    )

    expected_units = 26
    expected_rows = expected_units * len(EXPECTED_PREDICTORS)
    check(checks, "within_complete_status_grid", len(within) == expected_rows, len(within))
    scopes = sorted(among["scope"].unique().tolist())
    check(
        checks,
        "among_complete_status_grids",
        len(scopes) == 2 and all(len(part) == expected_rows for _, part in among.groupby("scope")),
        {scope: len(part) for scope, part in among.groupby("scope")},
    )
    check(checks, "cross_scale_complete_grid", len(cross) == expected_rows, len(cross))
    for label, frame in (("within", within), ("among", among)):
        check(
            checks,
            f"{label}_uses_exact_nine_predictors",
            set(frame["predictor"]) == set(EXPECTED_PREDICTORS),
            sorted(frame["predictor"].unique().tolist()),
        )
        forbidden = [
            column
            for column in frame.columns
            if column.lower() in {"evidence_grade", "evidence_tier", "abc_layer", "tier_family"}
        ]
        check(checks, f"{label}_has_no_abc_or_tier_family", not forbidden, forbidden)
        check(
            checks,
            f"{label}_retains_all_rows_regardless_of_q",
            frame["result_retained_regardless_of_q"].astype(str).str.lower().eq("true").all(),
            int(frame["result_retained_regardless_of_q"].astype(str).str.lower().eq("true").sum()),
        )

    fdr_frames = [("within", within)] + [(scope, part) for scope, part in among.groupby("scope")]
    for label, frame in fdr_frames:
        ok = frame[frame["status"].eq("ok") & frame["p_value"].notna()].copy()
        expected_q = bh_adjust(ok["p_value"].to_numpy(float)) if len(ok) else np.array([])
        actual_q = ok["q_fdr_bh_global_family"].to_numpy(float)
        check(
            checks,
            f"{label}_global_fdr_recomputed",
            np.allclose(expected_q, actual_q, rtol=1e-12, atol=1e-12),
            {"n_tests": len(ok), "max_abs_difference": float(np.max(np.abs(expected_q - actual_q))) if len(ok) else 0.0},
        )

    environment = pd.read_csv(
        args.environment,
        usecols=["obs_id", "taxon_name", *EXPECTED_PREDICTORS],
        low_memory=False,
    )
    environment["obs_id"] = environment["obs_id"].astype(str)
    check(
        checks,
        "environment_full_cohort_unique",
        len(environment) == 46276 and not environment["obs_id"].duplicated().any(),
        {"rows": len(environment), "taxa": environment["taxon_name"].nunique()},
    )
    coverage = {predictor: float(environment[predictor].notna().mean()) for predictor in EXPECTED_PREDICTORS}
    check(checks, "environment_coverage_at_least_0_98", min(coverage.values()) >= 0.98, coverage)
    if args.reference_environment is not None:
        reference = pd.read_csv(
            args.reference_environment,
            usecols=["obs_id", "taxon_name", *EXPECTED_PREDICTORS],
            low_memory=False,
        )
        reference["obs_id"] = reference["obs_id"].astype(str)
        overlap = reference.merge(
            environment,
            on="obs_id",
            how="left",
            suffixes=("_reference", "_full"),
            indicator=True,
            validate="one_to_one",
        )
        check(
            checks,
            "reference_environment_rows_all_recovered",
            overlap["_merge"].eq("both").all(),
            overlap["_merge"].value_counts().to_dict(),
        )
        variable_detail: dict[str, object] = {}
        all_equal = bool(
            overlap["taxon_name_reference"].fillna("<NA>").eq(
                overlap["taxon_name_full"].fillna("<NA>")
            ).all()
        )
        for predictor in EXPECTED_PREDICTORS:
            left = pd.to_numeric(overlap[f"{predictor}_reference"], errors="coerce")
            right = pd.to_numeric(overlap[f"{predictor}_full"], errors="coerce")
            finite = left.notna() & right.notna()
            difference = (left[finite] - right[finite]).abs()
            max_difference = float(difference.max()) if len(difference) else None
            missingness_mismatch = int((left.isna() != right.isna()).sum())
            predictor_equal = missingness_mismatch == 0 and (
                max_difference is None or max_difference <= 1e-12
            )
            all_equal = all_equal and predictor_equal
            variable_detail[predictor] = {
                "finite_overlap": int(finite.sum()),
                "maximum_absolute_difference": max_difference,
                "missingness_mismatch": missingness_mismatch,
            }
        check(
            checks,
            "reference_environment_values_exactly_reproduced",
            all_equal,
            {
                "reference_rows": len(reference),
                "taxon_name_mismatches": int(
                    (~overlap["taxon_name_reference"].fillna("<NA>").eq(
                        overlap["taxon_name_full"].fillna("<NA>")
                    )).sum()
                ),
                "variables": variable_detail,
            },
        )

    traits = pd.read_csv(
        args.traits_long,
        usecols=["obs_id", "taxon_name", "endpoint_id", "measurement_available", "value"],
        low_memory=False,
    )
    traits["obs_id"] = traits["obs_id"].astype(str)
    traits["value"] = pd.to_numeric(traits["value"], errors="coerce")
    available = traits["measurement_available"].astype(str).str.lower().isin({"true", "1", "yes"})
    traits = traits[available & traits["value"].notna()].copy()

    beta_differences: list[float] = []
    sample_within = within[
        within["status"].eq("ok") & within["inferential_unit"].eq("linear_endpoint")
    ].head(8)
    for row in sample_within.itertuples(index=False):
        part = traits[traits["endpoint_id"].eq(row.unit_id)][
            ["obs_id", "taxon_name", "value"]
        ].rename(columns={"value": row.unit_id})
        data = part.merge(environment[["obs_id", row.predictor]], on="obs_id", how="inner", validate="one_to_one")
        data = data.dropna(subset=[row.unit_id, row.predictor])
        counts = data.groupby("taxon_name").size()
        data = data[data["taxon_name"].isin(counts[counts >= 2].index)]
        y = demean_standardize(data, row.unit_id)
        x = demean_standardize(data, row.predictor)
        beta = float(np.dot(x, y) / np.dot(x, x))
        beta_differences.append(abs(beta - float(row.beta_std)))
    check(
        checks,
        "within_linear_betas_independently_recomputed",
        bool(beta_differences) and max(beta_differences) < 1e-10,
        {"n_spot_checks": len(beta_differences), "max_abs_difference": max(beta_differences) if beta_differences else None},
    )

    beta_differences = []
    primary_scope = sorted(scopes, key=lambda value: int(value.rsplit("min", 1)[1]), reverse=True)[0]
    sample_among = among[
        among["scope"].eq(primary_scope)
        & among["status"].eq("ok")
        & among["inferential_unit"].eq("linear_endpoint")
    ].head(8)
    taxon_env = environment.groupby("taxon_name")[EXPECTED_PREDICTORS].median(numeric_only=True).reset_index()
    for row in sample_among.itertuples(index=False):
        part = traits[traits["endpoint_id"].eq(row.unit_id)][["taxon_name", "value"]]
        counts = part.groupby("taxon_name").size().rename("n_trait_observations")
        medians = part.groupby("taxon_name")["value"].median().rename("trait_value")
        data = medians.to_frame().join(counts).reset_index()
        data = data[data["n_trait_observations"].ge(int(row.minimum_trait_observations_per_taxon))]
        data = data.merge(taxon_env[["taxon_name", row.predictor]], on="taxon_name", how="inner", validate="one_to_one")
        data = data.dropna(subset=["trait_value", row.predictor])
        beta = standardized_slope(data["trait_value"].to_numpy(float), data[row.predictor].to_numpy(float))
        beta_differences.append(abs(beta - float(row.beta_std)))
    check(
        checks,
        "among_linear_betas_independently_recomputed",
        bool(beta_differences) and max(beta_differences) < 1e-10,
        {"n_spot_checks": len(beta_differences), "max_abs_difference": max(beta_differences) if beta_differences else None},
    )

    check(
        checks,
        "report_confirms_no_v1_subset",
        run_report.get("old_v1_nine_endpoint_subset_used") is False,
        run_report.get("old_v1_nine_endpoint_subset_used"),
    )
    passed = all(item["passed"] for item in checks)
    output = {
        "status": "PASS" if passed else "FAIL",
        "n_checks": len(checks),
        "n_passed": sum(bool(item["passed"]) for item in checks),
        "checks": checks,
        "claim_ceiling": "exploratory trait-environment atlas; spatial and historical sensitivities are separate required gates",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

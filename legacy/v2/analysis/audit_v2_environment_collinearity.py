#!/usr/bin/env python3
"""Audit covariance among the nine fixed GEB-v2 environmental predictors.

This is a retrospective interpretation diagnostic.  It does not select
predictors, refit trait-environment models, or change the fixed multiplicity
family.  VIF is reported only to describe how poorly a hypothetical joint
nine-predictor model could separate the observed environmental gradients; the
primary atlas fits one predictor at a time.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PREDICTORS = [
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
EXPECTED_SHA256 = "1ab84254a80493776b4c435152ed3d2a1c1e68dd0e0342da0ea081eeb5cd3d9b"
EXPECTED_ROWS = 1_874
EXPECTED_TAXA = 124


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vif_table(values: pd.DataFrame) -> list[dict[str, float | str | int]]:
    standardized = (values - values.mean()) / values.std(ddof=0)
    rows: list[dict[str, float | str | int]] = []
    for predictor in PREDICTORS:
        response = standardized[predictor].to_numpy(dtype=float)
        others = [name for name in PREDICTORS if name != predictor]
        design = np.column_stack(
            [np.ones(len(standardized)), standardized[others].to_numpy(dtype=float)]
        )
        coefficients = np.linalg.lstsq(design, response, rcond=None)[0]
        residual = response - design @ coefficients
        total_ss = float(response @ response)
        residual_ss = float(residual @ residual)
        r_squared = 1.0 - residual_ss / total_ss
        vif = float("inf") if r_squared >= 1.0 else 1.0 / (1.0 - r_squared)
        rows.append(
            {
                "predictor": predictor,
                "vif": vif,
                "r_squared_against_other_eight": r_squared,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    input_hash = sha256(args.environment)
    if input_hash != EXPECTED_SHA256:
        raise SystemExit(
            f"Unexpected environment SHA-256: {input_hash}; expected {EXPECTED_SHA256}"
        )

    frame = pd.read_csv(
        args.environment,
        usecols=["taxon_name", *PREDICTORS],
        low_memory=False,
    )
    if len(frame) != EXPECTED_ROWS or frame["taxon_name"].nunique() != EXPECTED_TAXA:
        raise SystemExit(
            "Unexpected audit cohort: "
            f"rows={len(frame)}, taxa={frame['taxon_name'].nunique()}"
        )
    if frame[PREDICTORS].isna().any().any():
        raise SystemExit("The locked complete-18 overlap contains missing predictor values")

    scopes = {
        "observation_raw": frame[PREDICTORS],
        "within_taxon_demeaned": frame[PREDICTORS]
        - frame.groupby("taxon_name")[PREDICTORS].transform("mean"),
        "among_taxon_median": frame.groupby("taxon_name")[PREDICTORS].median(),
    }
    correlation_rows: list[dict[str, float | str | int]] = []
    vif_rows: list[dict[str, float | str | int]] = []
    summaries: dict[str, dict[str, float | int | str]] = {}

    for scope, values in scopes.items():
        pearson = values.corr(method="pearson")
        spearman = values.corr(method="spearman")
        for left_index, left in enumerate(PREDICTORS):
            for right in PREDICTORS[left_index + 1 :]:
                correlation_rows.append(
                    {
                        "scope": scope,
                        "n_units": len(values),
                        "predictor_1": left,
                        "predictor_2": right,
                        "pearson_r": float(pearson.loc[left, right]),
                        "spearman_rho": float(spearman.loc[left, right]),
                    }
                )
        scope_vif = vif_table(values)
        for row in scope_vif:
            vif_rows.append({"scope": scope, "n_units": len(values), **row})

        scope_pairs = [row for row in correlation_rows if row["scope"] == scope]
        largest_pair = max(scope_pairs, key=lambda row: abs(float(row["pearson_r"])))
        largest_vif = max(scope_vif, key=lambda row: float(row["vif"]))
        summaries[scope] = {
            "n_units": len(values),
            "maximum_absolute_pearson_r": abs(float(largest_pair["pearson_r"])),
            "maximum_pearson_pair": (
                f"{largest_pair['predictor_1']}__{largest_pair['predictor_2']}"
            ),
            "maximum_pair_pearson_r": float(largest_pair["pearson_r"]),
            "maximum_vif": float(largest_vif["vif"]),
            "maximum_vif_predictor": str(largest_vif["predictor"]),
        }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(correlation_rows).to_csv(
        args.out_dir / "environment_predictor_pairwise_correlations.csv",
        index=False,
        float_format="%.10g",
    )
    pd.DataFrame(vif_rows).to_csv(
        args.out_dir / "environment_predictor_vif.csv",
        index=False,
        float_format="%.10g",
    )
    report = {
        "analysis_id": "geb_v2_environment_collinearity_diagnostic_v1",
        "status": "retrospective_interpretation_diagnostic_not_model_selection",
        "source": {
            "github_actions_run": 33035785120,
            "artifact_id": 9632715852,
            "artifact_sha256": "51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6",
            "relative_path": "process_environment/complete18_chelsa_process.csv",
            "file_sha256": input_hash,
            "rows": len(frame),
            "taxa": int(frame["taxon_name"].nunique()),
            "relationship_to_primary_input": (
                "Exact all-nine-predictor overlap previously verified against 1,874 rows "
                "of the 46,276-observation primary environment table"
            ),
        },
        "predictors": PREDICTORS,
        "summaries": summaries,
        "interpretation": {
            "primary_models": "separate standardized univariate slopes",
            "vif_role": (
                "diagnostic for a hypothetical simultaneous nine-predictor model; "
                "VIF does not enter or filter the primary univariate fits"
            ),
            "claim_boundary": (
                "Environmental coefficients are marginal associations with correlated "
                "observational gradients, not effects independent of the other predictors"
            ),
        },
    }
    (args.out_dir / "environment_collinearity_report.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

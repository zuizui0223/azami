#!/usr/bin/env python3
"""Normalize, summarize and plot historical-constraint sensitivity outputs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PREDICTORS = {
    "env_chelsa_bio01_species_median": "Mean temperature",
    "env_chelsa_bio04_species_median": "Temperature seasonality",
    "env_chelsa_bio12_species_median": "Annual precipitation",
    "env_chelsa_bio15_species_median": "Precipitation seasonality",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree-audit-dir", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def short_endpoint(value: str) -> str:
    replacements = {
        "_species_median": "",
        "corolla_lab_": "corolla_",
        "orientation_angle_degrees": "orientation_angle",
        "spine_relative_length_max_proxy": "spine_length_proxy",
        "involucre_projection_roughness": "involucre_roughness",
        "involucre_spread_fraction": "involucre_spread",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value


def main() -> None:
    args = parse_args()
    tree_dir = Path(args.tree_audit_dir)
    model_dir = Path(args.model_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    audit = json.loads((tree_dir / "historical_tree_audit_report.json").read_text(encoding="utf-8"))
    coefficients = pd.read_csv(model_dir / "historical_constraint_model_coefficients.csv")
    statuses = pd.read_csv(model_dir / "historical_constraint_model_status.csv")
    random = pd.read_csv(model_dir / "historical_constraint_random_tree_summary.csv")

    # The R reference fit uses ordinary lm standard errors. Correct the output
    # label so it cannot be mistaken for the earlier Python HC3 screening model.
    coefficients["covariance_model"] = coefficients["covariance_model"].replace(
        {"OLS_HC3_reference": "OLS_standard_reference"}
    )
    statuses["covariance_model"] = statuses["covariance_model"].replace(
        {"OLS_HC3_reference": "OLS_standard_reference"}
    )
    coefficients.to_csv(out / "historical_constraint_model_coefficients_clean.csv", index=False)
    statuses.to_csv(out / "historical_constraint_model_status_clean.csv", index=False)
    random.to_csv(out / "historical_constraint_random_tree_summary.csv", index=False)

    # Figure 1: how much of the tree is direct evidence versus genus grafting.
    direct = int(audit["gbotb_lcvp"]["n_direct_backbone_tips"])
    grafted = int(audit["gbotb_lcvp"]["n_grafted_within_genus"])
    open_resolved = int(audit["opentree"]["n_resolved"])
    total = int(audit["n_input_taxa"])
    coverage = pd.DataFrame({
        "resource": ["GBOTB direct tips", "GBOTB genus-grafted", "OpenTree TNRS resolved"],
        "n_taxa": [direct, grafted, open_resolved],
        "fraction": [direct / total, grafted / total, open_resolved / total],
    })
    coverage.to_csv(out / "historical_tree_taxon_coverage.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(coverage["resource"], coverage["fraction"] * 100)
    ax.set_ylabel("Taxa (%)")
    ax.set_ylim(0, 105)
    ax.tick_params(axis="x", rotation=20)
    for index, row in coverage.iterrows():
        ax.text(index, row["fraction"] * 100 + 2, f"{int(row['n_taxa'])}/{total}", ha="center")
    fig.tight_layout()
    fig.savefig(out / "figure_historical_tree_coverage.png", dpi=300)
    fig.savefig(out / "figure_historical_tree_coverage.pdf")
    plt.close(fig)

    # Figure 2: deterministic Pagel-lambda PGLS coefficients for main traits.
    deterministic = coefficients.loc[
        coefficients["analysis_tier"].eq("main")
        & coefficients["covariance_model"].eq("PGLS_Pagel_lambda")
        & coefficients["scenario"].isin(["GBOTB_LCVP_S1", "GBOTB_LCVP_S3"])
    ].copy()
    deterministic["label"] = (
        deterministic["endpoint"].map(short_endpoint)
        + " — "
        + deterministic["predictor"].map(PREDICTORS)
        + " — "
        + deterministic["scenario"].str.replace("GBOTB_LCVP_", "", regex=False)
    )
    deterministic = deterministic.sort_values(["endpoint", "predictor", "scenario"]).reset_index(drop=True)
    y = np.arange(len(deterministic))
    fig, ax = plt.subplots(figsize=(10.0, max(7.0, len(deterministic) * 0.24 + 1.5)))
    ax.errorbar(
        deterministic["estimate_standardized"],
        y,
        xerr=[
            deterministic["estimate_standardized"] - deterministic["ci95_low"],
            deterministic["ci95_high"] - deterministic["estimate_standardized"],
        ],
        fmt="o",
        capsize=2,
    )
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y, deterministic["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Standardized PGLS coefficient (95% CI)")
    ax.set_title("Historical sensitivity: deterministic GBOTB grafting scenarios")
    fig.tight_layout()
    fig.savefig(out / "figure_pgls_deterministic_scenarios.png", dpi=300)
    fig.savefig(out / "figure_pgls_deterministic_scenarios.pdf")
    plt.close(fig)

    # Figure 3: coefficient ranges across randomized genus graftings.
    random["label"] = random["endpoint"].map(short_endpoint) + " — " + random["predictor"].map(PREDICTORS)
    random = random.sort_values(["historical_sensitivity_robust_candidate", "sign_stability"], ascending=[False, False]).reset_index(drop=True)
    y = np.arange(len(random))
    fig, ax = plt.subplots(figsize=(10.0, max(7.0, len(random) * 0.28 + 1.5)))
    ax.errorbar(
        random["estimate_median"],
        y,
        xerr=[
            random["estimate_median"] - random["estimate_q025"],
            random["estimate_q975"] - random["estimate_median"],
        ],
        fmt="o",
        capsize=2,
    )
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y, random["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Median standardized coefficient across randomized trees (95% tree range)")
    ax.set_title("Historical sensitivity: uncertainty in genus-level grafting")
    fig.tight_layout()
    fig.savefig(out / "figure_pgls_random_tree_ranges.png", dpi=300)
    fig.savefig(out / "figure_pgls_random_tree_ranges.pdf")
    plt.close(fig)

    robust = random.loc[random["historical_sensitivity_robust_candidate"].astype(bool)].copy()
    robust.to_csv(out / "historical_sensitivity_robust_candidates.csv", index=False)
    report = {
        "n_taxa": total,
        "n_direct_gbotb_tips": direct,
        "n_gbotb_genus_grafted": grafted,
        "n_opentree_resolved": open_resolved,
        "n_successful_model_fits": int(statuses["status"].eq("success").sum()),
        "n_failed_model_fits": int(statuses["status"].eq("failed").sum()),
        "n_random_tree_comparisons": int(len(random)),
        "n_robust_candidates_across_random_graftings": int(len(robust)),
        "semantic_status": (
            "Historical-constraint sensitivity only. Robust candidates are stable across randomized "
            "GBOTB genus graftings, but the analysis does not resolve hybridization, chloroplast capture "
            "or allopolyploid parentage and does not establish one true Cirsium species tree."
        ),
    }
    (out / "historical_sensitivity_final_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

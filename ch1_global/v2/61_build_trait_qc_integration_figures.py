#!/usr/bin/env python3
"""Build submission-ready QC/coverage figures for primary and auxiliary traits."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--primary-qc", required=True)
    p.add_argument("--environment-gee", required=True)
    p.add_argument("--involucre-status", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def save_primary_retention(table: pd.DataFrame, out: Path) -> None:
    ordered = table.set_index("trait").loc[["colour", "shape", "orientation"]].reset_index()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.bar(ordered["trait"], ordered["usable_fraction"] * 100)
    ax.set_ylabel("Heads passing QC (%)")
    ax.set_xlabel("Primary continuous trait")
    ax.set_ylim(0, 100)
    for index, row in ordered.iterrows():
        ax.text(index, row["usable_fraction"] * 100 + 2, f"{row['usable_fraction'] * 100:.1f}%", ha="center")
    fig.tight_layout()
    fig.savefig(out / "figure_primary_trait_qc_retention.png", dpi=300)
    fig.savefig(out / "figure_primary_trait_qc_retention.pdf")
    plt.close(fig)


def save_environment_bias(table: pd.DataFrame, out: Path) -> None:
    labels = {
        "env_chelsa_bio01_native": "Annual mean temperature",
        "env_chelsa_bio04_native": "Temperature seasonality",
        "env_chelsa_bio12_native": "Annual precipitation",
        "env_chelsa_bio15_native": "Precipitation seasonality",
    }
    frame = table.copy()
    frame["label"] = frame["trait"] + " — " + frame["predictor"].map(labels)
    frame = frame.sort_values(["trait", "predictor"]).reset_index(drop=True)
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(8.5, 6.4))
    ax.errorbar(
        frame["odds_ratio_per_predictor_sd"], y,
        xerr=[
            frame["odds_ratio_per_predictor_sd"] - frame["ci95_low_or"],
            frame["ci95_high_or"] - frame["odds_ratio_per_predictor_sd"],
        ],
        fmt="o", capsize=3,
    )
    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_yticks(y, frame["label"])
    ax.set_xlabel("Odds ratio for QC retention per predictor SD")
    ax.set_ylabel("")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(out / "figure_qc_environment_bias.png", dpi=300)
    fig.savefig(out / "figure_qc_environment_bias.pdf")
    plt.close(fig)


def save_auxiliary_coverage(table: pd.DataFrame, out: Path) -> None:
    frame = table.sort_values("n_heads", ascending=False).copy()
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.bar(frame["status"].astype(str), frame["fraction"] * 100)
    ax.set_ylabel("High-resolution heads (%)")
    ax.set_xlabel("Auxiliary involucre measurement status")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(out / "figure_involucre_auxiliary_qc_status.png", dpi=300)
    fig.savefig(out / "figure_involucre_auxiliary_qc_status.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    primary = pd.read_csv(args.primary_qc)
    gee = pd.read_csv(args.environment_gee)
    involucre = pd.read_csv(args.involucre_status)
    save_primary_retention(primary, out)
    save_environment_bias(gee, out)
    save_auxiliary_coverage(involucre, out)
    print({"n_figures": 3, "out_dir": str(out)})


if __name__ == "__main__":
    main()

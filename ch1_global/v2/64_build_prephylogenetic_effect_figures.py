#!/usr/bin/env python3
"""Plot standardized pre-phylogenetic continuous environment coefficients."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PREDICTORS = {
    "env_chelsa_bio01_native": "Mean temperature",
    "env_chelsa_bio04_native": "Temperature seasonality",
    "env_chelsa_bio12_native": "Annual precipitation",
    "env_chelsa_bio15_native": "Precipitation seasonality",
    "env_chelsa_bio01_species_median": "Mean temperature",
    "env_chelsa_bio04_species_median": "Temperature seasonality",
    "env_chelsa_bio12_species_median": "Annual precipitation",
    "env_chelsa_bio15_species_median": "Precipitation seasonality",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--coefficients", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def short_endpoint(value: str) -> str:
    replacements = {
        "_observation_median": "",
        "_species_median": "",
        "_median": "",
        "corolla_lab_": "corolla_",
        "orientation_angle_degrees": "orientation_angle",
        "spine_relative_length_max_proxy": "spine_length_proxy",
        "involucre_projection_roughness": "involucre_roughness",
        "involucre_spread_fraction": "involucre_spread",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value


def forest(frame: pd.DataFrame, title: str, filename: str, out: Path) -> None:
    frame = frame.sort_values(["trait_group", "endpoint", "predictor"]).reset_index(drop=True)
    if frame.empty:
        return
    frame["label"] = frame["endpoint"].map(short_endpoint) + " — " + frame["predictor"].map(PREDICTORS)
    y = np.arange(len(frame))
    height = max(5.0, 0.29 * len(frame) + 1.6)
    fig, ax = plt.subplots(figsize=(9.4, height))
    ax.errorbar(
        frame["estimate_standardized"], y,
        xerr=[
            frame["estimate_standardized"] - frame["ci95_low"],
            frame["ci95_high"] - frame["estimate_standardized"],
        ],
        fmt="o", capsize=2.5,
    )
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_yticks(y, frame["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Standardized coefficient (95% CI)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out / f"{filename}.png", dpi=300)
    fig.savefig(out / f"{filename}.pdf")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    table = pd.read_csv(args.coefficients)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    selections = [
        (
            table["analysis_tier"].eq("main")
            & table["model_scale"].eq("observation_within_species")
            & table["cohort"].eq("all_coordinate_eligible"),
            "Main traits: within-species environmental screening",
            "figure_main_within_species_all",
        ),
        (
            table["analysis_tier"].eq("main")
            & table["model_scale"].eq("observation_within_species")
            & table["cohort"].eq("positional_accuracy_le_10km"),
            "Main traits: within-species screening, positional accuracy ≤10 km",
            "figure_main_within_species_le10km",
        ),
        (
            table["analysis_tier"].eq("main")
            & table["model_scale"].eq("species_between"),
            "Main traits: between-species pre-PGLS screening",
            "figure_main_between_species_prephylogenetic",
        ),
        (
            table["analysis_tier"].eq("auxiliary")
            & table["model_scale"].eq("observation_within_species")
            & table["cohort"].eq("positional_accuracy_le_10km"),
            "Auxiliary involucre proxies: within-species screening, ≤10 km",
            "figure_auxiliary_within_species_le10km",
        ),
        (
            table["analysis_tier"].eq("auxiliary")
            & table["model_scale"].eq("species_between"),
            "Auxiliary involucre proxies: between-species pre-PGLS screening",
            "figure_auxiliary_between_species_prephylogenetic",
        ),
    ]
    created = 0
    for mask, title, filename in selections:
        frame = table.loc[mask].copy()
        if len(frame):
            forest(frame, title, filename, out)
            created += 1
    print({"n_figures": created, "out_dir": str(out)})


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Join primary continuous traits, auxiliary involucre proxies and environment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--primary-observation", required=True)
    p.add_argument("--primary-species", required=True)
    p.add_argument("--auxiliary-observation", required=True)
    p.add_argument("--auxiliary-species", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    primary_observation = pd.read_csv(args.primary_observation)
    primary_species = pd.read_csv(args.primary_species)
    auxiliary_observation = pd.read_csv(args.auxiliary_observation)
    auxiliary_species = pd.read_csv(args.auxiliary_species)
    for frame in (primary_observation, auxiliary_observation):
        frame["obs_id"] = frame["obs_id"].astype(str)
    if primary_observation["obs_id"].duplicated().any() or auxiliary_observation["obs_id"].duplicated().any():
        raise ValueError("Observation inputs must be unique by obs_id")
    if primary_species["taxon_name"].duplicated().any() or auxiliary_species["taxon_name"].duplicated().any():
        raise ValueError("Species inputs must be unique by taxon_name")

    auxiliary_observation_payload = auxiliary_observation.drop(
        columns=["taxon_name", "latitude", "longitude"], errors="ignore"
    )
    observation = primary_observation.merge(
        auxiliary_observation_payload, on="obs_id", how="left", validate="one_to_one"
    )
    species = primary_species.merge(
        auxiliary_species, on="taxon_name", how="left", validate="one_to_one"
    )
    if len(observation) != len(primary_observation) or len(species) != len(primary_species):
        raise ValueError("Integrated joins lost primary observations or taxa")
    observation.to_csv(out / "integrated_primary_auxiliary_observation.csv", index=False, encoding="utf-8-sig")
    species.to_csv(out / "integrated_primary_auxiliary_species.csv", index=False, encoding="utf-8-sig")

    endpoints = pd.DataFrame([
        ("main", "orientation", "orientation_angle_degrees_median", "orientation_angle_degrees_species_median", "degrees", "signed angle: 0 up, 90 horizontal, 180 down"),
        ("main", "colour", "corolla_lab_lightness_median", "corolla_lab_lightness_species_median", "Lab L*", "corolla lightness"),
        ("main", "colour", "corolla_lab_chroma_median", "corolla_lab_chroma_species_median", "Lab chroma", "corolla chroma"),
        ("main", "colour", "corolla_hue_sin_median", "corolla_hue_sin_species_median", "unitless", "circular hue sine"),
        ("main", "colour", "corolla_hue_cos_median", "corolla_hue_cos_species_median", "unitless", "circular hue cosine"),
        ("main", "shape", "shape_aspect_ratio_median", "shape_aspect_ratio_species_median", "ratio", "capitulum outline aspect ratio"),
        ("main", "shape", "shape_circularity_median", "shape_circularity_species_median", "unitless", "capitulum outline circularity"),
        ("main", "shape", "shape_solidity_median", "shape_solidity_species_median", "unitless", "capitulum outline solidity"),
        ("main", "shape", "shape_width_cv_median", "shape_width_cv_species_median", "unitless", "capitulum width-profile variation"),
        ("auxiliary", "involucre", "involucre_projection_roughness_observation_median", "involucre_projection_roughness_species_median", "relative", "exploratory contour projection roughness"),
        ("auxiliary", "involucre", "involucre_spread_fraction_observation_median", "involucre_spread_fraction_species_median", "fraction", "exploratory positive-projection contour fraction"),
        ("auxiliary", "spine", "spine_relative_length_max_proxy_observation_median", "spine_relative_length_max_proxy_species_median", "relative", "exploratory maximum spine-like contour projection"),
    ], columns=[
        "analysis_tier", "trait_group", "observation_variable", "species_variable",
        "unit", "interpretation",
    ])
    endpoints["phylogenetic_status"] = endpoints["analysis_tier"].map({
        "main": "species_between_requires_later_PGLS_sensitivity",
        "auxiliary": "exploratory_species_between_requires_later_PGLS_sensitivity",
    })
    endpoints.to_csv(out / "integrated_continuous_endpoint_registry.csv", index=False)

    report = {
        "n_observations": int(len(observation)),
        "n_species": int(len(species)),
        "n_observations_with_auxiliary": int(observation["n_usable_heads_observation"].fillna(0).gt(0).sum()),
        "n_species_with_auxiliary": int(species["n_usable_heads_species"].fillna(0).gt(0).sum()),
        "n_main_endpoints": int(endpoints["analysis_tier"].eq("main").sum()),
        "n_auxiliary_endpoints": int(endpoints["analysis_tier"].eq("auxiliary").sum()),
        "semantic_status": (
            "Integrated observation/species tables. Main and auxiliary endpoints remain tiered; "
            "species-between results are not phylogenetically corrected at this stage."
        ),
    }
    (out / "integrated_continuous_join_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

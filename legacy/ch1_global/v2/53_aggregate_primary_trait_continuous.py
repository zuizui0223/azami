#!/usr/bin/env python3
"""Aggregate continuous head measurements to observation and species levels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

STATUS_COLUMNS = ("colour_status", "shape_status", "orientation_status")
METRICS = {
    "colour": [
        "corolla_visible_fraction", "corolla_lab_lightness", "corolla_lab_chroma",
        "corolla_hue_sin", "corolla_hue_cos", "corolla_white_fraction",
        "corolla_redmagenta_fraction", "corolla_purple_fraction", "corolla_yellow_fraction",
    ],
    "shape": [
        "shape_aspect_ratio", "shape_circularity", "shape_solidity",
        "shape_width_cv", "shape_narrow_end_width", "shape_middle_width",
    ],
    "orientation": ["orientation_angle_degrees"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--head-measurements", required=True)
    p.add_argument("--head-metadata", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-heads-per-observation", type=int, default=1)
    p.add_argument("--min-observations-per-species", type=int, default=3)
    return p.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def mad(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if not len(numeric):
        return float("nan")
    median = np.median(numeric)
    return float(np.median(np.abs(numeric - median)))


def aggregate_group(
    frame: pd.DataFrame, group_columns: list[str], unit_name: str
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(group_columns, dropna=False, sort=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = {column: value for column, value in zip(group_columns, keys)}
        row[f"n_heads_{unit_name}"] = int(group["annotation_unit_id"].nunique())
        for trait, columns in METRICS.items():
            usable = group.loc[group[f"{trait}_status"].eq("usable")]
            row[f"n_usable_heads_{trait}"] = int(
                usable["annotation_unit_id"].nunique()
            )
            for column in columns:
                values = pd.to_numeric(usable[column], errors="coerce")
                row[f"{column}_median"] = (
                    float(values.median()) if values.notna().any() else float("nan")
                )
                row[f"{column}_mad"] = mad(values)
        hue_sin = row.get("corolla_hue_sin_median", float("nan"))
        hue_cos = row.get("corolla_hue_cos_median", float("nan"))
        if np.isfinite(hue_sin) and np.isfinite(hue_cos):
            row["corolla_hue_degrees_circular"] = float(
                (np.rad2deg(np.arctan2(hue_sin, hue_cos)) + 360.0) % 360.0
            )
        else:
            row["corolla_hue_degrees_circular"] = float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    heads = pd.read_csv(args.head_measurements, dtype=str, keep_default_na=False)
    metadata = pd.read_csv(args.head_metadata, dtype=str, keep_default_na=False)
    required_heads = {"annotation_unit_id", *STATUS_COLUMNS}
    required_metadata = {
        "annotation_unit_id", "obs_id", "photo_id", "taxon_name", "latitude", "longitude",
    }
    missing_heads = required_heads.difference(heads.columns)
    missing_metadata = required_metadata.difference(metadata.columns)
    if missing_heads or missing_metadata:
        raise ValueError(
            f"Missing columns: measurements={sorted(missing_heads)} "
            f"metadata={sorted(missing_metadata)}"
        )
    if heads.empty or heads["annotation_unit_id"].duplicated().any():
        raise ValueError(
            "Head measurements must contain one row per annotation_unit_id"
        )

    metadata = metadata.drop_duplicates("annotation_unit_id").copy()
    if metadata["annotation_unit_id"].duplicated().any():
        raise ValueError("Head metadata are not unique after deduplication")
    merged = heads.merge(
        metadata[[
            "annotation_unit_id", "obs_id", "photo_id", "taxon_name",
            "latitude", "longitude",
        ]],
        on="annotation_unit_id",
        how="left",
        validate="one_to_one",
    )
    if merged["obs_id"].map(text).eq("").any():
        raise ValueError("Some continuous measurements do not match head metadata")

    observation = aggregate_group(
        merged,
        ["obs_id", "photo_id", "taxon_name", "latitude", "longitude"],
        "per_observation",
    )
    observation = observation.loc[
        observation["n_heads_per_observation"].ge(args.min_heads_per_observation)
    ].copy()

    species_rows: list[dict[str, Any]] = []
    observation_metric_columns = [
        column
        for column in observation.columns
        if column.endswith("_median") and not column.startswith("n_")
    ]
    for taxon_name, group in observation.groupby("taxon_name", sort=True):
        row: dict[str, Any] = {
            "taxon_name": taxon_name,
            "n_observations_per_species": int(group["obs_id"].nunique()),
        }
        for trait in METRICS:
            row[f"n_observations_usable_{trait}"] = int(
                group[f"n_usable_heads_{trait}"].gt(0).sum()
            )
        for column in observation_metric_columns:
            values = pd.to_numeric(group[column], errors="coerce")
            base = column[:-7]
            row[f"{base}_species_median"] = (
                float(values.median()) if values.notna().any() else float("nan")
            )
            row[f"{base}_species_mad"] = mad(values)
        hue_sin = row.get("corolla_hue_sin_species_median", float("nan"))
        hue_cos = row.get("corolla_hue_cos_species_median", float("nan"))
        row["corolla_hue_degrees_species_circular"] = (
            float((np.rad2deg(np.arctan2(hue_sin, hue_cos)) + 360.0) % 360.0)
            if np.isfinite(hue_sin) and np.isfinite(hue_cos)
            else float("nan")
        )
        species_rows.append(row)
    species = pd.DataFrame(species_rows)
    species = species.loc[
        species["n_observations_per_species"].ge(
            args.min_observations_per_species
        )
    ].copy()

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    merged.to_csv(
        output / "primary_trait_continuous_head_with_metadata.csv",
        index=False,
        encoding="utf-8-sig",
    )
    observation.to_csv(
        output / "primary_trait_continuous_observation_level.csv",
        index=False,
        encoding="utf-8-sig",
    )
    species.to_csv(
        output / "primary_trait_continuous_species_level.csv",
        index=False,
        encoding="utf-8-sig",
    )

    report = {
        "n_head_measurements": int(len(heads)),
        "n_observations": int(len(observation)),
        "n_species": int(len(species)),
        "min_heads_per_observation": args.min_heads_per_observation,
        "min_observations_per_species": args.min_observations_per_species,
        "semantic_status": (
            "Observation/species summaries of automated continuous image "
            "measurements. Medians and MAD retain within-unit dispersion."
        ),
    }
    (
        output / "primary_trait_continuous_aggregation_report.json"
    ).write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Merge exhaustive continuous-trait shard outputs and build nested analysis units.

The merged head table is immutable source data. Photo- and observation-level
summaries are derived with medians across usable heads, preserving counts and QC
failure rates. No species or spatial thinning occurs here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

TRAITS = {
    "orientation_angle_degrees": "orientation_status",
    "corolla_lab_lightness": "colour_status",
    "corolla_lab_chroma": "colour_status",
    "corolla_hue_sin": "colour_status",
    "corolla_hue_cos": "colour_status",
    "shape_aspect_ratio": "shape_status",
    "shape_circularity": "shape_status",
    "shape_solidity": "shape_status",
    "shape_width_cv": "shape_status",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard-manifest", required=True)
    parser.add_argument("--shard-results-root", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_exact(root: Path, filename: str) -> list[Path]:
    return sorted(path for path in root.rglob(filename) if path.is_file())


def numeric_median(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    return float(numeric.median()) if numeric.notna().any() else np.nan


def aggregate_level(heads: pd.DataFrame, keys: list[str], level_name: str) -> pd.DataFrame:
    identity_candidates = [
        "taxon_name", "latitude", "longitude", "positional_accuracy", "observed_on",
        "coordinate_usable_for_environment", "geoprivacy", "obscured",
        "spatial_block_1deg", "spatial_block_2deg", "spatial_block_5deg",
    ]
    rows: list[dict[str, Any]] = []
    for key_values, group in heads.groupby(keys, sort=True, dropna=False):
        if not isinstance(key_values, tuple):
            key_values = (key_values,)
        row = dict(zip(keys, key_values))
        row[f"n_heads_{level_name}"] = int(len(group))
        row[f"n_photos_{level_name}"] = int(group["photo_id"].astype(str).nunique())
        for column in identity_candidates:
            if column in group.columns:
                nonempty = group[column].dropna().astype(str)
                nonempty = nonempty.loc[nonempty.ne("")]
                row[column] = nonempty.iloc[0] if not nonempty.empty else ""
        for trait, status in TRAITS.items():
            if trait not in group.columns or status not in group.columns:
                row[f"{trait}_median"] = np.nan
                row[f"{trait}_n_usable_heads"] = 0
                continue
            usable = group.loc[group[status].astype(str).eq("usable"), trait]
            row[f"{trait}_median"] = numeric_median(usable)
            row[f"{trait}_n_usable_heads"] = int(pd.to_numeric(usable, errors="coerce").notna().sum())
        for status in {value for value in TRAITS.values()}:
            if status in group.columns:
                stem = status.removesuffix("_status")
                row[f"{stem}_usable_fraction"] = float(group[status].astype(str).eq("usable").mean())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.shard_manifest)
    results_root = Path(args.shard_results_root)
    if not manifest_path.is_file() or not results_root.is_dir():
        raise FileNotFoundError("Shard manifest or result root is missing")
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    queue_paths = find_exact(results_root, "queue_shard.csv")
    crop_paths = find_exact(results_root, "yolo_crop_metadata.csv")
    trait_paths = find_exact(results_root, "primary_trait_continuous_head_measurements.csv")
    summary_paths = find_exact(results_root, "shard_run_summary.json")
    expected = len(manifest)
    if not (len(queue_paths) == len(crop_paths) == len(trait_paths) == len(summary_paths) == expected):
        raise ValueError(
            f"Incomplete shard collection: expected={expected}, queues={len(queue_paths)}, "
            f"crops={len(crop_paths)}, traits={len(trait_paths)}, summaries={len(summary_paths)}"
        )

    queue = pd.concat([pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False) for path in queue_paths], ignore_index=True)
    crops = pd.concat([pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False) for path in crop_paths], ignore_index=True)
    traits = pd.concat([pd.read_csv(path, low_memory=False) for path in trait_paths], ignore_index=True)
    if queue["queue_id"].duplicated().any() or queue["photo_id"].duplicated().any():
        raise ValueError("Merged exhaustive queue overlaps in queue_id or photo_id")
    if crops.duplicated(["audit_id", "det_index"]).any():
        raise ValueError("Merged detector crops overlap")
    if traits["annotation_unit_id"].astype(str).duplicated().any():
        raise ValueError("Merged continuous measurements overlap")

    expected_photos = int(pd.to_numeric(manifest["n_photos"], errors="raise").sum())
    if len(queue) != expected_photos:
        raise ValueError(f"Queue coverage differs from manifest: {len(queue)} != {expected_photos}")

    metadata_columns = [column for column in queue.columns if column not in traits.columns or column in {"queue_id", "obs_id", "photo_id"}]
    head = traits.merge(
        crops.drop_duplicates("annotation_unit_id") if "annotation_unit_id" in crops.columns else crops,
        on="annotation_unit_id",
        how="left",
        suffixes=("", "_crop"),
        validate="one_to_one",
    )
    join_key = "queue_id" if "queue_id" in head.columns else "audit_id"
    if join_key != "queue_id":
        head["queue_id"] = head[join_key].astype(str)
    head = head.merge(queue[metadata_columns].drop_duplicates("queue_id"), on="queue_id", how="left", suffixes=("", "_queue"), validate="many_to_one")
    if head["photo_id"].astype(str).eq("").any() or head["obs_id"].astype(str).eq("").any():
        raise ValueError("At least one measured head could not be joined to photo/observation metadata")

    photo = aggregate_level(head, ["photo_id"], "photo")
    photo_to_obs = queue[["photo_id", "obs_id"]].drop_duplicates("photo_id")
    photo = photo.merge(photo_to_obs, on="photo_id", how="left", validate="one_to_one")
    observation = aggregate_level(head, ["obs_id"], "observation")

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    queue.sort_values(["taxon_name", "obs_id", "photo_id"]).to_csv(output / "exhaustive_photo_queue_merged.csv", index=False, encoding="utf-8-sig")
    crops.sort_values(["audit_id", "det_index"]).to_csv(output / "exhaustive_yolo_crop_metadata.csv", index=False, encoding="utf-8-sig")
    head.sort_values(["taxon_name", "obs_id", "photo_id", "annotation_unit_id"]).to_csv(output / "exhaustive_continuous_head_level.csv", index=False, encoding="utf-8-sig")
    photo.sort_values("photo_id").to_csv(output / "exhaustive_continuous_photo_level.csv", index=False, encoding="utf-8-sig")
    observation.sort_values("obs_id").to_csv(output / "exhaustive_continuous_observation_level.csv", index=False, encoding="utf-8-sig")

    shard_summaries = [json.loads(path.read_text(encoding="utf-8")) for path in summary_paths]
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_shards": expected,
        "n_queue_photos": int(len(queue)),
        "n_queue_observations": int(queue["obs_id"].nunique()),
        "n_queue_species": int(queue["taxon_name"].nunique()),
        "n_photos_with_detected_heads": int(head["photo_id"].nunique()),
        "n_observations_with_detected_heads": int(head["obs_id"].nunique()),
        "n_detected_heads": int(len(head)),
        "n_colour_usable": int(head["colour_status"].astype(str).eq("usable").sum()),
        "n_shape_usable": int(head["shape_status"].astype(str).eq("usable").sum()),
        "n_orientation_usable": int(head["orientation_status"].astype(str).eq("usable").sum()),
        "source_sha256": {
            "shard_manifest": sha256_file(manifest_path),
            "queues": {str(path.relative_to(results_root)): sha256_file(path) for path in queue_paths},
            "traits": {str(path.relative_to(results_root)): sha256_file(path) for path in trait_paths},
        },
        "shard_summaries": shard_summaries,
        "semantic_boundary": "All eligible photos were predicted before any analysis thinning. Derived cohorts must never replace these source tables.",
    }
    (output / "exhaustive_continuous_merge_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

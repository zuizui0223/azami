#!/usr/bin/env python3
"""Merge exhaustive all-photo predictions without pre-prediction thinning.

All queued photos remain in the detection-audit tables. Photos with no detected
capitulum (including leaf-only photographs) are excluded from continuous-trait
aggregation but are never silently deleted. Head, photo and observation trait
tables contain only detector-positive images.
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


def read_concat(paths: list[Path], **kwargs: Any) -> pd.DataFrame:
    frames = [pd.read_csv(path, **kwargs) for path in paths]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def numeric_median(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    return float(numeric.median()) if numeric.notna().any() else np.nan


def aggregate_level(heads: pd.DataFrame, keys: list[str], level_name: str) -> pd.DataFrame:
    identity = [
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
        for column in identity:
            if column in group.columns:
                values = group[column].dropna().astype(str)
                values = values.loc[values.ne("")]
                row[column] = values.iloc[0] if not values.empty else ""
        for trait, status in TRAITS.items():
            usable = group.loc[group[status].astype(str).eq("usable"), trait]
            row[f"{trait}_median"] = numeric_median(usable)
            row[f"{trait}_n_usable_heads"] = int(pd.to_numeric(usable, errors="coerce").notna().sum())
        for status in sorted(set(TRAITS.values())):
            stem = status.removesuffix("_status")
            row[f"{stem}_usable_fraction"] = float(group[status].astype(str).eq("usable").mean())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.shard_manifest)
    results_root = Path(args.shard_results_root)
    if not manifest_path.is_file() or not results_root.is_dir():
        raise FileNotFoundError("Shard manifest or results root is missing")

    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    expected = len(manifest)
    paths = {
        "queue": find_exact(results_root, "queue_shard.csv"),
        "screen": find_exact(results_root, "screening_results.csv"),
        "crop": find_exact(results_root, "yolo_crop_metadata.csv"),
        "key": find_exact(results_root, "trait_annotation_key_private.csv"),
        "trait": find_exact(results_root, "primary_trait_continuous_head_measurements.csv"),
        "summary": find_exact(results_root, "shard_run_summary.json"),
    }
    for label in ("queue", "screen", "summary"):
        if len(paths[label]) != expected:
            raise ValueError(f"Expected {expected} {label} files, found {len(paths[label])}")
    # Detector-positive files can be absent only when an entire shard has no heads.
    if not (len(paths["crop"]) == len(paths["key"]) == len(paths["trait"])):
        raise ValueError("Crop, private-key and continuous-trait shard counts differ")

    queue = read_concat(paths["queue"], dtype=str, keep_default_na=False, low_memory=False)
    screens = read_concat(paths["screen"], dtype=str, keep_default_na=False, low_memory=False)
    crops = read_concat(paths["crop"], dtype=str, keep_default_na=False, low_memory=False)
    keys = read_concat(paths["key"], dtype=str, keep_default_na=False, low_memory=False)
    traits = read_concat(paths["trait"], low_memory=False)

    if queue.empty or queue["queue_id"].duplicated().any() or queue["photo_id"].duplicated().any():
        raise ValueError("Merged exhaustive queue must contain unique photo rows")
    if screens.empty or screens["queue_id"].duplicated().any():
        raise ValueError("Every photo must have exactly one terminal screening result")
    if set(queue["queue_id"]) != set(screens["queue_id"]):
        missing = len(set(queue["queue_id"]).difference(set(screens["queue_id"])))
        extra = len(set(screens["queue_id"]).difference(set(queue["queue_id"])))
        raise ValueError(f"Screening coverage mismatch: missing={missing}, extra={extra}")
    allowed_screen = {"detected", "no_detection", "missing_image", "error"}
    unexpected = sorted(set(screens["screen_status"]).difference(allowed_screen))
    if unexpected:
        raise ValueError(f"Unexpected screening statuses: {unexpected}")

    expected_photos = int(pd.to_numeric(manifest["n_photos"], errors="raise").sum())
    if len(queue) != expected_photos:
        raise ValueError(f"Queue coverage differs from manifest: {len(queue)} != {expected_photos}")

    detected_ids = set(screens.loc[screens["screen_status"].eq("detected"), "queue_id"])
    if detected_ids:
        if crops.empty or keys.empty or traits.empty:
            raise ValueError("Detected photos exist but detector/key/trait tables are missing")
        if crops.duplicated(["queue_id", "det_index"]).any():
            raise ValueError("Merged detector crops overlap")
        if keys["annotation_unit_id"].duplicated().any() or traits["annotation_unit_id"].astype(str).duplicated().any():
            raise ValueError("Merged private keys or continuous measurements overlap")
        if set(crops["queue_id"]).difference(detected_ids):
            raise ValueError("Crop metadata references non-detected photos")
        head = traits.merge(keys, on="annotation_unit_id", how="left", validate="one_to_one")
        head = head.merge(
            crops,
            on=["queue_id", "audit_id", "photo_id", "det_index"],
            how="left",
            suffixes=("", "_crop"),
            validate="one_to_one",
        )
        metadata = queue.drop_duplicates("queue_id")
        duplicate_columns = [column for column in metadata.columns if column in head.columns and column != "queue_id"]
        metadata = metadata.drop(columns=duplicate_columns)
        head = head.merge(metadata, on="queue_id", how="left", validate="many_to_one")
        if head["obs_id"].astype(str).eq("").any() or head["taxon_name"].astype(str).eq("").any():
            raise ValueError("At least one measured head could not be joined to observation metadata")
        photo = aggregate_level(head, ["photo_id"], "photo")
        photo = photo.merge(queue[["photo_id", "obs_id"]].drop_duplicates("photo_id"), on="photo_id", how="left", validate="one_to_one")
        observation = aggregate_level(head, ["obs_id"], "observation")
    else:
        head = pd.DataFrame()
        photo = pd.DataFrame()
        observation = pd.DataFrame()

    photo_audit = queue.merge(
        screens[["queue_id", "screen_status", "n_detections", "max_yolo_conf", "message"]],
        on="queue_id",
        how="left",
        validate="one_to_one",
    )
    photo_audit["eligible_for_trait_analysis"] = photo_audit["screen_status"].eq("detected")
    photo_audit["exclusion_reason"] = np.where(
        photo_audit["screen_status"].eq("detected"),
        "",
        np.where(photo_audit["screen_status"].eq("no_detection"), "no_visible_capitulum_detected", photo_audit["screen_status"]),
    )
    species_audit = (
        photo_audit.groupby("taxon_name", sort=True)
        .agg(
            n_all_photos=("photo_id", "size"),
            n_all_observations=("obs_id", "nunique"),
            n_detected_photos=("eligible_for_trait_analysis", "sum"),
            n_detected_observations=("obs_id", lambda values: values[photo_audit.loc[values.index, "eligible_for_trait_analysis"]].nunique()),
            n_no_detection_photos=("screen_status", lambda values: values.eq("no_detection").sum()),
            n_error_photos=("screen_status", lambda values: values.isin(["error", "missing_image"]).sum()),
        )
        .reset_index()
    )
    species_audit["detected_photo_fraction"] = species_audit["n_detected_photos"] / species_audit["n_all_photos"]

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    queue.sort_values(["taxon_name", "obs_id", "photo_id"]).to_csv(output / "exhaustive_photo_queue_merged.csv", index=False, encoding="utf-8-sig")
    photo_audit.sort_values(["taxon_name", "obs_id", "photo_id"]).to_csv(output / "exhaustive_photo_detection_audit.csv", index=False, encoding="utf-8-sig")
    species_audit.to_csv(output / "exhaustive_species_detection_audit.csv", index=False, encoding="utf-8-sig")
    screens.sort_values("queue_id").to_csv(output / "exhaustive_screening_results.csv", index=False, encoding="utf-8-sig")
    if not head.empty:
        crops.sort_values(["queue_id", "det_index"]).to_csv(output / "exhaustive_yolo_crop_metadata.csv", index=False, encoding="utf-8-sig")
        head.sort_values(["taxon_name", "obs_id", "photo_id", "annotation_unit_id"]).to_csv(output / "exhaustive_continuous_head_level.csv", index=False, encoding="utf-8-sig")
        photo.sort_values("photo_id").to_csv(output / "exhaustive_continuous_photo_level.csv", index=False, encoding="utf-8-sig")
        observation.sort_values("obs_id").to_csv(output / "exhaustive_continuous_observation_level.csv", index=False, encoding="utf-8-sig")

    summaries = [json.loads(path.read_text(encoding="utf-8")) for path in paths["summary"]]
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_shards": expected,
        "n_queue_photos": int(len(queue)),
        "n_queue_observations": int(queue["obs_id"].nunique()),
        "n_queue_species": int(queue["taxon_name"].nunique()),
        "n_detected_photos": int(photo_audit["eligible_for_trait_analysis"].sum()),
        "n_no_detection_photos": int(photo_audit["screen_status"].eq("no_detection").sum()),
        "n_failed_screening_photos": int(photo_audit["screen_status"].isin(["error", "missing_image"]).sum()),
        "n_observations_with_detected_heads": int(head["obs_id"].nunique()) if not head.empty else 0,
        "n_detected_heads": int(len(head)),
        "n_colour_usable": int(head["colour_status"].astype(str).eq("usable").sum()) if not head.empty else 0,
        "n_shape_usable": int(head["shape_status"].astype(str).eq("usable").sum()) if not head.empty else 0,
        "n_orientation_usable": int(head["orientation_status"].astype(str).eq("usable").sum()) if not head.empty else 0,
        "source_sha256": {
            "shard_manifest": sha256_file(manifest_path),
            "queues": {str(path.relative_to(results_root)): sha256_file(path) for path in paths["queue"]},
            "screens": {str(path.relative_to(results_root)): sha256_file(path) for path in paths["screen"]},
            "traits": {str(path.relative_to(results_root)): sha256_file(path) for path in paths["trait"]},
        },
        "shard_summaries": summaries,
        "semantic_boundary": {
            "all_photos": "Every eligible photo is retained in detection-audit outputs.",
            "leaf_only": "Photos with no detected visible capitulum are excluded from trait tables and labelled no_visible_capitulum_detected.",
            "thinning": "No species, spatial or environmental thinning occurs before prediction or merge.",
        },
    }
    (output / "exhaustive_continuous_merge_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Join automated head-level traits to iNaturalist metadata without pseudoreplication.

The source AI measurements are head-crop level. This script preserves that
level, then derives photo and observation level categorical summaries. The
observation table is the recommended unit for geographic/environmental analysis;
head and photo tables remain as provenance and sensitivity layers.

No human labels enter this workflow. All trait states remain explicitly marked
as fully automated model-derived measurements, together with their ensemble
agreement metadata.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


AI_REQUIRED = {
    "annotation_unit_id",
    "trait_id",
    "analysis_state_ai_all",
    "analysis_state_ai_conservative",
    "ai_candidate_state",
    "ai_measurement_status",
    "measurement_basis",
}
CROP_REQUIRED = {"audit_id", "queue_id", "obs_id", "photo_id", "det_index", "yolo_conf"}
META_REQUIRED = {
    "obs_id", "photo_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess",
    "observed_on", "quality_grade", "captive", "latitude", "longitude",
    "positional_accuracy", "geoprivacy", "obscured", "coordinate_usable_for_environment",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build photo and observation analysis units from automated Cirsium trait measurements.")
    parser.add_argument("--ai-measurements-long", required=True)
    parser.add_argument("--detector-crop-metadata", required=True)
    parser.add_argument("--inat-photo-metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--metadata-chunk-size", type=int, default=200_000)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def parse_nonnegative_int(value: Any, field: str) -> int:
    raw = text(value)
    try:
        numeric = int(float(raw))
    except ValueError as error:
        raise ValueError(f"{field} must be an integer, got {raw!r}") from error
    if numeric < 0:
        raise ValueError(f"{field} must be nonnegative, got {numeric}")
    return numeric


def expected_annotation_unit_id(audit_id: Any, det_index: Any) -> str:
    return f"{text(audit_id)}_head_{parse_nonnegative_int(det_index, 'det_index') + 1:02d}"


def unique_or_empty(values: Iterable[Any], field: str) -> str:
    observed = sorted({text(value) for value in values if text(value)})
    if len(observed) > 1:
        raise ValueError(f"Inconsistent {field} values within aggregation group: {observed[:5]}")
    return observed[0] if observed else ""


def state_summary(group: pd.DataFrame, state_column: str, prefix: str) -> dict[str, Any]:
    """Summarize categorical states without arbitrarily resolving ties."""
    total = int(len(group))
    statuses = group["ai_measurement_status"].map(text)
    nonauto = total > 0 and statuses.str.startswith("not_auto_estimated").all()
    candidates = [text(value) for value in group[state_column] if text(value)]
    if nonauto:
        return {
            f"{prefix}_state": "",
            f"{prefix}_n_available": 0,
            f"{prefix}_n_distinct_states": 0,
            f"{prefix}_mode_fraction": "",
            f"{prefix}_state_counts_json": "{}",
            f"{prefix}_aggregation_status": "not_auto_estimated",
        }
    if not candidates:
        return {
            f"{prefix}_state": "",
            f"{prefix}_n_available": 0,
            f"{prefix}_n_distinct_states": 0,
            f"{prefix}_mode_fraction": "",
            f"{prefix}_state_counts_json": "{}",
            f"{prefix}_aggregation_status": "no_available_state",
        }
    counts = Counter(candidates)
    maximum = max(counts.values())
    tied = sorted(state for state, count in counts.items() if count == maximum)
    state = tied[0] if len(tied) == 1 else ""
    if len(tied) > 1:
        status = "tie_unresolved"
    elif total == 1:
        status = "single_source_unit"
    elif len(counts) == 1:
        status = "multi_source_unanimous"
    else:
        status = "multi_source_majority"
    return {
        f"{prefix}_state": state,
        f"{prefix}_n_available": int(len(candidates)),
        f"{prefix}_n_distinct_states": int(len(counts)),
        f"{prefix}_mode_fraction": float(maximum / len(candidates)),
        f"{prefix}_state_counts_json": json.dumps(dict(sorted(counts.items())), ensure_ascii=False, sort_keys=True),
        f"{prefix}_aggregation_status": status,
    }


def read_selected_metadata(path: Path, wanted_photo_ids: set[str], chunk_size: int) -> pd.DataFrame:
    if chunk_size < 1:
        raise ValueError("metadata-chunk-size must be >= 1")
    available = pd.read_csv(path, dtype=str, keep_default_na=False, nrows=0).columns.tolist()
    require_columns(pd.DataFrame(columns=available), META_REQUIRED, "iNaturalist metadata header")
    optional = ["preferred_common_name", "source_file_stem", "inat_flowers_annotation", "inat_annotation_summary"]
    usecols = sorted(META_REQUIRED.union(set(optional).intersection(available)))
    selected: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, dtype=str, keep_default_na=False, usecols=usecols, chunksize=chunk_size):
        matched = chunk.loc[chunk["photo_id"].map(text).isin(wanted_photo_ids)].copy()
        if not matched.empty:
            selected.append(matched)
    if not selected:
        raise ValueError("No detector photo IDs were found in the supplied iNaturalist metadata")
    metadata = pd.concat(selected, ignore_index=True)
    if metadata["photo_id"].duplicated().any():
        duplicate = metadata.loc[metadata["photo_id"].duplicated(keep=False), "photo_id"].head(10).tolist()
        raise ValueError(f"Merged iNaturalist metadata contains duplicate matched photo IDs: {duplicate}")
    return metadata


def build_head_table(ai: pd.DataFrame, crops: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    crops = crops.copy()
    crops["annotation_unit_id"] = [
        expected_annotation_unit_id(audit_id, det_index)
        for audit_id, det_index in zip(crops["audit_id"], crops["det_index"])
    ]
    if crops["annotation_unit_id"].duplicated().any():
        duplicate = crops.loc[crops["annotation_unit_id"].duplicated(keep=False), "annotation_unit_id"].head(10).tolist()
        raise ValueError(f"Detector crop rows do not map one-to-one to annotation units: {duplicate}")
    crop_columns = [
        "annotation_unit_id", "audit_id", "queue_id", "obs_id", "photo_id", "det_index", "yolo_conf",
    ]
    if "source_file_stem" in crops.columns:
        crop_columns.append("source_file_stem")
    crop_join = crops[crop_columns].copy()
    head = ai.merge(crop_join, on="annotation_unit_id", how="left", validate="many_to_one")
    if head["photo_id"].map(text).eq("").any():
        missing = head.loc[head["photo_id"].map(text).eq(""), "annotation_unit_id"].head(10).tolist()
        raise ValueError(f"AI measurements have no detector-crop mapping: {missing}")

    metadata_join = metadata.copy()
    metadata_join = metadata_join.rename(columns={"obs_id": "metadata_obs_id"})
    head = head.merge(metadata_join, on="photo_id", how="left", validate="many_to_one")
    if head["taxon_name"].map(text).eq("").any():
        missing = head.loc[head["taxon_name"].map(text).eq(""), "photo_id"].drop_duplicates().head(10).tolist()
        raise ValueError(f"Detector photos could not be joined to iNaturalist metadata: {missing}")
    obs_mismatch = head.loc[head["obs_id"].map(text) != head["metadata_obs_id"].map(text)]
    if not obs_mismatch.empty:
        example = obs_mismatch[["annotation_unit_id", "obs_id", "metadata_obs_id"]].head(5).to_dict("records")
        raise ValueError(f"Detector crop obs_id and iNaturalist metadata obs_id disagree: {example}")
    head = head.drop(columns=["metadata_obs_id"])
    head["coordinate_usable_for_environment"] = head["coordinate_usable_for_environment"].map(as_bool)
    head["recommended_analysis_unit"] = "observation"
    head["sampling_scope"] = "stratified_detector_audit_positive_head_subset_not_global_occurrence_sample"
    return head


def aggregate_photo(head: pd.DataFrame) -> pd.DataFrame:
    metadata_columns = [
        "obs_id", "photo_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess", "preferred_common_name",
        "observed_on", "quality_grade", "captive", "latitude", "longitude", "positional_accuracy", "geoprivacy",
        "obscured", "coordinate_usable_for_environment", "inat_flowers_annotation", "inat_annotation_summary",
    ]
    metadata_columns = [column for column in metadata_columns if column in head.columns]
    rows: list[dict[str, Any]] = []
    for (obs_id, photo_id, trait_id), group in head.groupby(["obs_id", "photo_id", "trait_id"], sort=True):
        row: dict[str, Any] = {"obs_id": text(obs_id), "photo_id": text(photo_id), "trait_id": text(trait_id)}
        for column in metadata_columns:
            if column not in row:
                row[column] = unique_or_empty(group[column], column)
        row["photo_n_detected_heads"] = int(group["annotation_unit_id"].nunique())
        row["photo_detector_confidence_mean"] = float(pd.to_numeric(group["yolo_conf"], errors="raise").mean())
        row.update(state_summary(group, "analysis_state_ai_all", "photo_ai_all"))
        row.update(state_summary(group, "analysis_state_ai_conservative", "photo_ai_conservative"))
        row["measurement_basis"] = "fully_automated_zero_shot_clip_ensemble_aggregated_from_head_crops"
        row["recommended_analysis_unit"] = "observation"
        row["sampling_scope"] = "stratified_detector_audit_positive_head_subset_not_global_occurrence_sample"
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["obs_id", "photo_id", "trait_id"]).reset_index(drop=True)


def aggregate_observation(photo: pd.DataFrame) -> pd.DataFrame:
    metadata_columns = [
        "obs_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess", "preferred_common_name",
        "observed_on", "quality_grade", "captive", "latitude", "longitude", "positional_accuracy", "geoprivacy",
        "obscured", "coordinate_usable_for_environment",
    ]
    metadata_columns = [column for column in metadata_columns if column in photo.columns]
    rows: list[dict[str, Any]] = []
    for (obs_id, trait_id), group in photo.groupby(["obs_id", "trait_id"], sort=True):
        row: dict[str, Any] = {"obs_id": text(obs_id), "trait_id": text(trait_id)}
        for column in metadata_columns:
            if column not in row:
                row[column] = unique_or_empty(group[column], column)
        row["observation_n_source_photos"] = int(group["photo_id"].nunique())
        row["observation_n_detected_heads"] = int(pd.to_numeric(group["photo_n_detected_heads"], errors="raise").sum())
        row.update(state_summary(group.rename(columns={"photo_ai_all_state": "_state"}), "_state", "observation_ai_all"))
        row.update(state_summary(group.rename(columns={"photo_ai_conservative_state": "_state"}), "_state", "observation_ai_conservative"))
        row["coordinate_usable_for_environment"] = as_bool(row.get("coordinate_usable_for_environment", ""))
        row["recommended_analysis_unit"] = "observation"
        row["measurement_basis"] = "fully_automated_zero_shot_clip_ensemble_aggregated_to_observation"
        row["sampling_scope"] = "stratified_detector_audit_positive_head_subset_not_global_occurrence_sample"
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["obs_id", "trait_id"]).reset_index(drop=True)


def build_observation_wide(observation: pd.DataFrame) -> pd.DataFrame:
    base_columns = [
        "obs_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess", "preferred_common_name", "observed_on",
        "quality_grade", "captive", "latitude", "longitude", "positional_accuracy", "geoprivacy", "obscured",
        "coordinate_usable_for_environment", "observation_n_source_photos", "observation_n_detected_heads",
        "recommended_analysis_unit", "measurement_basis", "sampling_scope",
    ]
    base_columns = [column for column in base_columns if column in observation.columns]
    base = observation[base_columns].drop_duplicates(subset=["obs_id"]).sort_values("obs_id")
    attributes = [
        "observation_ai_all_state", "observation_ai_all_aggregation_status", "observation_ai_all_mode_fraction",
        "observation_ai_conservative_state", "observation_ai_conservative_aggregation_status",
        "observation_ai_conservative_mode_fraction",
    ]
    output = base.copy()
    for attribute in attributes:
        pivot = observation.pivot(index="obs_id", columns="trait_id", values=attribute)
        pivot.columns = [f"ai__{trait_id}__{attribute}" for trait_id in pivot.columns]
        output = output.merge(pivot.reset_index(), on="obs_id", how="left", validate="one_to_one")
    return output.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    ai_path = Path(args.ai_measurements_long)
    crops_path = Path(args.detector_crop_metadata)
    metadata_path = Path(args.inat_photo_metadata)
    for path in (ai_path, crops_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    ai = pd.read_csv(ai_path, dtype=str, keep_default_na=False)
    crops = pd.read_csv(crops_path, dtype=str, keep_default_na=False)
    require_columns(ai, AI_REQUIRED, "AI measurement table")
    require_columns(crops, CROP_REQUIRED, "detector crop metadata")
    if ai.empty or ai.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("AI measurement table must have unique annotation_unit_id/trait_id rows")
    if crops.empty:
        raise ValueError("Detector crop metadata is empty")

    expected_units = set(ai["annotation_unit_id"].map(text))
    crop_units = {expected_annotation_unit_id(audit_id, det_index) for audit_id, det_index in zip(crops["audit_id"], crops["det_index"])}
    if expected_units != crop_units:
        only_ai = sorted(expected_units.difference(crop_units))[:10]
        only_crop = sorted(crop_units.difference(expected_units))[:10]
        raise ValueError(f"AI and detector crop annotation units differ; only_ai={only_ai}, only_crop={only_crop}")

    wanted_photo_ids = set(crops["photo_id"].map(text))
    metadata = read_selected_metadata(metadata_path, wanted_photo_ids, args.metadata_chunk_size)
    if set(metadata["photo_id"].map(text)) != wanted_photo_ids:
        missing = sorted(wanted_photo_ids.difference(set(metadata["photo_id"].map(text))))[:10]
        raise ValueError(f"Metadata coverage incomplete for detector photos: {missing}")

    head = build_head_table(ai, crops, metadata)
    photo = aggregate_photo(head)
    observation = aggregate_observation(photo)
    observation_wide = build_observation_wide(observation)

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    head.to_csv(output / "ai_trait_head_level_with_metadata.csv", index=False, encoding="utf-8-sig")
    photo.to_csv(output / "ai_trait_photo_level.csv", index=False, encoding="utf-8-sig")
    observation.to_csv(output / "ai_trait_observation_level_long.csv", index=False, encoding="utf-8-sig")
    observation_wide.to_csv(output / "ai_trait_observation_level_wide.csv", index=False, encoding="utf-8-sig")

    trait_counts = observation.groupby("trait_id").size().astype(int).to_dict()
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "semantic_status": "Fully automated model-derived trait measurements aggregated from detected head crops. No human annotations were used.",
        "recommended_analysis_unit": "observation",
        "sampling_scope": "Stratified 1,000-photo detector audit subset restricted to photos with one or more bootstrap-detector head proposals; not a global occurrence sample.",
        "n_head_trait_rows": int(len(head)),
        "n_detected_head_units": int(head["annotation_unit_id"].nunique()),
        "n_source_photos": int(photo["photo_id"].nunique()),
        "n_observations": int(observation["obs_id"].nunique()),
        "n_observation_trait_rows": int(len(observation)),
        "n_observations_coordinate_usable": int(observation_wide["coordinate_usable_for_environment"].map(as_bool).sum()),
        "n_traits_per_observation": trait_counts,
        "input_sha256": {
            "ai_measurements_long": sha256_file(ai_path),
            "detector_crop_metadata": sha256_file(crops_path),
            "inat_photo_metadata": sha256_file(metadata_path),
        },
        "aggregation_rules": {
            "head": "One row per detected head crop and trait; preserves raw automated ensemble result.",
            "photo": "One row per source photo and trait; categorical mode across head crops, with ties retained as unresolved rather than arbitrarily assigned.",
            "observation": "One row per iNaturalist observation and trait; categorical mode across photo-level states, with ties retained as unresolved. This is the recommended geographic/environmental analysis unit.",
            "conservative_state": "Only strict ensemble consensus with above-median trait-specific margin at head level. Empty conservative states are not imputed.",
        },
    }
    (output / "ai_trait_analysis_units_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

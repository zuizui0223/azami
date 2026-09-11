#!/usr/bin/env python3
"""Join automated head traits to iNaturalist metadata and remove nested pseudoreplication.

Input traits are detected-head crops. Outputs retain head-level provenance and
also derive photo-level and observation-level categorical summaries. The
observation table is the recommended geographic/environmental analysis unit.
No human annotations are read or inferred anywhere in this workflow.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


AI_REQUIRED = {
    "annotation_unit_id", "trait_id", "analysis_state_ai_all",
    "analysis_state_ai_conservative", "ai_candidate_state",
    "ai_measurement_status", "measurement_basis",
}
CROP_REQUIRED = {"audit_id", "queue_id", "obs_id", "photo_id", "det_index", "yolo_conf"}
META_REQUIRED = {
    "obs_id", "photo_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess",
    "observed_on", "quality_grade", "captive", "latitude", "longitude",
    "positional_accuracy", "geoprivacy", "obscured", "coordinate_usable_for_environment",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build AI trait head/photo/observation tables.")
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


def require_columns(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def det_index(value: Any) -> int:
    try:
        parsed = int(float(text(value)))
    except ValueError as error:
        raise ValueError(f"det_index must be an integer, got {value!r}") from error
    if parsed < 0:
        raise ValueError(f"det_index must be nonnegative, got {parsed}")
    return parsed


def annotation_unit_id(audit_id: Any, index: Any) -> str:
    return f"{text(audit_id)}_head_{det_index(index) + 1:02d}"


def unique_or_empty(values, field: str) -> str:
    observed = sorted({text(value) for value in values if text(value)})
    if len(observed) > 1:
        raise ValueError(f"Inconsistent {field} inside aggregation group: {observed[:5]}")
    return observed[0] if observed else ""


def categorical_summary(group: pd.DataFrame, state_col: str, prefix: str, status_col: str | None = None) -> dict[str, Any]:
    """Mode aggregation that never resolves ties arbitrarily."""
    if status_col and status_col in group.columns:
        statuses = group[status_col].map(text)
        if len(statuses) and statuses.str.startswith("not_auto_estimated").all():
            return {
                f"{prefix}_state": "", f"{prefix}_n_available": 0,
                f"{prefix}_n_distinct_states": 0, f"{prefix}_mode_fraction": "",
                f"{prefix}_state_counts_json": "{}", f"{prefix}_aggregation_status": "not_auto_estimated",
            }
    states = [text(value) for value in group[state_col] if text(value)]
    if not states:
        return {
            f"{prefix}_state": "", f"{prefix}_n_available": 0,
            f"{prefix}_n_distinct_states": 0, f"{prefix}_mode_fraction": "",
            f"{prefix}_state_counts_json": "{}", f"{prefix}_aggregation_status": "no_available_state",
        }
    counts = Counter(states)
    maximum = max(counts.values())
    tied = sorted(state for state, count in counts.items() if count == maximum)
    if len(tied) > 1:
        state, status = "", "tie_unresolved"
    else:
        state = tied[0]
        status = (
            "single_source_unit" if len(group) == 1 else
            "multi_source_unanimous" if len(counts) == 1 else
            "multi_source_majority"
        )
    return {
        f"{prefix}_state": state,
        f"{prefix}_n_available": int(len(states)),
        f"{prefix}_n_distinct_states": int(len(counts)),
        f"{prefix}_mode_fraction": float(maximum / len(states)),
        f"{prefix}_state_counts_json": json.dumps(dict(sorted(counts.items())), ensure_ascii=False, sort_keys=True),
        f"{prefix}_aggregation_status": status,
    }


def read_matched_metadata(path: Path, wanted_photo_ids: set[str], chunk_size: int) -> pd.DataFrame:
    if chunk_size < 1:
        raise ValueError("metadata-chunk-size must be >= 1")
    header = pd.read_csv(path, dtype=str, keep_default_na=False, nrows=0).columns.tolist()
    require_columns(pd.DataFrame(columns=header), META_REQUIRED, "iNaturalist metadata header")
    optional = ["preferred_common_name", "source_file_stem", "inat_flowers_annotation", "inat_annotation_summary"]
    usecols = sorted(META_REQUIRED.union(set(optional).intersection(header)))
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, dtype=str, keep_default_na=False, usecols=usecols, chunksize=chunk_size):
        selected = chunk.loc[chunk["photo_id"].map(text).isin(wanted_photo_ids)].copy()
        if not selected.empty:
            pieces.append(selected)
    if not pieces:
        raise ValueError("No detector photo IDs matched the full iNaturalist metadata")
    matched = pd.concat(pieces, ignore_index=True)
    if matched["photo_id"].duplicated().any():
        example = matched.loc[matched["photo_id"].duplicated(keep=False), "photo_id"].head(10).tolist()
        raise ValueError(f"Matched metadata photo_id must be unique: {example}")
    return matched


def aggregate_photo(head: pd.DataFrame) -> pd.DataFrame:
    metadata_fields = [
        "obs_id", "photo_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess", "preferred_common_name",
        "observed_on", "quality_grade", "captive", "latitude", "longitude", "positional_accuracy", "geoprivacy",
        "obscured", "coordinate_usable_for_environment", "inat_flowers_annotation", "inat_annotation_summary",
    ]
    metadata_fields = [field for field in metadata_fields if field in head.columns]
    rows: list[dict[str, Any]] = []
    for (obs_id, photo_id, trait_id), group in head.groupby(["obs_id", "photo_id", "trait_id"], sort=True):
        row: dict[str, Any] = {"obs_id": text(obs_id), "photo_id": text(photo_id), "trait_id": text(trait_id)}
        for field in metadata_fields:
            if field not in row:
                row[field] = unique_or_empty(group[field], field)
        row["photo_n_detected_heads"] = int(group["annotation_unit_id"].nunique())
        row["photo_detector_confidence_mean"] = float(pd.to_numeric(group["yolo_conf"], errors="raise").mean())
        row.update(categorical_summary(group, "analysis_state_ai_all", "photo_ai_all", "ai_measurement_status"))
        row.update(categorical_summary(group, "analysis_state_ai_conservative", "photo_ai_conservative", "ai_measurement_status"))
        row["measurement_basis"] = "fully_automated_zero_shot_clip_ensemble_aggregated_from_head_crops"
        row["recommended_analysis_unit"] = "observation"
        row["sampling_scope"] = "stratified_detector_audit_positive_head_subset_not_global_occurrence_sample"
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["obs_id", "photo_id", "trait_id"]).reset_index(drop=True)


def aggregate_observation(photo: pd.DataFrame) -> pd.DataFrame:
    metadata_fields = [
        "obs_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess", "preferred_common_name",
        "observed_on", "quality_grade", "captive", "latitude", "longitude", "positional_accuracy", "geoprivacy",
        "obscured", "coordinate_usable_for_environment",
    ]
    metadata_fields = [field for field in metadata_fields if field in photo.columns]
    rows: list[dict[str, Any]] = []
    for (obs_id, trait_id), group in photo.groupby(["obs_id", "trait_id"], sort=True):
        row: dict[str, Any] = {"obs_id": text(obs_id), "trait_id": text(trait_id)}
        for field in metadata_fields:
            if field not in row:
                row[field] = unique_or_empty(group[field], field)
        row["observation_n_source_photos"] = int(group["photo_id"].nunique())
        row["observation_n_detected_heads"] = int(pd.to_numeric(group["photo_n_detected_heads"], errors="raise").sum())
        all_view = group.copy()
        all_view["_status"] = all_view["photo_ai_all_aggregation_status"]
        row.update(categorical_summary(all_view, "photo_ai_all_state", "observation_ai_all", "_status"))
        strict_view = group.copy()
        strict_view["_status"] = strict_view["photo_ai_conservative_aggregation_status"]
        row.update(categorical_summary(strict_view, "photo_ai_conservative_state", "observation_ai_conservative", "_status"))
        row["coordinate_usable_for_environment"] = as_bool(row.get("coordinate_usable_for_environment", ""))
        row["measurement_basis"] = "fully_automated_zero_shot_clip_ensemble_aggregated_to_observation"
        row["recommended_analysis_unit"] = "observation"
        row["sampling_scope"] = "stratified_detector_audit_positive_head_subset_not_global_occurrence_sample"
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["obs_id", "trait_id"]).reset_index(drop=True)


def observation_wide(observation: pd.DataFrame) -> pd.DataFrame:
    base_fields = [
        "obs_id", "taxon_id", "taxon_name", "taxon_rank", "species_guess", "preferred_common_name", "observed_on",
        "quality_grade", "captive", "latitude", "longitude", "positional_accuracy", "geoprivacy", "obscured",
        "coordinate_usable_for_environment", "observation_n_source_photos", "observation_n_detected_heads",
        "recommended_analysis_unit", "measurement_basis", "sampling_scope",
    ]
    base = observation[[field for field in base_fields if field in observation.columns]].drop_duplicates("obs_id").sort_values("obs_id")
    output = base.copy()
    for attribute in [
        "observation_ai_all_state", "observation_ai_all_aggregation_status", "observation_ai_all_mode_fraction",
        "observation_ai_conservative_state", "observation_ai_conservative_aggregation_status",
        "observation_ai_conservative_mode_fraction",
    ]:
        pivot = observation.pivot(index="obs_id", columns="trait_id", values=attribute)
        pivot.columns = [f"ai__{trait_id}__{attribute}" for trait_id in pivot.columns]
        output = output.merge(pivot.reset_index(), on="obs_id", how="left", validate="one_to_one")
    return output.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    ai_path = Path(args.ai_measurements_long)
    crop_path = Path(args.detector_crop_metadata)
    metadata_path = Path(args.inat_photo_metadata)
    for path in (ai_path, crop_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    ai = pd.read_csv(ai_path, dtype=str, keep_default_na=False)
    crop = pd.read_csv(crop_path, dtype=str, keep_default_na=False)
    require_columns(ai, AI_REQUIRED, "AI measurements")
    require_columns(crop, CROP_REQUIRED, "detector crop metadata")
    if ai.empty or ai.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("AI measurements must have unique annotation_unit_id/trait_id rows")

    crop = crop.copy()
    crop["annotation_unit_id"] = [annotation_unit_id(audit, index) for audit, index in zip(crop["audit_id"], crop["det_index"])]
    if crop["annotation_unit_id"].duplicated().any():
        raise ValueError("Detector crop metadata maps multiple rows to one annotation unit")
    ai_units = set(ai["annotation_unit_id"].map(text))
    crop_units = set(crop["annotation_unit_id"].map(text))
    if ai_units != crop_units:
        raise ValueError(f"AI/detector annotation unit mismatch; only_ai={sorted(ai_units - crop_units)[:10]} only_detector={sorted(crop_units - ai_units)[:10]}")

    metadata = read_matched_metadata(metadata_path, set(crop["photo_id"].map(text)), args.metadata_chunk_size)
    if set(metadata["photo_id"].map(text)) != set(crop["photo_id"].map(text)):
        raise ValueError("Full iNaturalist metadata does not cover every detector source photo")

    crop_fields = ["annotation_unit_id", "audit_id", "queue_id", "obs_id", "photo_id", "det_index", "yolo_conf"]
    if "source_file_stem" in crop.columns:
        crop_fields.append("source_file_stem")
    head = ai.merge(crop[crop_fields], on="annotation_unit_id", how="left", validate="many_to_one")
    metadata = metadata.rename(columns={"obs_id": "metadata_obs_id"})
    head = head.merge(metadata, on="photo_id", how="left", validate="many_to_one", suffixes=("", "_metadata"))
    if head["taxon_name"].map(text).eq("").any():
        raise ValueError("Some detector photos failed to join to taxonomic/geographic metadata")
    if (head["obs_id"].map(text) != head["metadata_obs_id"].map(text)).any():
        raise ValueError("Detector and merged metadata observation IDs disagree")
    head = head.drop(columns=["metadata_obs_id"])
    head["coordinate_usable_for_environment"] = head["coordinate_usable_for_environment"].map(as_bool)
    head["recommended_analysis_unit"] = "observation"
    head["sampling_scope"] = "stratified_detector_audit_positive_head_subset_not_global_occurrence_sample"

    photo = aggregate_photo(head)
    observation = aggregate_observation(photo)
    wide = observation_wide(observation)

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    head.to_csv(output / "ai_trait_head_level_with_metadata.csv", index=False, encoding="utf-8-sig")
    photo.to_csv(output / "ai_trait_photo_level.csv", index=False, encoding="utf-8-sig")
    observation.to_csv(output / "ai_trait_observation_level_long.csv", index=False, encoding="utf-8-sig")
    wide.to_csv(output / "ai_trait_observation_level_wide.csv", index=False, encoding="utf-8-sig")
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
        "n_observations_coordinate_usable": int(wide["coordinate_usable_for_environment"].map(as_bool).sum()),
        "input_sha256": {
            "ai_measurements_long": sha256_file(ai_path),
            "detector_crop_metadata": sha256_file(crop_path),
            "inat_photo_metadata": sha256_file(metadata_path),
        },
        "aggregation_rules": {
            "head": "One detected head crop per trait.",
            "photo": "Categorical mode across detected heads; ties are retained unresolved.",
            "observation": "Categorical mode across source-photo states; ties are retained unresolved. This is the recommended analysis unit.",
            "conservative_state": "Only strict ensemble consensus with above-median trait-specific margin at head level; empty conservative states are not imputed.",
        },
    }
    (output / "ai_trait_analysis_units_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

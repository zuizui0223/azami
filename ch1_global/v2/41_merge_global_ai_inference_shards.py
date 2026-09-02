#!/usr/bin/env python3
"""Merge compact outputs from all automated global AI inference shards.

This script validates queue coverage against the original shard manifest and
requires matching ensemble configuration fingerprints before concatenating
head-level AI measurements and detector crop metadata.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


QUEUE_REQUIRED = {"queue_id", "obs_id", "photo_id", "taxon_name"}
CROP_REQUIRED = {"audit_id", "queue_id", "obs_id", "photo_id", "det_index", "yolo_conf"}
AI_REQUIRED = {"annotation_unit_id", "trait_id", "analysis_state_ai_all", "ai_measurement_status"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge compact global automated AI inference shard outputs.")
    parser.add_argument("--shard-manifest", required=True)
    parser.add_argument("--shard-results-root", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {missing}")


def find_exact(root: Path, suffix: str) -> list[Path]:
    matches = sorted(root.rglob(suffix))
    if not matches:
        raise FileNotFoundError(f"No {suffix} under {root}")
    return matches


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.shard_manifest)
    results_root = Path(args.shard_results_root)
    if not manifest_path.is_file() or not results_root.is_dir():
        raise FileNotFoundError("Shard manifest or shard-results root is missing")
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    if not {"shard_id", "filename", "n_observations"}.issubset(manifest.columns):
        raise ValueError("Shard manifest has an invalid schema")

    queue_paths = find_exact(results_root, "queue_shard.csv")
    crop_paths = find_exact(results_root, "yolo_crop_metadata.csv")
    ai_paths = find_exact(results_root, "ai_ensemble_trait_measurements_long.csv")
    provenance_paths = find_exact(results_root, "ai_ensemble_measurement_provenance.json")
    expected = len(manifest)
    if not (len(queue_paths) == len(crop_paths) == len(ai_paths) == len(provenance_paths) == expected):
        raise ValueError(
            f"Incomplete shard result collection; expected={expected}, queues={len(queue_paths)}, crops={len(crop_paths)}, ai={len(ai_paths)}, provenance={len(provenance_paths)}"
        )

    queue_parts = [pd.read_csv(path, dtype=str, keep_default_na=False) for path in queue_paths]
    crop_parts = [pd.read_csv(path, dtype=str, keep_default_na=False) for path in crop_paths]
    ai_parts = [pd.read_csv(path, dtype=str, keep_default_na=False) for path in ai_paths]
    for index, frame in enumerate(queue_parts, start=1):
        require(frame, QUEUE_REQUIRED, f"queue shard {index}")
    for index, frame in enumerate(crop_parts, start=1):
        require(frame, CROP_REQUIRED, f"crop shard {index}")
    for index, frame in enumerate(ai_parts, start=1):
        require(frame, AI_REQUIRED, f"AI shard {index}")

    queue = pd.concat(queue_parts, ignore_index=True)
    crops = pd.concat(crop_parts, ignore_index=True)
    ai = pd.concat(ai_parts, ignore_index=True)
    if queue["queue_id"].duplicated().any() or queue["obs_id"].duplicated().any():
        raise ValueError("Shard queues overlap in queue_id or observation ID")
    if crops.duplicated(["audit_id", "det_index"]).any():
        raise ValueError("Shard crop metadata overlaps in audit_id/det_index")
    if ai.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Shard AI measurements overlap in annotation_unit_id/trait_id")
    if set(crops["queue_id"].map(text)).difference(set(queue["queue_id"].map(text))):
        raise ValueError("Detector crop results reference queue IDs outside merged shard queues")

    expected_n_queue = int(pd.to_numeric(manifest["n_observations"], errors="raise").sum())
    if len(queue) != expected_n_queue:
        raise ValueError(f"Merged queue coverage differs from shard manifest: got={len(queue)} expected={expected_n_queue}")

    configurations = []
    for path in provenance_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        configurations.append({
            "model_locks": payload.get("model_locks"),
            "prompt_spec_sha256": payload.get("prompt_spec_sha256"),
            "ontology_sha256": payload.get("ontology_sha256"),
            "seed": payload.get("seed"),
            "image_variants": payload.get("image_variants"),
        })
    fingerprints = {json.dumps(configuration, sort_keys=True, ensure_ascii=False) for configuration in configurations}
    if len(fingerprints) != 1:
        raise ValueError("Shard AI ensemble configurations differ and cannot be merged as one measurement set")

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    queue.sort_values("queue_id").to_csv(output / "global_ai_merged_queue.csv", index=False, encoding="utf-8-sig")
    crops.sort_values(["audit_id", "det_index"]).to_csv(output / "global_ai_merged_yolo_crop_metadata.csv", index=False, encoding="utf-8-sig")
    ai.sort_values(["annotation_unit_id", "trait_id"]).to_csv(output / "global_ai_merged_trait_measurements_long.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_shards": expected,
        "n_queue_observations": int(len(queue)),
        "n_detected_heads": int(crops[["audit_id", "det_index"]].drop_duplicates().shape[0]),
        "n_ai_trait_rows": int(len(ai)),
        "ensemble_configuration": configurations[0],
        "source_files_sha256": {
            "shard_manifest": sha256_file(manifest_path),
            "queue_shards": {str(path.relative_to(results_root)): sha256_file(path) for path in queue_paths},
            "crop_shards": {str(path.relative_to(results_root)): sha256_file(path) for path in crop_paths},
            "ai_shards": {str(path.relative_to(results_root)): sha256_file(path) for path in ai_paths},
        },
        "semantic_status": "Fully automated model-derived measurements merged from independently processed queue shards. No human annotations were used.",
    }
    (output / "global_ai_shard_merge_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build an outcome-blind repeat-photo cohort from public photo metadata."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


URL_COLUMNS = (
    "medium_image_url",
    "large_image_url",
    "small_image_url",
    "raw_image_url",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict-observations", required=True)
    parser.add_argument("--photo-metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--chunksize", type=int, default=100_000)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def first_url(frame: pd.DataFrame) -> pd.Series:
    result = pd.Series("", index=frame.index, dtype="string")
    for column in URL_COLUMNS:
        values = frame[column].fillna("").astype(str).str.strip()
        take = result.eq("") & values.ne("")
        result.loc[take] = values.loc[take]
    return result


def build_repeat_cohort(
    strict: pd.DataFrame,
    metadata_path: Path,
    chunksize: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_strict = {"obs_id", "taxon_name", "observed_on"}
    missing = required_strict.difference(strict.columns)
    if missing:
        raise ValueError(f"Strict observation table is missing: {sorted(missing)}")
    strict = strict[["obs_id", "taxon_name", "observed_on"]].copy()
    strict["obs_id"] = strict["obs_id"].astype(str)
    if strict["obs_id"].duplicated().any():
        raise ValueError("Strict observation table must be unique by obs_id")
    strict_ids = set(strict["obs_id"])

    selected_chunks: list[pd.DataFrame] = []
    required_metadata = {
        "obs_id", "photo_id", "photo_index", "photo_license_code", *URL_COLUMNS
    }
    for chunk in pd.read_csv(
        metadata_path,
        dtype=str,
        keep_default_na=False,
        chunksize=chunksize,
        low_memory=False,
    ):
        missing = required_metadata.difference(chunk.columns)
        if missing:
            raise ValueError(f"Photo metadata are missing: {sorted(missing)}")
        chunk = chunk.loc[chunk["obs_id"].isin(strict_ids)].copy()
        if not chunk.empty:
            selected_chunks.append(chunk[[
                "obs_id", "photo_id", "photo_index", "photo_license_code", *URL_COLUMNS
            ]])
    if not selected_chunks:
        raise ValueError("No strict-cohort observations matched the photo metadata")

    photos = pd.concat(selected_chunks, ignore_index=True)
    photos = photos.drop_duplicates(["obs_id", "photo_id"], keep="first")
    photos["selected_image_url"] = first_url(photos)
    photos = photos.loc[photos["selected_image_url"].ne("")].copy()
    photos["photo_index_numeric"] = pd.to_numeric(
        photos["photo_index"], errors="coerce"
    )
    photos = photos.sort_values(
        ["obs_id", "photo_index_numeric", "photo_id"], kind="stable"
    )
    photos["photo_rank_within_observation"] = (
        photos.groupby("obs_id").cumcount() + 1
    )
    counts = photos.groupby("obs_id").size().rename("n_available_photos")
    repeat_ids = counts.loc[counts.ge(2)].index
    manifest = photos.loc[photos["obs_id"].isin(repeat_ids)].copy()
    manifest = manifest.merge(strict, on="obs_id", how="left", validate="many_to_one")
    manifest["n_available_photos"] = manifest["obs_id"].map(counts).astype(int)
    manifest = manifest[[
        "obs_id", "taxon_name", "observed_on", "photo_id", "photo_index",
        "photo_rank_within_observation", "n_available_photos",
        "photo_license_code", "selected_image_url",
    ]]

    summary = (
        manifest.groupby(["obs_id", "taxon_name", "observed_on"], as_index=False)
        .agg(n_available_photos=("photo_id", "nunique"))
        .sort_values(["n_available_photos", "obs_id"], ascending=[False, True])
    )
    summary["n_additional_photos"] = summary["n_available_photos"] - 1
    return manifest, summary


def main() -> None:
    args = parse_args()
    strict_path = Path(args.strict_observations)
    metadata_path = Path(args.photo_metadata)
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    strict = pd.read_csv(strict_path, low_memory=False)
    manifest, summary = build_repeat_cohort(strict, metadata_path, args.chunksize)
    manifest.to_csv(output / "repeat_photo_manifest.csv", index=False)
    summary.to_csv(output / "repeat_photo_observation_summary.csv", index=False)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "repeat_photo_cohort_frozen_measurements_pending",
        "selection_rule": "strict-cohort observations with at least two distinct public photo IDs having a usable image URL",
        "selection_is_outcome_blind": True,
        "n_repeat_observations": int(len(summary)),
        "n_repeat_taxa": int(summary["taxon_name"].nunique()),
        "n_photo_rows": int(len(manifest)),
        "n_additional_photos": int(summary["n_additional_photos"].sum()),
        "max_photos_per_observation": int(summary["n_available_photos"].max()),
        "interpretation_ceiling": (
            "Between-photo variance includes camera, illumination, viewpoint and possibly subject differences; "
            "it is not pure measurement error without same-individual adjudication."
        ),
        "source_sha256": {
            "strict_observations": sha256_file(strict_path),
            "photo_metadata": sha256_file(metadata_path),
        },
    }
    (output / "repeat_photo_cohort_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build a leakage-free blinded source-image audit for the frozen capitulum detector.

The audit sample is selected before detector inference. All photos and observations
used by earlier detector proposal, pseudo-label or training workflows are excluded.
The public packet omits taxonomy and coordinates; a separate key supports post-audit
stratified diagnostics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

REQUIRED = {
    "obs_id", "photo_id", "taxon_name", "taxon_rank", "source_file_stem",
}
URL_COLUMNS = ("medium_image_url", "small_image_url", "large_image_url", "raw_image_url")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--exclude-queue", action="append", default=[])
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-images", type=int, default=1000)
    parser.add_argument("--target-images-per-species", type=int, default=3)
    parser.add_argument("--max-images-per-species", type=int, default=10)
    parser.add_argument("--double-label-fraction", type=float, default=0.25)
    parser.add_argument("--spatial-block-deg", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=20260722)
    parser.add_argument("--audit-prefix", default="det_ind")
    parser.add_argument("--metadata-chunksize", type=int, default=100000)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes", "y"}


def stable_score(seed: int, *parts: Any) -> float:
    payload = "|".join([str(seed), *(text(part) for part in parts)]).encode("utf-8")
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return value / float(2**64 - 1)


def stable_scores(values: pd.Series, seed: int) -> pd.Series:
    payload = values.fillna("").astype(str) + f"|{seed}"
    hashed = pd.util.hash_pandas_object(payload, index=False).astype("uint64")
    return hashed.astype("float64") / float(np.iinfo("uint64").max)


def coordinate_block(lat: Any, lon: Any, width: float) -> str:
    try:
        lat_value, lon_value = float(lat), float(lon)
    except (TypeError, ValueError):
        return "missing"
    if not np.isfinite(lat_value) or not np.isfinite(lon_value):
        return "missing"
    return (
        f"lat{int(np.floor((lat_value + 90) / width)):03d}_"
        f"lon{int(np.floor((lon_value + 180) / width)):03d}"
    )


def latitude_band(lat: Any) -> str:
    try:
        value = float(lat)
    except (TypeError, ValueError):
        return "missing"
    if not np.isfinite(value):
        return "missing"
    if value >= 45:
        return "north_high"
    if value >= 20:
        return "north_temperate"
    if value > -20:
        return "tropical"
    if value > -45:
        return "south_temperate"
    return "south_high"


def first_url(frame: pd.DataFrame) -> pd.Series:
    values = pd.Series("", index=frame.index, dtype="object")
    for column in URL_COLUMNS:
        if column not in frame.columns:
            continue
        candidate = frame[column].fillna("").astype(str).str.strip()
        values = values.mask(values.eq("") & candidate.ne(""), candidate)
    return values


def load_exclusions(paths: Iterable[str]) -> tuple[set[str], set[str], list[dict[str, Any]]]:
    photo_ids: set[str] = set()
    obs_ids: set[str] = set()
    sources: list[dict[str, Any]] = []
    for raw_path in paths:
        path = Path(raw_path)
        table = pd.read_csv(path, dtype=str, keep_default_na=False)
        if "photo_id" not in table.columns and "obs_id" not in table.columns:
            raise ValueError(f"Exclusion table lacks photo_id/obs_id: {path}")
        before_photo, before_obs = len(photo_ids), len(obs_ids)
        if "photo_id" in table.columns:
            photo_ids.update(value for value in table["photo_id"].map(text) if value)
        if "obs_id" in table.columns:
            obs_ids.update(value for value in table["obs_id"].map(text) if value)
        sources.append(
            {
                "path": str(path),
                "rows": int(len(table)),
                "new_photo_ids": len(photo_ids) - before_photo,
                "new_obs_ids": len(obs_ids) - before_obs,
            }
        )
    return photo_ids, obs_ids, sources


def read_candidates(
    metadata: Path,
    excluded_photo_ids: set[str],
    excluded_obs_ids: set[str],
    chunksize: int,
    seed: int,
    spatial_block_deg: float,
) -> tuple[pd.DataFrame, dict[str, int]]:
    header = pd.read_csv(metadata, nrows=0)
    missing = REQUIRED.difference(header.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {sorted(missing)}")
    wanted = sorted((REQUIRED | {
        "captive", "latitude", "longitude", "quality_grade", "observed_on",
        "user_id", "user_login", "photo_license_code", "photo_attribution",
        "geoprivacy", "obscured", "coordinate_usable_for_environment",
        *URL_COLUMNS,
    }) & set(header.columns))
    raw = pd.read_csv(
        metadata, usecols=wanted, dtype=str, keep_default_na=False, low_memory=False
    )
    counts = {"metadata_rows": int(len(raw))}
    raw["taxon_name"] = raw["taxon_name"].map(text)
    raw["taxon_rank"] = raw["taxon_rank"].map(text).str.lower()
    captive = raw.get("captive", pd.Series("", index=raw.index)).map(as_bool)
    raw["selected_image_url"] = first_url(raw)
    base = (
        raw["taxon_rank"].eq("species")
        & raw["taxon_name"].ne("")
        & ~captive
        & raw["selected_image_url"].ne("")
    )
    counts["species_rank_non_captive_photo_rows"] = int(base.sum())
    excluded_photo = base & raw["photo_id"].map(text).isin(excluded_photo_ids)
    counts["excluded_photo_rows"] = int(excluded_photo.sum())
    base &= ~excluded_photo
    excluded_obs = base & raw["obs_id"].map(text).isin(excluded_obs_ids)
    counts["excluded_observation_rows"] = int(excluded_obs.sum())
    base &= ~excluded_obs
    candidates = raw.loc[base].copy()
    candidates["_score"] = stable_scores(candidates["photo_id"], seed).to_numpy()
    candidates = (
        candidates.sort_values(["obs_id", "_score"], kind="stable")
        .groupby("obs_id", sort=False)
        .head(1)
        .copy()
    )
    candidates["spatial_block"] = [
        coordinate_block(lat, lon, spatial_block_deg)
        for lat, lon in zip(
            candidates.get("latitude", pd.Series("", index=candidates.index)),
            candidates.get("longitude", pd.Series("", index=candidates.index)),
        )
    ]
    candidates["latitude_band"] = [
        latitude_band(value)
        for value in candidates.get("latitude", pd.Series("", index=candidates.index))
    ]
    counts["candidate_observations_after_exclusion"] = int(len(candidates))
    counts["candidate_species_after_exclusion"] = int(candidates["taxon_name"].nunique())
    return candidates, counts


def order_species(part: pd.DataFrame, seed: int) -> list[int]:
    work = part.copy()
    payload = (
        work["taxon_name"].astype(str) + "|" + work["spatial_block"].astype(str)
        + "|" + work["photo_id"].astype(str)
    )
    work["_within_score"] = stable_scores(payload, seed).to_numpy()
    first = (
        work.sort_values(["spatial_block", "_within_score"], kind="stable")
        .groupby("spatial_block", sort=False)
        .head(1)
        .sort_values("_within_score", kind="stable")
    )
    remaining = work.drop(index=first.index).sort_values("_within_score", kind="stable")
    return [int(index) for index in pd.concat([first, remaining]).index]


def choose_sample(
    candidates: pd.DataFrame,
    n_images: int,
    target_per_species: int,
    max_per_species: int,
    seed: int,
) -> pd.DataFrame:
    species = sorted(candidates["taxon_name"].unique())
    if n_images < len(species):
        raise ValueError(
            f"n-images={n_images} is smaller than the {len(species)} remaining species; "
            "increase n-images or define an explicit taxonomic subset"
        )
    queues = {
        name: order_species(candidates.loc[candidates["taxon_name"].eq(name)], seed)
        for name in species
    }
    selected: list[int] = []
    counts = {name: 0 for name in species}
    for _ in range(target_per_species):
        for name in species:
            if queues[name] and counts[name] < max_per_species:
                selected.append(queues[name].pop(0))
                counts[name] += 1
    capacity = sum(min(max_per_species, int((candidates["taxon_name"] == name).sum())) for name in species)
    target = min(n_images, capacity)
    while len(selected) < target:
        progressed = False
        for name in species:
            if len(selected) >= target:
                break
            if queues[name] and counts[name] < max_per_species:
                selected.append(queues[name].pop(0))
                counts[name] += 1
                progressed = True
        if not progressed:
            break
    return candidates.loc[selected].copy().reset_index(drop=True)


def main() -> None:
    args = parse_args()
    if args.n_images < 1 or args.target_images_per_species < 1:
        raise ValueError("n-images and target-images-per-species must be positive")
    if args.max_images_per_species < args.target_images_per_species:
        raise ValueError("max-images-per-species must be >= target-images-per-species")
    if not 0 <= args.double_label_fraction <= 1:
        raise ValueError("double-label-fraction must be in [0,1]")
    if args.spatial_block_deg <= 0 or args.metadata_chunksize < 1:
        raise ValueError("spatial-block-deg and metadata-chunksize must be positive")
    if not args.audit_prefix.replace("_", "").replace("-", "").isalnum():
        raise ValueError("audit-prefix may contain only letters, digits, hyphens and underscores")

    excluded_photos, excluded_obs, exclusion_sources = load_exclusions(args.exclude_queue)
    candidates, counts = read_candidates(
        Path(args.metadata), excluded_photos, excluded_obs, args.metadata_chunksize,
        args.seed, args.spatial_block_deg,
    )
    selected = choose_sample(
        candidates, args.n_images, args.target_images_per_species,
        args.max_images_per_species, args.seed,
    )
    selected.insert(0, "audit_id", [f"{args.audit_prefix}_{index:05d}" for index in range(1, len(selected) + 1)])
    selected.insert(1, "queue_id", selected["audit_id"])
    selected["screen_download_filename"] = [
        f"{audit_id}_photo_{text(photo_id)}.jpg"
        for audit_id, photo_id in zip(selected["audit_id"], selected["photo_id"])
    ]
    selected["audit_unit"] = "source_image"
    selected["annotation_status"] = "not_started"
    selected["double_label"] = False
    n_double = int(round(len(selected) * args.double_label_fraction))
    if n_double:
        ordered = sorted(
            selected.index,
            key=lambda index: stable_score(args.seed + 1, selected.at[index, "audit_id"]),
        )
        selected.loc[ordered[:n_double], "double_label"] = True

    overlap_photo = set(selected["photo_id"].map(text)) & excluded_photos
    overlap_obs = set(selected["obs_id"].map(text)) & excluded_obs
    if overlap_photo or overlap_obs:
        raise AssertionError("Independent audit selection overlaps prior proposal/training source records")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    queue_columns = [
        "queue_id", "audit_id", "obs_id", "photo_id", "source_file_stem",
        "screen_download_filename", "small_image_url", "medium_image_url",
        "large_image_url", "raw_image_url", "selected_image_url",
    ]
    for column in queue_columns:
        if column not in selected:
            selected[column] = ""
    selected[queue_columns].to_csv(out / "detector_independent_audit_queue.csv", index=False)

    blinded_columns = [
        "audit_id", "queue_id", "screen_download_filename", "audit_unit",
        "annotation_status", "double_label",
    ]
    selected[blinded_columns].to_csv(out / "detector_independent_audit_blinded_manifest.csv", index=False)

    key_columns = [
        "audit_id", "queue_id", "obs_id", "photo_id", "taxon_name", "spatial_block",
        "latitude_band", "latitude", "longitude", "quality_grade", "observed_on",
        "user_id", "user_login", "photo_license_code", "photo_attribution", "geoprivacy",
        "obscured", "coordinate_usable_for_environment", "selected_image_url",
    ]
    for column in key_columns:
        if column not in selected:
            selected[column] = ""
    selected[key_columns].to_csv(out / "detector_independent_audit_private_key.csv", index=False)

    assignments: list[dict[str, Any]] = []
    for _, row in selected.iterrows():
        assignments.append({
            "audit_id": row["audit_id"], "annotator_id": "annotator_1",
            "double_label": bool(row["double_label"]), "status": "not_started",
        })
        if bool(row["double_label"]):
            assignments.append({
                "audit_id": row["audit_id"], "annotator_id": "annotator_2",
                "double_label": True, "status": "not_started",
            })
    pd.DataFrame(assignments).to_csv(out / "detector_independent_audit_assignments.csv", index=False)
    pd.DataFrame(columns=[
        "audit_id", "annotator_id", "annotation_status", "assessability",
        "image_quality", "image_width", "image_height", "gt_index", "gt_label",
        "life_stage", "occlusion", "edge_truncated", "gt_x1", "gt_y1", "gt_x2",
        "gt_y2", "notes",
    ]).to_csv(out / "detector_independent_audit_annotation_template.csv", index=False)

    availability = candidates.groupby("taxon_name").size().rename("n_candidate_images").reset_index()
    selected_qc = selected.groupby("taxon_name").agg(
        n_selected=("audit_id", "size"), n_spatial_blocks=("spatial_block", "nunique"),
        n_latitude_bands=("latitude_band", "nunique"),
    ).reset_index()
    availability.merge(selected_qc, on="taxon_name", how="left").fillna(0).to_csv(
        out / "detector_independent_audit_species_qc.csv", index=False
    )

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "metadata": str(Path(args.metadata)),
        "selection_stage": "before frozen detector inference",
        "independence_definition": "exclude every photo_id and obs_id in all prior audit/proposal/training source queues",
        "exclusion_sources": exclusion_sources,
        "n_unique_excluded_photo_ids": len(excluded_photos),
        "n_unique_excluded_obs_ids": len(excluded_obs),
        **counts,
        "n_selected_images": int(len(selected)),
        "n_selected_species": int(selected["taxon_name"].nunique()),
        "n_double_label_images": int(selected["double_label"].sum()),
        "target_images_per_species": args.target_images_per_species,
        "max_images_per_species": args.max_images_per_species,
        "seed": args.seed,
        "photo_overlap_with_exclusions": len(overlap_photo),
        "observation_overlap_with_exclusions": len(overlap_obs),
        "annotation_blinding": "taxon, coordinates and detector predictions are omitted from the public annotation packet",
        "audit_status": "selected and blinded; independent manual annotations still required",
    }
    (out / "detector_independent_audit_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

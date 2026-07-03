#!/usr/bin/env python3
"""Create a blinded, taxonomically and spatially spread detector-audit sample.

This sample is selected *before* YOLO screening. Human annotators draw every
visible capitulum box (or record none/unassessable) on source images. The result
is a detector audit, not a trait-label dataset.
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

REQUIRED = {"obs_id", "photo_id", "taxon_name", "taxon_rank", "source_file_stem"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a blinded Cirsium YOLO detector-audit manifest.")
    parser.add_argument("--metadata", required=True, help="Merged photo_metadata.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-images", type=int, default=800)
    parser.add_argument("--min-images-per-species", type=int, default=3)
    parser.add_argument("--max-images-per-species", type=int, default=12)
    parser.add_argument("--double-label-fraction", type=float, default=0.25)
    parser.add_argument("--spatial-block-deg", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=20260703)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def coordinate_block(lat: Any, lon: Any, width: float) -> str:
    try:
        lat_value, lon_value = float(lat), float(lon)
    except (TypeError, ValueError):
        return "missing"
    if not np.isfinite(lat_value) or not np.isfinite(lon_value):
        return "missing"
    return f"lat{int(np.floor((lat_value + 90) / width)):03d}_lon{int(np.floor((lon_value + 180) / width)):03d}"


def usable_url(row: pd.Series) -> str:
    for field in ("medium_image_url", "small_image_url", "large_image_url", "raw_image_url"):
        value = text(row.get(field, ""))
        if value:
            return value
    return ""


def choose_balanced(candidates: pd.DataFrame, n_images: int, max_per_species: int, rng: np.random.Generator) -> pd.DataFrame:
    """Spread across species × spatial block first, then round-robin species."""
    if candidates.empty or n_images <= 0:
        return candidates.iloc[0:0].copy()
    work = candidates.copy()
    work["_random"] = rng.random(len(work))
    work = work.sort_values(["taxon_name", "spatial_block", "_random"])
    first = work.groupby(["taxon_name", "spatial_block"], sort=False).head(1).copy()
    first_by_species = {name: list(part.index) for name, part in first.groupby("taxon_name", sort=True)}
    selected: list[int] = []
    counts: dict[str, int] = {}
    while len(selected) < n_images and any(first_by_species.values()):
        for species in sorted(first_by_species):
            if len(selected) >= n_images:
                break
            if first_by_species[species] and counts.get(species, 0) < max_per_species:
                index = first_by_species[species].pop(0)
                selected.append(index)
                counts[species] = counts.get(species, 0) + 1
    remaining = work.drop(index=selected, errors="ignore")
    while len(selected) < n_images and not remaining.empty:
        remaining = remaining.assign(_selected_n=remaining["taxon_name"].map(counts).fillna(0))
        eligible = remaining.loc[remaining["_selected_n"] < max_per_species]
        if eligible.empty:
            break
        index = eligible.sort_values(["_selected_n", "_random"]).index[0]
        species = text(remaining.loc[index, "taxon_name"])
        selected.append(index)
        counts[species] = counts.get(species, 0) + 1
        remaining = remaining.drop(index=index)
    return work.loc[selected].drop(columns=["_random"])


def main() -> None:
    args = parse_args()
    if args.n_images < 1 or args.min_images_per_species < 1 or args.max_images_per_species < args.min_images_per_species:
        raise ValueError("Require n_images >= 1 and 1 <= min-images-per-species <= max-images-per-species")
    if not 0 <= args.double_label_fraction <= 1 or args.spatial_block_deg <= 0:
        raise ValueError("double-label-fraction must be in [0,1] and spatial-block-deg > 0")
    raw = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    missing = REQUIRED.difference(raw.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {sorted(missing)}")
    raw["taxon_name"] = raw["taxon_name"].map(text)
    raw["taxon_rank"] = raw["taxon_rank"].map(text).str.lower()
    raw["captive"] = raw.get("captive", pd.Series("", index=raw.index)).map(as_bool)
    raw["selected_image_url"] = raw.apply(usable_url, axis=1)
    candidates = raw.loc[
        raw["taxon_rank"].eq("species")
        & raw["taxon_name"].ne("")
        & ~raw["captive"]
        & raw["selected_image_url"].ne("")
    ].copy()
    if candidates.empty:
        raise ValueError("No species-rank, non-captive, photo-bearing candidates")
    candidates["spatial_block"] = [
        coordinate_block(row.get("latitude", ""), row.get("longitude", ""), args.spatial_block_deg)
        for _, row in candidates.iterrows()
    ]
    rng = np.random.default_rng(args.seed)
    candidates["_random"] = rng.random(len(candidates))
    # One photograph per observation prevents serial-photo inflation in detector evaluation.
    candidates = candidates.sort_values(["obs_id", "_random"]).groupby("obs_id", sort=False).head(1).drop(columns="_random")
    availability = candidates.groupby("taxon_name").size().rename("n_candidate_images").reset_index()
    eligible_species = set(availability.loc[availability["n_candidate_images"] >= args.min_images_per_species, "taxon_name"])
    candidates = candidates.loc[candidates["taxon_name"].isin(eligible_species)].copy()
    if candidates.empty:
        raise ValueError("No species meet --min-images-per-species")
    selected = choose_balanced(candidates, min(args.n_images, len(candidates)), args.max_images_per_species, rng).reset_index(drop=True)
    selected.insert(0, "audit_id", [f"det_audit_{i:05d}" for i in range(1, len(selected) + 1)])
    selected.insert(1, "queue_id", selected["audit_id"])
    selected["screen_download_filename"] = [f"detector_audit_{i:05d}_photo_{text(photo_id)}.jpg" for i, photo_id in enumerate(selected["photo_id"], start=1)]
    selected["audit_unit"] = "source_image"
    selected["annotation_status"] = "not_started"
    selected["needs_all_visible_capitulum_boxes"] = True
    selected["double_label"] = False
    n_double = int(round(len(selected) * args.double_label_fraction))
    if n_double:
        double_indices = rng.choice(selected.index.to_numpy(), size=n_double, replace=False)
        selected.loc[double_indices, "double_label"] = True

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queue_columns = [
        "queue_id", "audit_id", "obs_id", "photo_id", "source_file_stem", "screen_download_filename",
        "small_image_url", "medium_image_url", "large_image_url", "raw_image_url", "selected_image_url",
    ]
    for column in queue_columns:
        if column not in selected.columns:
            selected[column] = ""
    selected[queue_columns].to_csv(out_dir / "detector_audit_queue.csv", index=False, encoding="utf-8-sig")
    blinded = selected[["audit_id", "queue_id", "photo_id", "source_file_stem", "screen_download_filename", "audit_unit", "annotation_status", "double_label"]].copy()
    blinded.to_csv(out_dir / "detector_audit_blinded_manifest.csv", index=False, encoding="utf-8-sig")
    key_columns = [
        "audit_id", "queue_id", "obs_id", "photo_id", "taxon_name", "spatial_block", "latitude", "longitude",
        "geoprivacy", "obscured", "coordinate_usable_for_environment", "selected_image_url",
    ]
    for column in key_columns:
        if column not in selected.columns:
            selected[column] = ""
    selected[key_columns].to_csv(out_dir / "detector_audit_key.csv", index=False, encoding="utf-8-sig")

    assignments: list[dict[str, object]] = []
    for _, row in blinded.iterrows():
        assignments.append({"audit_id": row["audit_id"], "annotator_id": "annotator_1", "double_label": row["double_label"], "status": "not_started"})
        if bool(row["double_label"]):
            assignments.append({"audit_id": row["audit_id"], "annotator_id": "annotator_2", "double_label": True, "status": "not_started"})
    pd.DataFrame(assignments).to_csv(out_dir / "detector_audit_assignments.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(columns=[
        "audit_id", "annotator_id", "annotation_status", "assessability", "gt_index", "gt_label",
        "gt_x1", "gt_y1", "gt_x2", "gt_y2", "notes",
    ]).to_csv(out_dir / "detector_audit_annotation_template.csv", index=False, encoding="utf-8-sig")
    selected_qc = selected.groupby("taxon_name").agg(n_selected=("audit_id", "size"), n_spatial_blocks=("spatial_block", "nunique")).reset_index()
    availability.merge(selected_qc, on="taxon_name", how="left").fillna({"n_selected": 0, "n_spatial_blocks": 0}).sort_values("n_selected", ascending=False).to_csv(out_dir / "detector_audit_species_qc.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "metadata": str(Path(args.metadata).resolve()),
        "n_metadata_rows": int(len(raw)),
        "n_eligible_species": int(len(eligible_species)),
        "n_audit_images": int(len(selected)),
        "n_double_label_images": int(selected["double_label"].sum()),
        "seed": args.seed,
        "note": "The blinded manifest omits taxon and coordinate fields. Annotators must mark all visible capitula, none, or unassessable; detector predictions are not shown.",
    }
    (out_dir / "detector_audit_manifest_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

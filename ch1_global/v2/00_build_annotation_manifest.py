#!/usr/bin/env python3
"""Create the image-level annotation manifest for Chapter 1 v2.

Input: v1 YOLO crop metadata, one detected capitulum per row.
Output: an observation-deduplicated, species-balanced annotation pool plus
explicit observation/species/spatial grouping columns for grouped CV and LOSO.

Example
-------
python ch1_global/v2/00_build_annotation_manifest.py \
  --input <flowering_annotation_yolo_crop_metadata.csv> \
  --out-dir <annotation_dir> \
  --n-calibration 720 --n-main 3000 --n-segmentation 1000
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


CANDIDATES = {
    "observation_group": ["obs_id", "observation_id", "inat_obs_id"],
    "photo_group": ["photo_id", "image_id", "inat_photo_id"],
    "species": ["taxon_name", "scientific_name", "species", "taxon"],
    "latitude": ["latitude", "lat", "decimalLatitude"],
    "longitude": ["longitude", "lon", "decimalLongitude"],
    "detector_confidence": ["yolo_conf", "detector_confidence", "confidence"],
    "source_image": ["source_image", "flower_local_path", "image_path", "local_path"],
    "crop_path": ["crop_path", "head_crop_path"],
    "raw_image_url": ["raw_image_url", "keep_image_url", "image_url"],
    "positional_accuracy": ["positional_accuracy", "coordinate_uncertainty_m"],
    "geoprivacy": ["geoprivacy"],
    "obscured": ["obscured"],
    "photo_index": ["photo_index"],
    "det_index": ["det_index", "detection_index"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build a Ch.1 v2 annotation manifest.")
    p.add_argument("--input", required=True, help="v1 YOLO crop metadata CSV")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--seed", type=int, default=20260703)
    p.add_argument("--min-yolo-conf", type=float, default=0.25)
    p.add_argument("--min-heads-per-species", type=int, default=30)
    p.add_argument("--max-heads-per-species", type=int, default=80)
    p.add_argument("--n-calibration", type=int, default=720)
    p.add_argument("--n-main", type=int, default=3000)
    p.add_argument("--n-segmentation", type=int, default=1000)
    p.add_argument("--spatial-block-deg", type=float, default=10.0)
    p.add_argument("--n-cv-folds", type=int, default=5)
    p.add_argument("--allow-obscured", action="store_true")
    return p.parse_args()


def find_column(df: pd.DataFrame, logical_name: str, required: bool = False) -> str | None:
    for name in CANDIDATES[logical_name]:
        if name in df.columns:
            return name
    if required:
        raise ValueError(
            f"Missing '{logical_name}'. Expected one of {CANDIDATES[logical_name]}; "
            f"available columns: {list(df.columns)}"
        )
    return None


def text(x: object) -> str:
    return "" if pd.isna(x) else str(x).strip()


def fold(value: object, n_folds: int, seed: int) -> int:
    digest = hashlib.sha256(f"{seed}|{text(value)}".encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % n_folds + 1


def make_spatial_block(lat: object, lon: object, width: float) -> str:
    try:
        lat_f, lon_f = float(lat), float(lon)
    except (TypeError, ValueError):
        return "missing"
    if not np.isfinite(lat_f) or not np.isfinite(lon_f):
        return "missing"
    return f"lat{int(np.floor((lat_f + 90) / width)):03d}_lon{int(np.floor((lon_f + 180) / width)):03d}"


def choose_spatially_spread(df: pd.DataFrame, n: int, rng: np.random.Generator) -> pd.DataFrame:
    """One record per species × spatial block first, then fill with low-count species."""
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()
    work = df.copy()
    work["_random"] = rng.random(len(work))
    seed_rows = (
        work.sort_values(["species", "spatial_block", "_random"])
        .groupby(["species", "spatial_block"], sort=False)
        .head(1)
    )
    selected = []
    by_species = {
        species: list(part.index)
        for species, part in seed_rows.groupby("species", sort=True)
    }
    while len(selected) < n and any(by_species.values()):
        for species in sorted(by_species):
            if len(selected) >= n:
                break
            if by_species[species]:
                selected.append(by_species[species].pop(0))
    remaining = work.drop(index=selected, errors="ignore")
    counts = work.loc[selected, "species"].value_counts().to_dict() if selected else {}
    while len(selected) < n and not remaining.empty:
        remaining = remaining.assign(_count=remaining["species"].map(counts).fillna(0))
        pick = remaining.sort_values(["_count", "_random"]).index[0]
        species = remaining.loc[pick, "species"]
        selected.append(pick)
        counts[species] = int(counts.get(species, 0)) + 1
        remaining = remaining.drop(index=pick)
    return work.loc[selected].drop(columns="_random")


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.input, low_memory=False)
    if raw.empty:
        raise ValueError("Input CSV is empty")

    resolved = {
        key: find_column(raw, key, required=key in {"observation_group", "species", "detector_confidence"})
        for key in CANDIDATES
    }
    df = pd.DataFrame(index=raw.index)
    for key, source in resolved.items():
        df[key] = raw[source] if source else pd.NA

    df["species"] = df["species"].map(text)
    df["detector_confidence"] = pd.to_numeric(df["detector_confidence"], errors="coerce")
    df = df.loc[df["species"].ne("") & df["detector_confidence"].ge(args.min_yolo_conf)].copy()
    if df.empty:
        raise ValueError("No candidate heads remain after species/confidence filtering")

    df["observation_group"] = df["observation_group"].map(text)
    df["photo_group"] = df["photo_group"].map(text)
    missing_obs = df["observation_group"].eq("")
    df.loc[missing_obs, "observation_group"] = "photo_" + df.loc[missing_obs, "photo_group"].replace("", "missing").astype(str)
    df.loc[df["photo_group"].eq(""), "photo_group"] = df.loc[df["photo_group"].eq(""), "observation_group"]

    df["spatial_block"] = [
        make_spatial_block(lat, lon, args.spatial_block_deg)
        for lat, lon in zip(df["latitude"], df["longitude"])
    ]
    geoprivacy = df["geoprivacy"].map(text).str.lower().isin({"obscured", "private"})
    obscured = df["obscured"].astype(str).str.lower().isin({"true", "1", "yes"})
    df["coordinate_usable_for_environment"] = ~(geoprivacy | obscured | df["spatial_block"].eq("missing"))
    if not args.allow_obscured:
        df = df.loc[df["coordinate_usable_for_environment"]].copy()

    # One selected head per observation for human trait labelling: no serial copies of the same individual.
    df["_random"] = rng.random(len(df))
    df = (
        df.sort_values(["observation_group", "detector_confidence", "_random"], ascending=[True, False, True])
        .groupby("observation_group", sort=False)
        .head(1)
        .drop(columns="_random")
        .copy()
    )
    species_counts = df["species"].value_counts()
    df = df.loc[df["species"].map(species_counts).ge(args.min_heads_per_species)].copy()
    if df.empty:
        raise ValueError("No species meet --min-heads-per-species after observation deduplication")

    # Cap common species while retaining their geographic breadth.
    capped = []
    for _, part in df.groupby("species", sort=True):
        capped.append(choose_spatially_spread(part, min(len(part), args.max_heads_per_species), rng))
    pool = pd.concat(capped, ignore_index=True)
    main = choose_spatially_spread(pool, min(args.n_main, len(pool)), rng).copy()
    main["species_cv_fold"] = main["species"].map(lambda x: fold(x, args.n_cv_folds, args.seed))
    main["spatial_cv_fold"] = main["spatial_block"].map(lambda x: fold(x, args.n_cv_folds, args.seed + 1))

    calibration = choose_spatially_spread(main, min(args.n_calibration, len(main)), rng)
    calibration_signature = set((calibration["observation_group"].astype(str) + "|" + calibration["photo_group"].astype(str)).tolist())
    main_signature = main["observation_group"].astype(str) + "|" + main["photo_group"].astype(str)
    main["annotation_batch"] = np.where(main_signature.isin(calibration_signature), "calibration_double_label", "main_single_label")

    segmentation = choose_spatially_spread(main, min(args.n_segmentation, len(main)), rng)
    segmentation_signature = set((segmentation["observation_group"].astype(str) + "|" + segmentation["photo_group"].astype(str)).tolist())
    main["needs_segmentation"] = main_signature.isin(segmentation_signature)
    main["needs_detector_box"] = True
    main["needs_keypoints"] = True
    main["needs_trait_labels"] = True
    main["camera_tilt_flag"] = "unlabelled"
    main["view_angle_class"] = "unlabelled"
    main["flowering_qc"] = "unlabelled"
    main["annotation_status"] = "not_started"
    main = main.reset_index(drop=True)
    main.insert(0, "annotation_unit_id", [f"ch1_{i:06d}" for i in range(1, len(main) + 1)])

    desired = [
        "annotation_unit_id", "annotation_batch", "annotation_status",
        "observation_group", "photo_group", "species", "species_cv_fold", "spatial_block", "spatial_cv_fold",
        "latitude", "longitude", "coordinate_usable_for_environment", "positional_accuracy", "geoprivacy", "obscured",
        "source_image", "crop_path", "raw_image_url", "detector_confidence", "photo_index", "det_index",
        "needs_detector_box", "needs_keypoints", "needs_segmentation", "needs_trait_labels",
        "camera_tilt_flag", "view_angle_class", "flowering_qc",
    ]
    for name in desired:
        if name not in main.columns:
            main[name] = pd.NA
    manifest = main[desired].sort_values(["annotation_batch", "species", "spatial_block", "annotation_unit_id"])

    assignments = []
    for _, row in manifest.iterrows():
        annotators = ("annotator_1", "annotator_2") if row["annotation_batch"] == "calibration_double_label" else ("annotator_1",)
        for annotator in annotators:
            assignments.append({
                "annotation_unit_id": row["annotation_unit_id"],
                "annotator_id": annotator,
                "annotation_batch": row["annotation_batch"],
                "image_to_annotate": row["source_image"],
                "crop_to_annotate": row["crop_path"],
                "blind_taxonomy": True,
                "blind_coordinates": True,
                "status": "not_started",
                "adjudication_required": row["annotation_batch"] == "calibration_double_label",
                "notes": "",
            })
    assignments = pd.DataFrame(assignments)

    species_qc = (
        manifest.groupby(["species", "annotation_batch"], dropna=False)
        .agg(n_units=("annotation_unit_id", "size"), n_spatial_blocks=("spatial_block", "nunique"))
        .reset_index()
        .sort_values(["annotation_batch", "n_units", "species"], ascending=[True, False, True])
    )
    run_info = pd.DataFrame([{
        "input": str(Path(args.input).resolve()),
        "rows_input": len(raw),
        "manifest_units": len(manifest),
        "calibration_units": int((manifest["annotation_batch"] == "calibration_double_label").sum()),
        "segmentation_units": int(manifest["needs_segmentation"].sum()),
        "n_species": manifest["species"].nunique(),
        "seed": args.seed,
        "min_yolo_conf": args.min_yolo_conf,
        "min_heads_per_species": args.min_heads_per_species,
        "max_heads_per_species": args.max_heads_per_species,
        "spatial_block_deg": args.spatial_block_deg,
    }])

    manifest.to_csv(out_dir / "annotation_manifest.csv", index=False, encoding="utf-8-sig")
    assignments.to_csv(out_dir / "annotation_assignments.csv", index=False, encoding="utf-8-sig")
    species_qc.to_csv(out_dir / "annotation_species_qc.csv", index=False, encoding="utf-8-sig")
    run_info.to_csv(out_dir / "annotation_manifest_run_info.csv", index=False, encoding="utf-8-sig")
    print("[OK] Wrote Ch.1 annotation manifest to", out_dir)
    print(run_info.to_string(index=False))
    print("Use observation_group, species_cv_fold, and spatial_cv_fold for validation; never headline a random image split.")


if __name__ == "__main__":
    main()

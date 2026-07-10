#!/usr/bin/env python3
"""Build a simple app-ready representative human holdout.

For each requested trait, select one proportional sample across taxon x 10-degree
spatial cells. The primary annotator receives the full sample; the secondary
annotator receives a blinded subset for inter-annotator agreement. No additional
species top-up or extra annotation stratum is created.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from PIL import Image

DEFAULT_TRAITS = ("capitulum_orientation", "corolla_colour_class", "capitulum_shape")
UNASSESSABLE = {"", "unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--head-traits", required=True)
    parser.add_argument("--yolo-crops", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--ontology", default=None)
    parser.add_argument("--traits", nargs="+", default=list(DEFAULT_TRAITS))
    parser.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative")
    parser.add_argument("--n-tasks-per-trait", type=int, default=120)
    parser.add_argument("--double-label-fraction", type=float, default=0.25)
    parser.add_argument("--context-pad-ratio", type=float, default=1.5)
    parser.add_argument("--seed", default="20260710")
    parser.add_argument("--sleep-sec", type=float, default=0.15)
    parser.add_argument("--timeout-sec", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=4)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def stable_rank(value: str, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def block_10deg(latitude: Any, longitude: Any) -> str:
    try:
        latitude = float(latitude)
        longitude = float(longitude)
    except (TypeError, ValueError):
        return ""
    if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
        return ""
    lat_index = min(17, max(0, math.floor((latitude + 90.0) / 10.0)))
    lon_index = min(35, max(0, math.floor((longitude + 180.0) / 10.0)))
    return f"lat{lat_index:02d}_lon{lon_index:02d}"


def representative_sample(part: pd.DataFrame, n_target: int, seed: str) -> pd.DataFrame:
    """Deterministically allocate the sample in proportion to each taxon x cell."""
    if part.empty or n_target <= 0:
        return part.head(0)
    part = part.copy()
    limit = min(n_target, len(part))
    part["_rank"] = part["annotation_unit_id"].map(lambda value: stable_rank(text(value), seed))
    groups = ["taxon_name", "spatial_block_10deg"]
    counts = part.groupby(groups, sort=True).size().rename("cell_size").reset_index()
    counts["exact"] = counts["cell_size"] * limit / len(part)
    counts["allocation"] = counts["exact"].map(math.floor).astype(int)
    counts["remainder"] = counts["exact"] - counts["allocation"]
    counts["_tie"] = [
        stable_rank(f"{taxon}|{block}", f"{seed}:allocation")
        for taxon, block in zip(counts["taxon_name"], counts["spatial_block_10deg"])
    ]
    missing = limit - int(counts["allocation"].sum())
    if missing:
        indices = counts.sort_values(["remainder", "_tie"], ascending=[False, True]).index[:missing]
        counts.loc[indices, "allocation"] += 1
    if int(counts["allocation"].sum()) != limit:
        raise ValueError("Representative allocation did not reach the requested sample size")

    allocation = counts.set_index(groups)["allocation"]
    selected: list[pd.DataFrame] = []
    for keys, group in part.groupby(groups, sort=True):
        n = int(allocation.loc[keys])
        if n:
            selected.append(group.sort_values("_rank").head(n))
    result = pd.concat(selected, ignore_index=True).drop(columns="_rank", errors="ignore")
    if len(result) != limit:
        raise ValueError("Representative sample has an unexpected size")
    return result


def require(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = columns.difference(frame.columns)
    if missing:
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def annotation_unit_from_yolo(row: pd.Series) -> str:
    queue = text(row.get("audit_id")) or text(row.get("queue_id"))
    return f"{queue}_head_{int(float(row['det_index'])) + 1:02d}"


def download_image(
    session: requests.Session,
    url: str,
    target: Path,
    timeout: float,
    retries: int,
) -> None:
    last_error = ""
    for attempt in range(1, retries + 1):
        temporary: Path | None = None
        try:
            response = session.get(url, timeout=timeout, stream=True)
            response.raise_for_status()
            target.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile("wb", delete=False, dir=target.parent, suffix=".part") as handle:
                temporary = Path(handle.name)
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
            with Image.open(temporary) as image:
                image.verify()
            temporary.replace(target)
            return
        except Exception as error:
            last_error = f"{type(error).__name__}: {error}"
            if temporary is not None:
                temporary.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(8.0, 2.0 ** (attempt - 1)))
    raise RuntimeError(last_error or "Image download failed")


def crop_box(image: Image.Image, row: pd.Series, pad_ratio: float) -> tuple[Image.Image, Image.Image]:
    x1, y1, x2, y2 = [float(row[name]) for name in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")]
    width, height = image.size
    left = max(0, min(width - 1, math.floor(x1)))
    top = max(0, min(height - 1, math.floor(y1)))
    right = max(left + 1, min(width, math.ceil(x2)))
    bottom = max(top + 1, min(height, math.ceil(y2)))
    head = image.crop((left, top, right, bottom)).convert("RGB")
    box_width, box_height = right - left, bottom - top
    context = image.crop((
        max(0, math.floor(left - pad_ratio * box_width)),
        max(0, math.floor(top - pad_ratio * box_height)),
        min(width, math.ceil(right + pad_ratio * box_width)),
        min(height, math.ceil(bottom + pad_ratio * box_height)),
    )).convert("RGB")
    return head, context


def main() -> None:
    args = parse_args()
    if args.n_tasks_per_trait < 1:
        raise ValueError("--n-tasks-per-trait must be positive")
    if not 0 <= args.double_label_fraction <= 1:
        raise ValueError("--double-label-fraction must be in [0,1]")
    if not 0 <= args.context_pad_ratio <= 5:
        raise ValueError("--context-pad-ratio must be in [0,5]")

    state_column = f"analysis_state_ai_{args.measurement_mode}"
    head = pd.read_csv(args.head_traits, dtype=str, keep_default_na=False)
    yolo = pd.read_csv(args.yolo_crops, dtype=str, keep_default_na=False)
    require(head, {"annotation_unit_id", "trait_id", "taxon_name", "obs_id", state_column}, "head traits")
    require(yolo, {
        "queue_id", "obs_id", "photo_id", "det_index", "latitude", "longitude",
        "medium_image_url", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
    }, "YOLO crops")
    if head.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Head-trait rows are not unique")

    yolo = yolo.copy()
    yolo["annotation_unit_id"] = yolo.apply(annotation_unit_from_yolo, axis=1)
    if yolo["annotation_unit_id"].duplicated().any():
        raise ValueError("YOLO annotation-unit IDs are not unique")
    yolo["spatial_block_10deg"] = [
        block_10deg(latitude, longitude)
        for latitude, longitude in zip(yolo["latitude"], yolo["longitude"])
    ]

    trait_rows = head.loc[head["trait_id"].isin(args.traits)].copy()
    trait_rows["ai_candidate_state"] = trait_rows[state_column].map(text)
    trait_rows = trait_rows.loc[~trait_rows["ai_candidate_state"].isin(UNASSESSABLE)].copy()
    metadata_columns = [
        "annotation_unit_id", "obs_id", "photo_id", "latitude", "longitude",
        "spatial_block_10deg", "medium_image_url", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
    ]
    candidates = trait_rows.merge(
        yolo[metadata_columns], on=["annotation_unit_id", "obs_id"], how="inner", validate="many_to_one"
    )
    candidates = candidates.loc[candidates["spatial_block_10deg"].ne("")].copy()

    selected_parts: list[pd.DataFrame] = []
    per_trait: dict[str, Any] = {}
    for trait in args.traits:
        part = candidates.loc[candidates["trait_id"].eq(trait)].copy()
        if len(part) < args.n_tasks_per_trait:
            raise ValueError(f"Not enough candidates for {trait}: {len(part)}")
        selected = representative_sample(part, args.n_tasks_per_trait, f"{args.seed}:{trait}")
        selected_parts.append(selected)
        per_trait[trait] = {
            "n_candidates": int(len(part)),
            "n_selected": int(len(selected)),
            "n_taxa": int(selected["taxon_name"].nunique()),
            "n_spatial_blocks": int(selected["spatial_block_10deg"].nunique()),
        }

    selected = pd.concat(selected_parts, ignore_index=True)
    selected["task_id"] = [f"representative_holdout_{index:04d}" for index in range(1, len(selected) + 1)]
    selected["double_label"] = False
    for trait, group in selected.groupby("trait_id", sort=True):
        n_double = int(round(len(group) * args.double_label_fraction))
        ranked = group.assign(
            _rank=group["task_id"].map(lambda value: stable_rank(value, f"{args.seed}:{trait}:double"))
        ).sort_values("_rank")
        selected.loc[ranked.head(n_double).index, "double_label"] = True

    out = Path(args.out_dir)
    primary_root = out / "representative_trait_holdout_primary"
    secondary_root = out / "representative_trait_holdout_secondary"
    private_root = out / "representative_trait_holdout_key"
    for root in (primary_root, secondary_root):
        for directory in (root / "source_images", root / "head_crops", root / "context_crops"):
            directory.mkdir(parents=True, exist_ok=True)
    private_root.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "azami-representative-holdout/2.0"})
    paths_by_unit: dict[str, tuple[str, str, str]] = {}
    statuses: list[dict[str, str]] = []
    for _, row in selected.drop_duplicates("annotation_unit_id").iterrows():
        unit = text(row["annotation_unit_id"])
        photo_id = text(row["photo_id"])
        source_relative = f"source_images/inat_photo_{photo_id}.jpg"
        head_relative = f"head_crops/{unit}.jpg"
        context_relative = f"context_crops/{unit}.jpg"
        try:
            source_path = primary_root / source_relative
            if not source_path.exists():
                download_image(session, text(row["medium_image_url"]), source_path, args.timeout_sec, args.retries)
                if args.sleep_sec:
                    time.sleep(args.sleep_sec)
            with Image.open(source_path) as image:
                image.load()
                head_crop, context_crop = crop_box(image, row, args.context_pad_ratio)
            head_crop.save(primary_root / head_relative, format="JPEG", quality=94)
            context_crop.save(primary_root / context_relative, format="JPEG", quality=94)
            paths_by_unit[unit] = (source_relative, head_relative, context_relative)
            statuses.append({"annotation_unit_id": unit, "photo_id": photo_id, "status": "success"})
        except Exception as error:
            statuses.append({
                "annotation_unit_id": unit,
                "photo_id": photo_id,
                "status": "failed",
                "message": f"{type(error).__name__}: {error}",
            })

    status = pd.DataFrame(statuses)
    status.to_csv(out / "representative_holdout_image_build_status.csv", index=False, encoding="utf-8-sig")
    failed = status.loc[~status["status"].eq("success")]
    if not failed.empty:
        raise RuntimeError(f"Image reconstruction failed for {len(failed)} selected units")

    public_rows = []
    for row in selected.to_dict("records"):
        source_relative, head_relative, context_relative = paths_by_unit[text(row["annotation_unit_id"])]
        public_rows.append({
            "task_id": text(row["task_id"]),
            "annotation_unit_id": text(row["annotation_unit_id"]),
            "trait_id": text(row["trait_id"]),
            "source_image": source_relative,
            "crop_path": head_relative,
            "context_crop_path": context_relative,
        })
    public = pd.DataFrame(public_rows)
    for filename in ("representative_trait_holdout_tasks.csv", "blinded_trait_audit_tasks.csv"):
        public.to_csv(primary_root / filename, index=False, encoding="utf-8-sig")

    secondary_ids = set(selected.loc[selected["double_label"], "task_id"])
    secondary = public.loc[public["task_id"].isin(secondary_ids)].copy()
    for filename in ("representative_trait_holdout_tasks.csv", "blinded_trait_audit_tasks.csv"):
        secondary.to_csv(secondary_root / filename, index=False, encoding="utf-8-sig")
    for relative_path in sorted(set(secondary[["source_image", "crop_path", "context_crop_path"]].stack())):
        source = primary_root / relative_path
        target = secondary_root / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    private = selected[[
        "task_id", "annotation_unit_id", "trait_id", "taxon_name", "obs_id", "photo_id",
        "spatial_block_10deg", "double_label", "ai_candidate_state",
    ]].copy()
    private["double_label"] = private["double_label"].map(lambda value: str(bool(value)).lower())
    private.to_csv(private_root / "representative_trait_holdout_key.csv", index=False, encoding="utf-8-sig")

    readme = (
        "Chapter 1 representative trait audit.\n"
        "Primary packet: score every task. Secondary packet: score every task in that smaller packet.\n"
        "Use unassessable rather than guessing. Do not open or upload the private key.\n"
    )
    guide = (
        "capitulum_orientation: upright / ascending / horizontal / nodding / pendant / unassessable\n"
        "corolla_colour_class: white_or_cream / pale_pink / pink / purple_blue / yellow_or_other / mixed / unassessable\n"
        "capitulum_shape: globose / hemispherical / ovoid / cylindrical / urceolate / unassessable\n"
    )
    for root in (primary_root, secondary_root):
        (root / "README.txt").write_text(readme, encoding="utf-8")
        (root / "ANNOTATION_GUIDE.txt").write_text(guide, encoding="utf-8")
        if args.ontology:
            shutil.copy2(args.ontology, root / "ch1_trait_ontology.csv")

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "traits": args.traits,
        "n_tasks_per_trait": args.n_tasks_per_trait,
        "double_label_fraction": args.double_label_fraction,
        "n_primary_tasks": int(len(public)),
        "n_secondary_tasks": int(len(secondary)),
        "n_unique_heads": int(public["annotation_unit_id"].nunique()),
        "n_unique_photos": int(selected["photo_id"].nunique()),
        "per_trait": per_trait,
        "selection_design": "One proportional sample across taxon x 10-degree spatial cells; no species top-up.",
    }
    (out / "representative_trait_holdout_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build an app-ready representative human holdout from merged global AI outputs.

The input tables contain no retained image bytes, but do retain the exact medium
image URL and YOLO box coordinates used for inference. This script selects a
representative taxon x 10-degree-block holdout, re-downloads only selected source
images, recreates head/context crops, and writes public and private outputs to
separate directories.
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
    p = argparse.ArgumentParser()
    p.add_argument("--head-traits", required=True)
    p.add_argument("--yolo-crops", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--ontology", default=None, help="Optional ontology CSV copied into each public packet")
    p.add_argument("--traits", nargs="+", default=list(DEFAULT_TRAITS))
    p.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative")
    p.add_argument("--n-tasks-per-trait", type=int, default=120)
    p.add_argument("--double-label-fraction", type=float, default=0.25)
    p.add_argument("--context-pad-ratio", type=float, default=1.5)
    p.add_argument("--seed", default="20260710")
    p.add_argument("--sleep-sec", type=float, default=0.15)
    p.add_argument("--timeout-sec", type=float, default=60.0)
    p.add_argument("--retries", type=int, default=4)
    return p.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def stable_rank(value: str, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def block_10deg(latitude: Any, longitude: Any) -> str:
    try:
        lat = float(latitude)
        lon = float(longitude)
    except (TypeError, ValueError):
        return ""
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        return ""
    lat_index = min(17, max(0, math.floor((lat + 90.0) / 10.0)))
    lon_index = min(35, max(0, math.floor((lon + 180.0) / 10.0)))
    return f"lat{lat_index:02d}_lon{lon_index:02d}"


def representative_sample(part: pd.DataFrame, n_target: int, seed: str) -> pd.DataFrame:
    if part.empty or n_target <= 0:
        return part.head(0)
    part = part.copy()
    part["_rank"] = part["annotation_unit_id"].map(lambda v: stable_rank(text(v), seed))
    strata = [
        group.sort_values("_rank")
        for _, group in sorted(part.groupby(["taxon_name", "spatial_block_10deg"], sort=True), key=lambda kv: str(kv[0]))
    ]
    picked: list[pd.Series] = []
    cursor = 0
    limit = min(n_target, len(part))
    while len(picked) < limit:
        advanced = False
        for group in strata:
            if cursor < len(group):
                picked.append(group.iloc[cursor])
                advanced = True
                if len(picked) >= limit:
                    break
        if not advanced:
            break
        cursor += 1
    return pd.DataFrame(picked).drop(columns="_rank", errors="ignore")


def require(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    if missing := columns.difference(frame.columns):
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def annotation_unit_from_yolo(row: pd.Series) -> str:
    queue = text(row.get("audit_id")) or text(row.get("queue_id"))
    return f"{queue}_head_{int(float(row['det_index'])) + 1:02d}"


def download_image(session: requests.Session, urls: list[str], target: Path, timeout: float, retries: int, sleep_sec: float) -> str:
    last_error = ""
    for url in [u for u in urls if u]:
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
                return url
            except Exception as error:
                last_error = f"{type(error).__name__}: {error}"
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
                if attempt < retries:
                    time.sleep(min(8.0, 2.0 ** (attempt - 1)))
    raise RuntimeError(last_error or "No usable image URL")


def crop_box(image: Image.Image, row: pd.Series, pad_ratio: float) -> tuple[Image.Image, Image.Image]:
    x1, y1, x2, y2 = [float(row[name]) for name in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")]
    width, height = image.size
    left = max(0, min(width - 1, math.floor(x1)))
    top = max(0, min(height - 1, math.floor(y1)))
    right = max(left + 1, min(width, math.ceil(x2)))
    bottom = max(top + 1, min(height, math.ceil(y2)))
    head = image.crop((left, top, right, bottom)).convert("RGB")

    box_w, box_h = right - left, bottom - top
    c_left = max(0, math.floor(left - pad_ratio * box_w))
    c_top = max(0, math.floor(top - pad_ratio * box_h))
    c_right = min(width, math.ceil(right + pad_ratio * box_w))
    c_bottom = min(height, math.ceil(bottom + pad_ratio * box_h))
    context = image.crop((c_left, c_top, c_right, c_bottom)).convert("RGB")
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
        "medium_image_url", "small_image_url", "source_file_stem",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"
    }, "YOLO crops")
    if head.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Head-trait rows must be unique by annotation_unit_id/trait_id")

    yolo = yolo.copy()
    yolo["annotation_unit_id"] = yolo.apply(annotation_unit_from_yolo, axis=1)
    if yolo["annotation_unit_id"].duplicated().any():
        raise ValueError("Derived YOLO annotation_unit_id values are not unique")
    yolo["spatial_block_10deg"] = [block_10deg(a, b) for a, b in zip(yolo["latitude"], yolo["longitude"])]

    selected_head = head.loc[head["trait_id"].isin(args.traits)].copy()
    selected_head["ai_candidate_state"] = selected_head[state_column].map(text)
    selected_head = selected_head[["annotation_unit_id", "trait_id", "taxon_name", "obs_id", "ai_candidate_state"]]
    selected_head = selected_head.loc[~selected_head["ai_candidate_state"].isin(UNASSESSABLE)].copy()
    meta_columns = [
        "annotation_unit_id", "queue_id", "obs_id", "photo_id", "det_index", "latitude", "longitude",
        "spatial_block_10deg", "medium_image_url", "small_image_url", "source_file_stem",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"
    ]
    candidates = selected_head.merge(yolo[meta_columns], on=["annotation_unit_id", "obs_id"], how="inner", validate="many_to_one")
    candidates = candidates.loc[candidates["spatial_block_10deg"].ne("")].copy()
    if candidates.empty:
        raise ValueError("No representative holdout candidates remain")

    chosen_parts: list[pd.DataFrame] = []
    per_trait: dict[str, Any] = {}
    for trait in args.traits:
        part = candidates.loc[candidates["trait_id"].eq(trait)].copy()
        chosen = representative_sample(part, args.n_tasks_per_trait, f"{args.seed}:{trait}")
        chosen_parts.append(chosen)
        per_trait[trait] = {
            "n_candidates": int(len(part)),
            "n_selected": int(len(chosen)),
            "n_taxa": int(chosen["taxon_name"].nunique()),
            "n_blocks_10deg": int(chosen["spatial_block_10deg"].nunique()),
        }
    selected = pd.concat(chosen_parts, ignore_index=True)
    selected["task_id"] = [f"representative_holdout_{i:05d}" for i in range(1, len(selected) + 1)]
    selected["selection_stratum"] = "representative_taxon_x_10deg_block"
    selected["_double_rank"] = selected["task_id"].map(lambda v: stable_rank(v, args.seed))
    n_double = int(round(len(selected) * args.double_label_fraction))
    double_ids = set(selected.sort_values("_double_rank").head(n_double)["task_id"])
    selected["double_label"] = selected["task_id"].isin(double_ids)

    out = Path(args.out_dir)
    public_root = out / "representative_trait_holdout_primary"
    secondary_root = out / "representative_trait_holdout_secondary"
    private_root = out / "representative_trait_holdout_key"
    for root in (public_root, secondary_root):
        for directory in (root / "source_images", root / "head_crops", root / "context_crops"):
            directory.mkdir(parents=True, exist_ok=True)
    private_root.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "azami-representative-holdout/1.0"})
    unique_units = selected.drop_duplicates("annotation_unit_id")
    source_by_photo: dict[str, str] = {}
    image_records: list[dict[str, Any]] = []
    public_paths: dict[str, tuple[str, str, str]] = {}

    for _, row in unique_units.iterrows():
        photo_id = text(row["photo_id"])
        source_rel = f"source_images/inat_photo_{photo_id}.jpg"
        source_path = public_root / source_rel
        try:
            if not source_path.exists():
                used_url = download_image(
                    session,
                    [text(row["medium_image_url"])],
                    source_path,
                    args.timeout_sec,
                    args.retries,
                    args.sleep_sec,
                )
                source_by_photo[photo_id] = used_url
                if args.sleep_sec:
                    time.sleep(args.sleep_sec)
            with Image.open(source_path) as image:
                image.load()
                head_crop, context_crop = crop_box(image, row, args.context_pad_ratio)
            unit = text(row["annotation_unit_id"])
            head_rel = f"head_crops/{unit}.jpg"
            context_rel = f"context_crops/{unit}.jpg"
            head_crop.save(public_root / head_rel, format="JPEG", quality=94)
            context_crop.save(public_root / context_rel, format="JPEG", quality=94)
            public_paths[unit] = (source_rel, head_rel, context_rel)
            image_records.append({
                "annotation_unit_id": unit, "photo_id": photo_id, "status": "success",
                "source_url_used": source_by_photo.get(photo_id, "cached"),
                "source_image": source_rel, "crop_path": head_rel, "context_crop_path": context_rel,
            })
        except Exception as error:
            image_records.append({
                "annotation_unit_id": text(row["annotation_unit_id"]), "photo_id": photo_id,
                "status": "failed", "source_url_used": "", "source_image": source_rel,
                "crop_path": "", "context_crop_path": "", "message": f"{type(error).__name__}: {error}",
            })

    image_status = pd.DataFrame(image_records)
    image_status.to_csv(out / "representative_holdout_image_build_status.csv", index=False, encoding="utf-8-sig")
    failed_units = set(image_status.loc[~image_status["status"].eq("success"), "annotation_unit_id"])
    if failed_units:
        raise RuntimeError(
            f"Image reconstruction failed for {len(failed_units)} selected units; "
            "refusing a post-selection-biased partial holdout. See representative_holdout_image_build_status.csv."
        )

    public_rows: list[dict[str, str]] = []
    for row in selected.to_dict("records"):
        source_rel, head_rel, context_rel = public_paths[text(row["annotation_unit_id"])]
        public_rows.append({
            "task_id": text(row["task_id"]),
            "annotation_unit_id": text(row["annotation_unit_id"]),
            "trait_id": text(row["trait_id"]),
            "source_image": source_rel,
            "crop_path": head_rel,
            "context_crop_path": context_rel,
        })
    public = pd.DataFrame(public_rows)
    public.to_csv(public_root / "representative_trait_holdout_tasks.csv", index=False, encoding="utf-8-sig")
    public.to_csv(public_root / "blinded_trait_audit_tasks.csv", index=False, encoding="utf-8-sig")

    secondary_ids = set(selected.loc[selected["double_label"], "task_id"])
    secondary = public.loc[public["task_id"].isin(secondary_ids)].copy()
    secondary.to_csv(secondary_root / "representative_trait_holdout_tasks.csv", index=False, encoding="utf-8-sig")
    secondary.to_csv(secondary_root / "blinded_trait_audit_tasks.csv", index=False, encoding="utf-8-sig")
    for relative_path in sorted(set(secondary[["source_image", "crop_path", "context_crop_path"]].stack())):
        source = public_root / relative_path
        target = secondary_root / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)

    private_columns = [
        "task_id", "annotation_unit_id", "trait_id", "taxon_name", "obs_id", "photo_id",
        "spatial_block_10deg", "selection_stratum", "double_label", "ai_candidate_state"
    ]
    private = selected[private_columns].copy()
    private["double_label"] = private["double_label"].map(lambda v: str(bool(v)).lower())
    private.to_csv(private_root / "representative_trait_holdout_key.csv", index=False, encoding="utf-8-sig")

    readme = (
        "Representative blinded holdout for Cirsium Chapter 1.\n"
        "Score visible evidence only. AI candidate, taxon, locality and selection stratum are hidden.\n"
        "Use unassessable whenever the required visual context does not support a decision.\n"
        "Do not combine this public directory with representative_trait_holdout_key.\n"
    )
    guide_ja = (
        "主要形質の判定ガイド\n\n"
        "capitulum_orientation（頭花の向き）: upright=直立、ascending=斜上、"
        "horizontal=水平、nodding=点頭、pendant=下垂。茎・花柄が見えなければunassessable。\n"
        "corolla_colour_class（花色）: white_or_cream=白〜クリーム、pale_pink=淡桃、"
        "pink=桃、purple_blue=紫〜青紫、yellow_or_other=黄その他、mixed=明瞭な混在。"
        "逆光や色かぶりが強ければunassessable。\n"
        "capitulum_shape（頭花外形）: globose=球形、hemispherical=半球形、"
        "ovoid=卵形、cylindrical=円筒形、urceolate=壺形。正面視や輪郭不明ならunassessable。\n"
    )
    for root in (public_root, secondary_root):
        (root / "README.txt").write_text(readme, encoding="utf-8")
        (root / "ANNOTATION_GUIDE_JA.txt").write_text(guide_ja, encoding="utf-8")
        if args.ontology:
            shutil.copy2(args.ontology, root / "ch1_trait_ontology.csv")

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "state_column": state_column,
        "traits": args.traits,
        "seed": args.seed,
        "n_tasks_requested_per_trait": args.n_tasks_per_trait,
        "double_label_fraction": args.double_label_fraction,
        "context_pad_ratio": args.context_pad_ratio,
        "n_candidate_rows": int(len(candidates)),
        "n_tasks_selected_before_image_build": int(sum(v["n_selected"] for v in per_trait.values())),
        "n_tasks_public": int(len(public)),
        "n_unique_units_public": int(public["annotation_unit_id"].nunique()),
        "n_unique_photos_public": int(selected["photo_id"].nunique()),
        "n_secondary_double_label_tasks": int(len(secondary)),
        "n_image_units_failed": int(len(failed_units)),
        "per_trait": {
            trait: {
                **per_trait[trait],
                "n_public_after_image_build": int(public["trait_id"].eq(trait).sum()),
                "n_public_taxa": int(selected.loc[selected["trait_id"].eq(trait), "taxon_name"].nunique()),
                "n_public_blocks_10deg": int(selected.loc[selected["trait_id"].eq(trait), "spatial_block_10deg"].nunique()),
            }
            for trait in args.traits
        },
        "selection_design": "Deterministic round-robin sample across taxon x 10-degree spatial block; model margin is not used.",
        "image_reconstruction": "Re-downloaded the same iNaturalist medium image size used by global inference and recreated YOLO head/context crops from stored boxes.",
    }
    (out / "representative_trait_holdout_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

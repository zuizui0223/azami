#!/usr/bin/env python3
"""Rebuild all global head/context crops from retained YOLO metadata.

The original shard workflows intentionally uploaded compact metadata rather than
thousands of image files. This script re-downloads each retained iNaturalist
medium image once, recreates all stored detector boxes, writes only head/context
crops, and emits a manifest for continuous trait measurement.
"""
from __future__ import annotations

import argparse
import io
import json
import math
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from PIL import Image


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--yolo-metadata", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--context-pad-ratio", type=float, default=1.5)
    p.add_argument("--timeout-sec", type=float, default=60.0)
    p.add_argument("--retries", type=int, default=4)
    p.add_argument("--sleep-sec", type=float, default=0.10)
    p.add_argument("--max-photos", type=int, default=0)
    p.add_argument("--allow-failures", action="store_true")
    return p.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def annotation_unit_id(row: pd.Series | dict[str, Any]) -> str:
    queue = text(row.get("audit_id")) or text(row.get("queue_id"))
    if not queue:
        raise ValueError("YOLO row lacks audit_id/queue_id")
    return f"{queue}_head_{int(float(row['det_index'])) + 1:02d}"


def clamp_box(
    width: int,
    height: int,
    x1: Any,
    y1: Any,
    x2: Any,
    y2: Any,
) -> tuple[int, int, int, int]:
    left = max(0, min(width - 1, math.floor(float(x1))))
    top = max(0, min(height - 1, math.floor(float(y1))))
    right = max(left + 1, min(width, math.ceil(float(x2))))
    bottom = max(top + 1, min(height, math.ceil(float(y2))))
    return left, top, right, bottom


def crop_pair(
    image: Image.Image,
    row: pd.Series | dict[str, Any],
    pad_ratio: float,
) -> tuple[Image.Image, Image.Image]:
    width, height = image.size
    left, top, right, bottom = clamp_box(
        width,
        height,
        row["bbox_x1"],
        row["bbox_y1"],
        row["bbox_x2"],
        row["bbox_y2"],
    )
    head = image.crop((left, top, right, bottom)).convert("RGB")
    box_width = right - left
    box_height = bottom - top
    context_left = max(0, math.floor(left - pad_ratio * box_width))
    context_top = max(0, math.floor(top - pad_ratio * box_height))
    context_right = min(width, math.ceil(right + pad_ratio * box_width))
    context_bottom = min(height, math.ceil(bottom + pad_ratio * box_height))
    context = image.crop(
        (context_left, context_top, context_right, context_bottom)
    ).convert("RGB")
    return head, context


def download_image(
    session: requests.Session,
    url: str,
    timeout: float,
    retries: int,
) -> Image.Image:
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            with Image.open(io.BytesIO(response.content)) as opened:
                opened.verify()
            with Image.open(io.BytesIO(response.content)) as opened:
                return opened.convert("RGB")
        except Exception as error:
            last_error = f"{type(error).__name__}: {error}"
            if attempt < retries:
                time.sleep(min(8.0, 2.0 ** (attempt - 1)))
    raise RuntimeError(last_error or "Image download failed")


def write_checkpoint(path: Path, rows: list[dict[str, Any]]) -> None:
    frame = pd.DataFrame(rows)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, encoding="utf-8-sig")
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    if not 0 <= args.context_pad_ratio <= 5:
        raise ValueError("--context-pad-ratio must be in [0,5]")
    if args.retries < 1:
        raise ValueError("--retries must be positive")
    metadata = pd.read_csv(
        args.yolo_metadata, dtype=str, keep_default_na=False
    )
    required = {
        "queue_id", "obs_id", "photo_id", "det_index", "medium_image_url",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
    }
    missing = required.difference(metadata.columns)
    if missing or metadata.empty:
        raise ValueError(f"YOLO metadata are invalid; missing={sorted(missing)}")
    metadata = metadata.copy()
    metadata["annotation_unit_id"] = metadata.apply(annotation_unit_id, axis=1)
    if metadata["annotation_unit_id"].duplicated().any():
        raise ValueError("Derived annotation_unit_id values are not unique")

    photo_groups = list(metadata.groupby("photo_id", sort=True))
    if args.max_photos > 0:
        photo_groups = photo_groups[: args.max_photos]

    out = Path(args.out_dir)
    head_dir = out / "head_crops"
    context_dir = out / "context_crops"
    head_dir.mkdir(parents=True, exist_ok=True)
    context_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = out / "rebuild_status.csv"

    session = requests.Session()
    session.headers.update({"User-Agent": "azami-continuous-traits/1.0"})
    status_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []

    for photo_index, (photo_id, group) in enumerate(photo_groups, start=1):
        url_values = [
            value for value in group["medium_image_url"].map(text).unique()
            if value
        ]
        try:
            if len(url_values) != 1:
                raise ValueError(
                    f"photo_id={photo_id} has {len(url_values)} medium URLs"
                )
            image = download_image(
                session,
                url_values[0],
                args.timeout_sec,
                args.retries,
            )
            for _, row in group.iterrows():
                unit = row["annotation_unit_id"]
                head_relative = f"head_crops/{unit}.jpg"
                context_relative = f"context_crops/{unit}.jpg"
                head_path = out / head_relative
                context_path = out / context_relative
                if not head_path.is_file() or not context_path.is_file():
                    head, context = crop_pair(
                        image, row, args.context_pad_ratio
                    )
                    head.save(head_path, format="JPEG", quality=94)
                    context.save(context_path, format="JPEG", quality=94)
                manifest_rows.append({
                    "annotation_unit_id": unit,
                    "crop_path": head_relative,
                    "context_crop_path": context_relative,
                })
                status_rows.append({
                    "photo_id": photo_id,
                    "annotation_unit_id": unit,
                    "status": "success",
                    "message": "",
                })
        except Exception as error:
            for _, row in group.iterrows():
                status_rows.append({
                    "photo_id": photo_id,
                    "annotation_unit_id": row["annotation_unit_id"],
                    "status": "failed",
                    "message": f"{type(error).__name__}: {error}",
                })
        if args.sleep_sec:
            time.sleep(args.sleep_sec)
        if photo_index % 25 == 0 or photo_index == len(photo_groups):
            write_checkpoint(checkpoint_path, status_rows)
            print(
                f"[INFO] rebuilt photos {photo_index}/{len(photo_groups)}; "
                f"heads={len(manifest_rows)}"
            )

    status = pd.DataFrame(status_rows)
    failures = status.loc[~status["status"].eq("success")]
    manifest = (
        pd.DataFrame(manifest_rows)
        .drop_duplicates("annotation_unit_id")
        .sort_values("annotation_unit_id")
    )
    manifest.to_csv(
        out / "continuous_trait_manifest.csv",
        index=False,
        encoding="utf-8-sig",
    )
    report = {
        "n_input_heads": int(sum(len(group) for _, group in photo_groups)),
        "n_input_photos": int(len(photo_groups)),
        "n_rebuilt_heads": int(len(manifest)),
        "n_failed_heads": int(len(failures)),
        "context_pad_ratio": args.context_pad_ratio,
        "semantic_status": (
            "Reconstructed head/context crops from retained medium-image URLs "
            "and stored YOLO boxes; no detector was rerun."
        ),
    }
    (out / "continuous_trait_rebuild_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if len(failures) and not args.allow_failures:
        raise RuntimeError(
            f"Reconstruction failed for {len(failures)} heads; "
            "see rebuild_status.csv"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build wide head-plus-peduncle context crops from YOLO detection metadata.

Direction cannot be reliably scored from a tight head crop. This script keeps
existing detector crops unchanged and creates a second, wider crop for each
valid detection. It does not infer orientation.

When detector outputs were restored from an artifact, their original absolute
source paths are stale. Supply --images-dir to resolve them by the immutable
screen download filename (preferred) or source-image basename.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import pandas as pd
from PIL import Image

REQUIRED = {"queue_id", "source_image", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create head-plus-peduncle context crops from YOLO boxes.")
    parser.add_argument("--crops", required=True, help="yolo_crop_metadata.csv")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--images-dir", default="", help="Optional restored source-image directory.")
    parser.add_argument("--context-pad-ratio", type=float, default=0.8)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def resolve_source_image(record: dict[str, Any], images_dir: Path | None) -> Path:
    original = Path(text(record.get("source_image", "")))
    if original.is_file():
        return original
    if images_dir is None:
        return original
    candidates = []
    filename = text(record.get("screen_download_filename", ""))
    if filename:
        candidates.append(images_dir / filename)
    if original.name:
        candidates.append(images_dir / original.name)
    image_local = Path(text(record.get("image_local_path", "")))
    if image_local.name:
        candidates.append(images_dir / image_local.name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return original


def safe_box(x1: float, y1: float, x2: float, y2: float, width: int, height: int, pad_ratio: float) -> tuple[int, int, int, int]:
    box_width = max(0.0, x2 - x1)
    box_height = max(0.0, y2 - y1)
    left = max(0, int(x1 - box_width * pad_ratio))
    top = max(0, int(y1 - box_height * pad_ratio))
    right = min(width, int(x2 + box_width * pad_ratio))
    bottom = min(height, int(y2 + box_height * pad_ratio))
    if right <= left or bottom <= top:
        raise ValueError("non-positive crop extent")
    return left, top, right, bottom


def atomic_save(image: Image.Image, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.stem}.tmp{path.suffix}")
    image.save(temporary)
    os.replace(temporary, path)


def detection_number(value: Any) -> int:
    try:
        return int(float(text(value))) + 1
    except ValueError:
        return 1


def main() -> None:
    args = parse_args()
    if args.context_pad_ratio < 0:
        raise ValueError("--context-pad-ratio must be >= 0")
    images_dir = Path(args.images_dir) if text(args.images_dir) else None
    if images_dir is not None and not images_dir.is_dir():
        raise FileNotFoundError(f"--images-dir is not a directory: {images_dir}")
    crops = pd.read_csv(args.crops, dtype=str, keep_default_na=False)
    missing = REQUIRED.difference(crops.columns)
    if missing:
        raise ValueError(f"YOLO crop metadata missing required columns: {sorted(missing)}")
    out_dir = Path(args.out_dir)
    context_dir = out_dir / "head_peduncle_context_crops"
    rows: list[dict[str, object]] = []
    for _, row in crops.iterrows():
        record = row.to_dict()
        source = resolve_source_image(record, images_dir)
        identifier = text(row.get("source_file_stem", "")) or text(row["queue_id"])
        det_number = detection_number(row.get("det_index", "0"))
        target = context_dir / f"{identifier}_context_det_{det_number:02d}.jpg"
        try:
            coordinates = [float(row[name]) for name in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")]
            with Image.open(source) as opened:
                image = opened.convert("RGB")
            crop_box = safe_box(*coordinates, image.width, image.height, args.context_pad_ratio)
            if args.force or not target.exists() or target.stat().st_size < 1000:
                atomic_save(image.crop(crop_box), target)
            record["context_crop_status"] = "success"
            record["resolved_source_image"] = str(source.resolve())
            record["context_crop_path"] = str(target.resolve())
            record["context_pad_ratio"] = args.context_pad_ratio
            record["context_crop_box"] = ",".join(map(str, crop_box))
            record["context_crop_message"] = ""
        except Exception as error:
            record["context_crop_status"] = "failed"
            record["resolved_source_image"] = str(source)
            record["context_crop_path"] = ""
            record["context_pad_ratio"] = args.context_pad_ratio
            record["context_crop_box"] = ""
            record["context_crop_message"] = str(error)
        rows.append(record)
    result = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_dir / "head_context_crop_metadata.csv", index=False, encoding="utf-8-sig")
    n_success = int(result["context_crop_status"].eq("success").sum())
    print(f"[OK] wrote {len(result)} context-crop rows ({n_success} successful) to {out_dir}")


if __name__ == "__main__":
    main()

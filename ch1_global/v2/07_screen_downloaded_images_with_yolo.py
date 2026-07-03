#!/usr/bin/env python3
"""Run the existing Cirsium head detector on a v2 image-screening queue.

This is deliberately detector-only: it produces head crops and detector
metadata but does *not* run the legacy upward/nodding classifier. Human image-
level annotations, not inherited species labels, become the v2 trait ground
truth after this stage.

Example
-------
python ch1_global/v2/07_screen_downloaded_images_with_yolo.py \
  --queue "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_screen_queue\\image_screening_queue.csv" \
  --downloads "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_screen_images\\download_results.csv" \
  --weights "C:\\Users\\zuizui\\mc\\runs\\cirsium_detect\\weights\\best.pt" \
  --out-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_yolo"
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
from typing import Any

import pandas as pd
from PIL import Image


SCREEN_VERSION = "1.0.0"
RESULT_FIELDS = ["queue_id", "screen_status", "n_detections", "max_yolo_conf", "source_image", "message"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Screen downloaded Cirsium images using a detector only.")
    parser.add_argument("--queue", required=True)
    parser.add_argument("--downloads", required=True)
    parser.add_argument("--weights", required=True, help="Ultralytics YOLO head detector weights (.pt)")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--conf", type=float, default=0.25)
    parser.add_argument("--pad-ratio", type=float, default=0.12)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--device", default="", help="Optional Ultralytics device, e.g. 0 or cpu")
    parser.add_argument("--force", action="store_true", help="Re-screen terminal image rows already present in screening_results.csv")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def safe_crop_box(x1: float, y1: float, x2: float, y2: float, width: int, height: int, pad_ratio: float) -> tuple[int, int, int, int]:
    box_width, box_height = max(0.0, x2 - x1), max(0.0, y2 - y1)
    pad_x, pad_y = box_width * pad_ratio, box_height * pad_ratio
    left = max(0, int(x1 - pad_x))
    top = max(0, int(y1 - pad_y))
    right = min(width, int(x2 + pad_x))
    bottom = min(height, int(y2 + pad_y))
    if right <= left or bottom <= top:
        raise ValueError("Detector box has non-positive cropped extent")
    return left, top, right, bottom


def append_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def load_terminal_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    try:
        results = pd.read_csv(path, dtype=str)
    except pd.errors.EmptyDataError:
        return set()
    terminal = {"detected", "no_detection", "missing_image"}
    return set(results.loc[results["screen_status"].isin(terminal), "queue_id"].dropna().astype(str))


def latest_successful_downloads(path: Path) -> pd.DataFrame:
    downloads = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"queue_id", "download_status", "image_local_path"}
    missing = required.difference(downloads.columns)
    if missing:
        raise ValueError(f"Download results missing columns: {sorted(missing)}")
    downloads["_row_order"] = range(len(downloads))
    latest = downloads.sort_values("_row_order").groupby("queue_id", as_index=False).tail(1)
    return latest.loc[latest["download_status"].eq("success")].copy()


def atomic_save_crop(image: Image.Image, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.stem}.tmp{path.suffix}")
    image.save(temporary)
    os.replace(temporary, path)


def yolo_rows_for_image(base: dict[str, Any], result: Any, crop_dir: Path, pad_ratio: float, model_names: Any) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    image_path = Path(text(base["image_local_path"]))
    try:
        with Image.open(image_path) as opened:
            source = opened.convert("RGB")
    except Exception as error:
        return ({"queue_id": base["queue_id"], "screen_status": "missing_image", "n_detections": 0, "max_yolo_conf": "", "source_image": str(image_path), "message": str(error)}, [])

    boxes = getattr(result, "boxes", None)
    if boxes is None or len(boxes) == 0:
        return ({"queue_id": base["queue_id"], "screen_status": "no_detection", "n_detections": 0, "max_yolo_conf": "", "source_image": str(image_path), "message": ""}, [])
    xyxy = boxes.xyxy.cpu().tolist()
    confidences = boxes.conf.cpu().tolist()
    classes = boxes.cls.cpu().tolist()
    crop_rows: list[dict[str, Any]] = []
    for index, (coordinates, confidence, class_id) in enumerate(zip(xyxy, confidences, classes), start=1):
        x1, y1, x2, y2 = coordinates
        try:
            left, top, right, bottom = safe_crop_box(x1, y1, x2, y2, source.width, source.height, pad_ratio)
        except ValueError:
            continue
        crop_path = crop_dir / f"{text(base['source_file_stem'])}_det_{index:02d}.jpg"
        atomic_save_crop(source.crop((left, top, right, bottom)), crop_path)
        class_int = int(class_id)
        if isinstance(model_names, dict):
            class_name = text(model_names.get(class_int, class_int))
        elif isinstance(model_names, list) and 0 <= class_int < len(model_names):
            class_name = text(model_names[class_int])
        else:
            class_name = str(class_int)
        crop_rows.append({
            **base,
            "det_index": index - 1,
            "yolo_class_id": class_int,
            "yolo_class_name": class_name,
            "yolo_conf": float(confidence),
            "bbox_x1": float(x1), "bbox_y1": float(y1), "bbox_x2": float(x2), "bbox_y2": float(y2),
            "crop_path": str(crop_path.resolve()),
            "source_image": str(image_path.resolve()),
        })
    n_detections = len(crop_rows)
    return (
        {
            "queue_id": base["queue_id"],
            "screen_status": "detected" if n_detections else "no_detection",
            "n_detections": n_detections,
            "max_yolo_conf": max([row["yolo_conf"] for row in crop_rows], default=""),
            "source_image": str(image_path),
            "message": "",
        },
        crop_rows,
    )


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.conf <= 1.0:
        raise ValueError("--conf must be in [0, 1]")
    if args.pad_ratio < 0 or args.batch_size < 1:
        raise ValueError("--pad-ratio must be >= 0 and --batch-size must be >= 1")
    weights = Path(args.weights)
    if not weights.exists():
        raise FileNotFoundError(f"YOLO weights not found: {weights}")

    try:
        from ultralytics import YOLO
    except ImportError as error:
        raise SystemExit("Ultralytics is required: python -m pip install ultralytics") from error

    queue = pd.read_csv(args.queue, dtype=str, keep_default_na=False)
    required = {"queue_id", "photo_id", "source_file_stem"}
    missing = required.difference(queue.columns)
    if missing:
        raise ValueError(f"Queue missing columns: {sorted(missing)}")
    if queue["queue_id"].duplicated().any():
        raise ValueError("Queue has duplicate queue_id values")
    successes = latest_successful_downloads(Path(args.downloads))
    candidates = queue.merge(successes[["queue_id", "image_local_path"]], on="queue_id", how="inner", validate="one_to_one")

    out_dir = Path(args.out_dir)
    crop_dir = out_dir / "head_crops"
    results_path = out_dir / "screening_results.csv"
    crop_meta_path = out_dir / "yolo_crop_metadata.csv"
    terminal = set() if args.force else load_terminal_ids(results_path)
    candidates = candidates.loc[~candidates["queue_id"].astype(str).isin(terminal)].copy()
    if candidates.empty:
        print("[INFO] no unprocessed successful downloads to screen")
        return

    model = YOLO(str(weights))
    model_names = getattr(model, "names", {})
    print(f"[INFO] screening {len(candidates)} downloaded images with detector only")
    for start in range(0, len(candidates), args.batch_size):
        batch = candidates.iloc[start:start + args.batch_size]
        paths = [text(value) for value in batch["image_local_path"]]
        try:
            kwargs: dict[str, Any] = {"conf": args.conf, "verbose": False}
            if args.device:
                kwargs["device"] = args.device
            predictions = model.predict(paths, **kwargs)
        except Exception as error:
            failed = [{"queue_id": row["queue_id"], "screen_status": "error", "n_detections": 0, "max_yolo_conf": "", "source_image": row["image_local_path"], "message": str(error)} for _, row in batch.iterrows()]
            append_csv(results_path, failed, RESULT_FIELDS)
            continue

        screen_rows: list[dict[str, Any]] = []
        crop_rows: list[dict[str, Any]] = []
        for (_, row), prediction in zip(batch.iterrows(), predictions):
            screen, crops = yolo_rows_for_image(row.to_dict(), prediction, crop_dir, args.pad_ratio, model_names)
            screen_rows.append(screen)
            crop_rows.extend(crops)
        append_csv(results_path, screen_rows, RESULT_FIELDS)
        if crop_rows:
            append_csv(crop_meta_path, crop_rows, list(crop_rows[0].keys()))
        print(f"[INFO] screened {min(start + len(batch), len(candidates))}/{len(candidates)} images")

    print("[OK] wrote detector screening outputs", out_dir.resolve())


if __name__ == "__main__":
    main()

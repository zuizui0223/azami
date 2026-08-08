#!/usr/bin/env python3
"""Build a blinded independent validation packet for continuous image traits.

Selection is fixed before inspecting production trait values: one deterministic head
per photograph is chosen, candidates are spread across taxon x 10-degree blocks,
and download failures are replaced only by the next pre-ranked reserve candidate.
The public packet contains images and opaque task IDs only. Taxonomy, coordinates,
production measurements and source URLs remain in a separate private key.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import shutil
import time
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd
from PIL import Image, ImageOps

AUTO_COLUMNS = [
    "colour_status", "corolla_lab_lightness", "corolla_lab_chroma",
    "corolla_hue_degrees", "corolla_hue_sin", "corolla_hue_cos",
    "shape_status", "shape_aspect_ratio", "shape_circularity",
    "shape_solidity", "shape_width_cv", "orientation_status",
    "orientation_angle_degrees",
]


def stable_hash(text: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}|{text}".encode()).hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--head-table", required=True, type=Path)
    p.add_argument("--yolo-metadata", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--n-tasks", type=int, default=180)
    p.add_argument("--n-double", type=int, default=45)
    p.add_argument("--reserve", type=int, default=240)
    p.add_argument("--seed", type=int, default=20260807)
    p.add_argument("--context-pad-ratio", type=float, default=1.5)
    p.add_argument("--timeout", type=float, default=45)
    p.add_argument("--retries", type=int, default=4)
    return p.parse_args()


def annotation_id_from_yolo(row: pd.Series) -> str:
    return f"{row['queue_id']}_head_{int(row['det_index']) + 1:02d}"


def download(url: str, timeout: float, retries: int) -> bytes:
    last: Exception | None = None
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "azami-ch1-independent-trait-validation/1.0"},
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except Exception as error:
            last = error
            time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"download failed after {retries} attempts: {last}")


def crop_boxes(row: pd.Series, width: int, height: int, pad_ratio: float) -> tuple[tuple[int,int,int,int], tuple[int,int,int,int]]:
    x1, y1, x2, y2 = [float(row[c]) for c in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")]
    x1 = max(0.0, min(x1, width - 1.0)); x2 = max(1.0, min(x2, float(width)))
    y1 = max(0.0, min(y1, height - 1.0)); y2 = max(1.0, min(y2, float(height)))
    if x2 <= x1 + 2 or y2 <= y1 + 2:
        raise ValueError("invalid detector box")
    tight = (int(math.floor(x1)), int(math.floor(y1)), int(math.ceil(x2)), int(math.ceil(y2)))
    bw, bh = x2 - x1, y2 - y1
    cx, cy = (x1+x2)/2, (y1+y2)/2
    cw, ch = bw * pad_ratio, bh * pad_ratio
    context = (
        max(0, int(math.floor(cx - cw/2))),
        max(0, int(math.floor(cy - ch/2))),
        min(width, int(math.ceil(cx + cw/2))),
        min(height, int(math.ceil(cy + ch/2))),
    )
    return tight, context


def main() -> int:
    args = parse_args()
    if not (0 < args.n_double <= args.n_tasks <= args.reserve):
        raise ValueError("Require 0 < n_double <= n_tasks <= reserve")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    blinded = args.out_dir / "blinded_packet"
    private = args.out_dir / "private_key"
    images = blinded / "images"
    images.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)

    heads = pd.read_csv(args.head_table, low_memory=False)
    meta = pd.read_csv(args.yolo_metadata, low_memory=False)
    required_head = {"annotation_unit_id", "obs_id", "photo_id", "taxon_name", "latitude", "longitude", *AUTO_COLUMNS}
    required_meta = {"queue_id", "det_index", "obs_id", "photo_id", "taxon_name", "spatial_block", "medium_image_url", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "yolo_conf"}
    mh = required_head.difference(heads.columns); mm = required_meta.difference(meta.columns)
    if mh: raise ValueError(f"head table missing: {sorted(mh)}")
    if mm: raise ValueError(f"YOLO metadata missing: {sorted(mm)}")
    if heads["annotation_unit_id"].duplicated().any(): raise ValueError("duplicate head IDs")
    meta = meta.copy()
    meta["annotation_unit_id"] = meta.apply(annotation_id_from_yolo, axis=1)
    if meta["annotation_unit_id"].duplicated().any(): raise ValueError("duplicate YOLO head IDs")

    keep_meta = ["annotation_unit_id", "queue_id", "det_index", "spatial_block", "medium_image_url", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "yolo_conf"]
    merged = heads.merge(meta[keep_meta], on="annotation_unit_id", how="inner", validate="one_to_one")
    if len(merged) != len(heads):
        raise ValueError(f"YOLO metadata matched {len(merged)} / {len(heads)} heads")

    # One head per photograph BEFORE validation selection prevents pseudo-replication.
    merged["head_hash"] = merged["annotation_unit_id"].map(lambda x: stable_hash(str(x), args.seed + 1))
    one = merged.sort_values("head_hash").drop_duplicates("photo_id", keep="first").copy()
    one["stratum"] = one["taxon_name"].astype(str) + "|" + one["spatial_block"].astype(str)
    one["within_stratum_rank"] = one.groupby("stratum")["head_hash"].rank(method="first").astype(int)
    one["selection_hash"] = one["annotation_unit_id"].map(lambda x: stable_hash(str(x), args.seed + 2))
    one = one.sort_values(["within_stratum_rank", "selection_hash"]).head(args.reserve).reset_index(drop=True)
    one["pre_download_rank"] = range(1, len(one)+1)

    audit_rows: list[dict[str, Any]] = []
    successful: list[tuple[pd.Series, dict[str, Any]]] = []
    for _, row in one.iterrows():
        record: dict[str, Any] = {
            "annotation_unit_id": row["annotation_unit_id"],
            "pre_download_rank": int(row["pre_download_rank"]),
            "download_status": "pending",
            "download_error": "",
        }
        try:
            blob = download(str(row["medium_image_url"]), args.timeout, args.retries)
            image = Image.open(io.BytesIO(blob))
            image = ImageOps.exif_transpose(image).convert("RGB")
            width, height = image.size
            tight_box, context_box = crop_boxes(row, width, height, args.context_pad_ratio)
            tight = image.crop(tight_box)
            context = image.crop(context_box)
            if min(tight.size) < 20 or min(context.size) < 20:
                raise ValueError("crop too small")
            successful.append((row, {
                "source_image": image,
                "tight_image": tight,
                "context_image": context,
                "source_width": width,
                "source_height": height,
                "tight_box": tight_box,
                "context_box": context_box,
            }))
            record["download_status"] = "success"
        except Exception as error:
            record["download_status"] = "failed"
            record["download_error"] = f"{type(error).__name__}: {error}"
        audit_rows.append(record)
        if len(successful) >= args.n_tasks:
            break

    if len(successful) < args.n_tasks:
        raise RuntimeError(f"Only {len(successful)} of {args.n_tasks} required validation images reconstructed")

    public_rows = []
    private_rows = []
    for index, (row, payload) in enumerate(successful[:args.n_tasks], start=1):
        task_id = f"traitval_{index:04d}"
        source_name = f"{task_id}_source.jpg"
        context_name = f"{task_id}_context.jpg"
        tight_name = f"{task_id}_tight.jpg"
        payload["source_image"].save(images/source_name, quality=94)
        payload["context_image"].save(images/context_name, quality=94)
        payload["tight_image"].save(images/tight_name, quality=94)
        public_rows.append({
            "task_id": task_id,
            "source_image": f"images/{source_name}",
            "context_image": f"images/{context_name}",
            "tight_image": f"images/{tight_name}",
        })
        private_row = {
            "task_id": task_id,
            "annotation_unit_id": row["annotation_unit_id"],
            "obs_id": row["obs_id"], "photo_id": row["photo_id"],
            "taxon_name": row["taxon_name"], "latitude": row["latitude"], "longitude": row["longitude"],
            "spatial_block": row["spatial_block"], "medium_image_url": row["medium_image_url"],
            "yolo_conf": row["yolo_conf"], "pre_download_rank": int(row["pre_download_rank"]),
            "source_width": payload["source_width"], "source_height": payload["source_height"],
            "bbox_source_xyxy": json.dumps(payload["tight_box"]),
            "context_source_xyxy": json.dumps(payload["context_box"]),
        }
        for column in AUTO_COLUMNS:
            private_row[column] = row[column]
        private_rows.append(private_row)

    public = pd.DataFrame(public_rows)
    private_frame = pd.DataFrame(private_rows)
    public.to_csv(blinded / "annotator_1_tasks.csv", index=False, encoding="utf-8-sig")
    private_frame.to_csv(private / "trait_validation_private_key.csv", index=False, encoding="utf-8-sig")

    # Second annotator is an independent deterministic subset of the finalized packet.
    secondary = public.copy()
    secondary["_hash"] = secondary["task_id"].map(lambda x: stable_hash(x, args.seed + 3))
    secondary = secondary.sort_values("_hash").head(args.n_double).drop(columns="_hash")
    secondary.to_csv(blinded / "annotator_2_tasks.csv", index=False, encoding="utf-8-sig")

    pd.DataFrame(audit_rows).to_csv(private / "download_reconstruction_audit.csv", index=False)
    summary = {
        "n_tasks": int(len(public)),
        "n_double_labelled": int(len(secondary)),
        "n_unique_photos": int(private_frame["photo_id"].nunique()),
        "n_taxa": int(private_frame["taxon_name"].nunique()),
        "n_spatial_blocks": int(private_frame["spatial_block"].nunique()),
        "selection_rule": "one deterministic head per photo; taxon x 10-degree block spread by within-stratum rank; no production trait value/status used for selection",
        "blinding": "public packet omits taxonomy, coordinates, detector confidence, production statuses and automated trait values",
        "modules": {
            "orientation": "human long-axis endpoints plus camera-roll/assessability",
            "colour": "human visible-corolla polygon on tight crop",
            "outline": "human whole-capitulum silhouette polygon on tight crop",
        },
        "interpretation": "independent image-measurement validity; not biological calibration beyond the photographed visible phenotype",
    }
    (private / "trait_validation_packet_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (blinded / "README.txt").write_text(
        "Independent blinded Chapter 1 continuous-trait validation packet.\n"
        "Do not combine this packet with the private key before annotations are finalized.\n"
        "Run: python analysis/trait_validation_annotation_app.py --tasks annotator_1_tasks.csv --packet-root . --output annotator_1.csv\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

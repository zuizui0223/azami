#!/usr/bin/env python3
"""Build a self-contained blinded trait-annotation packet for detected heads.

One packet row corresponds to one detector box, not one plant or one photograph.
The packet intentionally excludes taxon, locality, detector confidence, and box
coordinates. A separate key preserves those values for later linkage after
annotation, but must not be used during blinded scoring.
"""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


YOLO_REQUIRED = {
    "queue_id", "audit_id", "photo_id", "source_file_stem", "screen_download_filename",
    "det_index", "crop_path", "source_image", "yolo_conf", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2",
}
CONTEXT_REQUIRED = {"queue_id", "det_index", "context_crop_status", "context_crop_path", "resolved_source_image"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a local, blinded Chapter 1 trait annotation packet.")
    parser.add_argument("--yolo-crops", required=True, help="yolo_crop_metadata.csv from one immutable inference run")
    parser.add_argument("--context-metadata", required=True, help="head_context_crop_metadata.csv")
    parser.add_argument("--source-images-dir", required=True)
    parser.add_argument("--head-crops-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--batch-name", required=True)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def det_number(value: Any) -> int:
    try:
        return int(float(text(value))) + 1
    except ValueError as error:
        raise ValueError(f"Invalid det_index: {value!r}") from error


def copy_once(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists() or target.stat().st_size != source.stat().st_size:
        shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    source_images_dir = Path(args.source_images_dir)
    head_crops_dir = Path(args.head_crops_dir)
    if not source_images_dir.is_dir():
        raise FileNotFoundError(f"--source-images-dir is not a directory: {source_images_dir}")
    if not head_crops_dir.is_dir():
        raise FileNotFoundError(f"--head-crops-dir is not a directory: {head_crops_dir}")

    yolo = pd.read_csv(args.yolo_crops, dtype=str, keep_default_na=False)
    context = pd.read_csv(args.context_metadata, dtype=str, keep_default_na=False)
    if missing := YOLO_REQUIRED.difference(yolo.columns):
        raise ValueError(f"YOLO crop metadata missing columns: {sorted(missing)}")
    if missing := CONTEXT_REQUIRED.difference(context.columns):
        raise ValueError(f"Context metadata missing columns: {sorted(missing)}")
    if yolo.duplicated(["queue_id", "det_index"]).any():
        raise ValueError("YOLO crop metadata has duplicate queue_id/det_index combinations")
    context = context.loc[context["context_crop_status"].map(text).eq("success")].copy()
    if context.duplicated(["queue_id", "det_index"]).any():
        raise ValueError("Successful context metadata has duplicate queue_id/det_index combinations")

    merged = yolo.merge(
        context[["queue_id", "det_index", "context_crop_path", "resolved_source_image"]],
        on=["queue_id", "det_index"],
        how="left",
        validate="one_to_one",
    )
    if merged["context_crop_path"].map(text).eq("").any():
        missing_rows = merged.loc[merged["context_crop_path"].map(text).eq(""), ["queue_id", "det_index"]].head(10).to_dict("records")
        raise ValueError(f"Missing successful context crop for detections: {missing_rows}")

    out_dir = Path(args.out_dir)
    packet_dir = out_dir / "trait_annotation_packet"
    private_dir = out_dir / "private_key"
    source_dir = packet_dir / "source_images"
    head_dir = packet_dir / "head_crops"
    context_dir = packet_dir / "context_crops"
    packet_rows: list[dict[str, str]] = []
    key_rows: list[dict[str, str]] = []
    missing_files: list[dict[str, str]] = []

    for _, row in merged.sort_values(["audit_id", "det_index"], key=lambda series: series.map(text)).iterrows():
        source_name = text(row["screen_download_filename"])
        head_name = Path(text(row["crop_path"])).name
        context_name = Path(text(row["context_crop_path"])).name
        source = source_images_dir / source_name
        head = head_crops_dir / head_name
        context_crop = Path(text(row["context_crop_path"]))
        if not source.is_file() or not head.is_file() or not context_crop.is_file():
            missing_files.append({
                "queue_id": text(row["queue_id"]),
                "det_index": text(row["det_index"]),
                "source_exists": str(source.is_file()),
                "head_exists": str(head.is_file()),
                "context_exists": str(context_crop.is_file()),
            })
            continue
        copy_once(source, source_dir / source_name)
        copy_once(head, head_dir / head_name)
        copy_once(context_crop, context_dir / context_name)
        number = det_number(row["det_index"])
        unit_id = f"{text(row['audit_id'])}_head_{number:02d}"
        packet_rows.append({
            "annotation_unit_id": unit_id,
            "annotation_batch": args.batch_name,
            "source_image": f"source_images/{source_name}",
            "crop_path": f"head_crops/{head_name}",
            "context_crop_path": f"context_crops/{context_name}",
        })
        key_rows.append({
            "annotation_unit_id": unit_id,
            "annotation_batch": args.batch_name,
            "queue_id": text(row["queue_id"]),
            "audit_id": text(row["audit_id"]),
            "photo_id": text(row["photo_id"]),
            "source_file_stem": text(row["source_file_stem"]),
            "det_index": text(row["det_index"]),
            "yolo_conf": text(row["yolo_conf"]),
            "bbox_x1": text(row["bbox_x1"]),
            "bbox_y1": text(row["bbox_y1"]),
            "bbox_x2": text(row["bbox_x2"]),
            "bbox_y2": text(row["bbox_y2"]),
            "source_image_original": text(row["resolved_source_image"]),
            "detector_head_crop_original": text(row["crop_path"]),
            "context_crop_original": text(row["context_crop_path"]),
        })

    if missing_files:
        raise FileNotFoundError(f"Packet build cannot continue; missing source/crop files: {missing_files[:10]}")
    if not packet_rows:
        raise ValueError("No detected-head annotation units were packaged")
    packet = pd.DataFrame(packet_rows)
    if packet["annotation_unit_id"].duplicated().any():
        raise ValueError("Generated annotation_unit_id values are not unique")
    packet_dir.mkdir(parents=True, exist_ok=True)
    private_dir.mkdir(parents=True, exist_ok=True)
    packet.to_csv(packet_dir / "blinded_trait_annotation_manifest.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(key_rows).to_csv(private_dir / "trait_annotation_key_private.csv", index=False, encoding="utf-8-sig")
    (packet_dir / "README.txt").write_text(
        "This packet contains one row per detected visible-capitulum candidate.\n"
        "It is blinded: no taxon, locality, detector confidence, or bbox coordinates are included.\n"
        "Run the ontology app from any directory; relative image paths resolve from this packet folder.\n"
        "Do not interpret the number of rows as a plant-level flower count.\n",
        encoding="utf-8",
    )
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "batch_name": args.batch_name,
        "n_detector_boxes_input": int(len(yolo)),
        "n_successful_context_crops": int(len(context)),
        "n_packaged_annotation_units": int(len(packet_rows)),
        "n_unique_source_images": int(packet["source_image"].nunique()),
        "annotation_unit": "one detector box / visible capitulum candidate in one photograph",
        "blinding": "Packet excludes taxon, locality, detector confidence, and bbox coordinates. These remain only in a separate key file.",
        "important_limit": "Detector boxes and crops are predictions. The packet does not itself validate detector correctness or define plant-level flower number.",
    }
    (out_dir / "trait_packet_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

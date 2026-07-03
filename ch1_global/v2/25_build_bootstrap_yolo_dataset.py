#!/usr/bin/env python3
"""Build a provisional single-class YOLO dataset from open-vocabulary proposals.

This is a bootstrap dataset only. Proposal boxes are explicitly retained as
pseudo-labels and must never be reported as human annotations or used for a
final accuracy estimate. The independent 1,000-image audit remains untouched.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a provisional visible-capitulum YOLO bootstrap dataset.")
    parser.add_argument("--proposals", required=True)
    parser.add_argument("--image-summary", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-score", type=float, default=0.35)
    parser.add_argument("--max-boxes-per-image", type=int, default=8)
    parser.add_argument("--validation-fraction", type=float, default=0.20)
    parser.add_argument("--seed", default="20260703")
    return parser.parse_args()


def stable_fraction(value: str, seed: str) -> float:
    digest = hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()
    return int(digest[:12], 16) / float(16**12)


def yolo_line(row: pd.Series) -> str:
    width, height = float(row["image_width"]), float(row["image_height"])
    x1, y1, x2, y2 = [float(row[column]) for column in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")]
    x_center = ((x1 + x2) / 2) / width
    y_center = ((y1 + y2) / 2) / height
    box_width = (x2 - x1) / width
    box_height = (y2 - y1) / height
    if not (0 < box_width <= 1 and 0 < box_height <= 1):
        raise ValueError("Invalid normalized proposal box")
    return f"0 {x_center:.8f} {y_center:.8f} {box_width:.8f} {box_height:.8f}"


def main() -> None:
    args = parse_args()
    if not 0 <= args.min_score <= 1 or args.max_boxes_per_image < 1 or not 0 < args.validation_fraction < 1:
        raise ValueError("Invalid bootstrap-dataset thresholds")
    proposals = pd.read_csv(args.proposals, dtype=str, keep_default_na=False)
    summary = pd.read_csv(args.image_summary, dtype=str, keep_default_na=False)
    required = {"queue_id", "proposal_score", "source_image", "image_width", "image_height", "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2"}
    missing = required.difference(proposals.columns)
    if missing:
        raise ValueError(f"Proposal table missing columns: {sorted(missing)}")
    work = proposals.copy()
    work["proposal_score"] = pd.to_numeric(work["proposal_score"], errors="coerce")
    work = work.loc[work["proposal_score"].ge(args.min_score)].sort_values(["queue_id", "proposal_score"], ascending=[True, False], kind="stable")
    work = work.groupby("queue_id", sort=False).head(args.max_boxes_per_image).copy()
    usable_ids = set(work["queue_id"])
    source_paths = summary.set_index("queue_id")["source_image"].to_dict()
    out_dir = Path(args.out_dir)
    images_root, labels_root = out_dir / "images", out_dir / "labels"
    records = []
    for queue_id, part in work.groupby("queue_id", sort=True):
        source = Path(str(source_paths.get(queue_id, "")))
        if not source.is_file():
            continue
        split = "val" if stable_fraction(str(queue_id), args.seed) < args.validation_fraction else "train"
        image_dest = images_root / split / f"{queue_id}{source.suffix.lower() or '.jpg'}"
        label_dest = labels_root / split / f"{queue_id}.txt"
        image_dest.parent.mkdir(parents=True, exist_ok=True)
        label_dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, image_dest)
        label_dest.write_text("\n".join(yolo_line(row) for _, row in part.iterrows()) + "\n", encoding="utf-8")
        records.append({"queue_id": queue_id, "split": split, "source_image": str(source), "dataset_image": str(image_dest), "n_pseudo_boxes": int(len(part)), "mean_proposal_score": float(part["proposal_score"].mean())})
    if not records:
        raise ValueError("No usable pseudo-labeled source images met the threshold")
    dataset_yaml = out_dir / "visible_capitulum_bootstrap.yaml"
    dataset_yaml.write_text(
        "path: " + str(out_dir.resolve()) + "\n"
        "train: images/train\n"
        "val: images/val\n"
        "names:\n  0: visible_capitulum\n",
        encoding="utf-8",
    )
    pd.DataFrame(records).to_csv(out_dir / "bootstrap_dataset_manifest.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_open_vocab_proposal_rows": int(len(proposals)),
        "n_pseudo_boxes_after_score_threshold": int(len(work)),
        "n_pseudo_labeled_images": int(len(records)),
        "n_train_images": int(sum(row["split"] == "train" for row in records)),
        "n_val_images": int(sum(row["split"] == "val" for row in records)),
        "min_score": args.min_score,
        "max_boxes_per_image": args.max_boxes_per_image,
        "semantic_status": "provisional pseudo-label dataset; not human-labeled ground truth",
    }
    (out_dir / "bootstrap_dataset_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

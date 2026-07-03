#!/usr/bin/env python3
"""Generate *provisional* visible-capitulum boxes using an open-vocabulary detector.

These are proposal labels for bootstrapping a new YOLO detector, not ground truth.
No accuracy claim may be made from this file. The model is queried with a narrow
botanical phrase and all boxes are retained with scores and original image IDs so
later human/adjudicated audit labels can supersede them.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from PIL import Image

PROMPT = "a thistle flower head ."


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate open-vocabulary Cirsium capitulum proposals.")
    parser.add_argument("--queue", required=True)
    parser.add_argument("--downloads", required=True)
    parser.add_argument("--images-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--model-id", default="IDEA-Research/grounding-dino-base")
    parser.add_argument("--box-threshold", type=float, default=0.20)
    parser.add_argument("--text-threshold", type=float, default=0.15)
    parser.add_argument("--max-images", type=int, default=0, help="0 means all successfully downloaded queue rows.")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def resolve_path(raw: str, images_dir: Path) -> Path:
    source = Path(text(raw))
    return source if source.exists() else images_dir / source.name


def main() -> None:
    args = parse_args()
    if not 0 <= args.box_threshold <= 1 or not 0 <= args.text_threshold <= 1:
        raise ValueError("thresholds must be in [0,1]")
    images_dir = Path(args.images_dir)
    if not images_dir.is_dir():
        raise FileNotFoundError(images_dir)
    try:
        import torch
        from transformers import AutoModelForZeroShotObjectDetection, AutoProcessor, __version__ as transformers_version
    except ImportError as error:
        raise SystemExit("Install torch and transformers before running open-vocabulary proposals.") from error

    queue = pd.read_csv(args.queue, dtype=str, keep_default_na=False)
    downloads = pd.read_csv(args.downloads, dtype=str, keep_default_na=False)
    latest = downloads.assign(_order=range(len(downloads))).sort_values("_order").groupby("queue_id", as_index=False).tail(1)
    candidates = queue.merge(latest.loc[latest["download_status"].eq("success"), ["queue_id", "image_local_path"]], on="queue_id", how="inner", validate="one_to_one")
    candidates["resolved_image_path"] = candidates["image_local_path"].map(lambda value: str(resolve_path(value, images_dir)))
    candidates = candidates.loc[candidates["resolved_image_path"].map(lambda value: Path(value).is_file())].copy()
    if args.max_images:
        candidates = candidates.head(args.max_images).copy()
    if candidates.empty:
        raise ValueError("No successfully downloaded images available for proposals.")

    processor = AutoProcessor.from_pretrained(args.model_id)
    model = AutoModelForZeroShotObjectDetection.from_pretrained(args.model_id)
    model.eval()
    proposal_rows: list[dict[str, Any]] = []
    image_rows: list[dict[str, Any]] = []
    for i, (_, row) in enumerate(candidates.iterrows(), start=1):
        image_path = Path(row["resolved_image_path"])
        try:
            with Image.open(image_path) as opened:
                image = opened.convert("RGB")
            inputs = processor(images=image, text=PROMPT, return_tensors="pt")
            with torch.no_grad():
                outputs = model(**inputs)
            target_sizes = torch.tensor([[image.height, image.width]])
            result = processor.post_process_grounded_object_detection(
                outputs,
                inputs.input_ids,
                box_threshold=args.box_threshold,
                text_threshold=args.text_threshold,
                target_sizes=target_sizes,
            )[0]
            boxes = result.get("boxes", [])
            scores = result.get("scores", [])
            labels = result.get("labels", [])
            n = 0
            for det_index, (box, score, label) in enumerate(zip(boxes, scores, labels)):
                x1, y1, x2, y2 = [float(value) for value in box.tolist()]
                proposal_rows.append({
                    "queue_id": row["queue_id"], "photo_id": row["photo_id"], "source_file_stem": row["source_file_stem"],
                    "proposal_index": det_index, "proposal_label": text(label), "proposal_score": float(score),
                    "bbox_x1": x1, "bbox_y1": y1, "bbox_x2": x2, "bbox_y2": y2,
                    "image_width": image.width, "image_height": image.height,
                    "source_image": str(image_path.resolve()),
                    "proposal_status": "unreviewed_open_vocab",
                })
                n += 1
            image_rows.append({"queue_id": row["queue_id"], "photo_id": row["photo_id"], "source_image": str(image_path.resolve()), "n_open_vocab_proposals": n, "proposal_error": ""})
        except Exception as error:
            image_rows.append({"queue_id": row["queue_id"], "photo_id": row["photo_id"], "source_image": str(image_path), "n_open_vocab_proposals": 0, "proposal_error": str(error)})
        print(f"[INFO] proposed {i}/{len(candidates)}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(proposal_rows, columns=[
        "queue_id", "photo_id", "source_file_stem", "proposal_index", "proposal_label", "proposal_score",
        "bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2", "image_width", "image_height", "source_image", "proposal_status",
    ]).to_csv(out_dir / "open_vocab_capitulum_proposals.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(image_rows).to_csv(out_dir / "open_vocab_proposal_image_summary.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model_id": args.model_id,
        "prompt": PROMPT,
        "box_threshold": args.box_threshold,
        "text_threshold": args.text_threshold,
        "n_images_attempted": int(len(candidates)),
        "n_images_with_errors": int(sum(bool(row["proposal_error"]) for row in image_rows)),
        "n_proposal_boxes": int(len(proposal_rows)),
        "semantic_status": "provisional pseudo-label proposals; not human ground truth and not an accuracy-evaluation set",
    }
    (out_dir / "open_vocab_proposal_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

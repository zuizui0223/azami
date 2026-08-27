#!/usr/bin/env python3
"""Harmonize EAzami direct-image colour proxies with the production Azami detector/crop stage.

This script receives licensed direct images that have already been screened by the
frozen visible-capitulum detector. It applies the same pinned ``colour_measurement``
function to every detector crop, then aggregates without pseudo-replication:

head -> source image -> locality/date series -> paper concept.

Multiple detected heads and multiple photos within one observation never count as
independent populations. Promotion rules are fixed before execution:

- JPN_30 promotes when >=2 independent locality/date series contain >=1 image with
  >=1 colour-usable detector crop.
- JPN_05 promotes under the same >=2-locality rule.

The colour confidence floor remains 0.55. This is a source-harmonization sensitivity,
not a discrete colour-state or evolutionary-transition analysis.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import statistics
from pathlib import Path

import cv2

COLOUR_CONFIDENCE_FLOOR = 0.55
TARGETS = {"JPN_05": "Cirsium aomorense", "JPN_30": "Cirsium yezoense"}


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_measurement(path: Path):
    spec = importlib.util.spec_from_file_location("azami_pinned_measurement", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load pinned measurement module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "colour_measurement"):
        raise RuntimeError("pinned module lacks colour_measurement")
    return module


def finite(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def median(values):
    vals = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    return statistics.median(vals) if vals else None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sources", type=Path, required=True)
    p.add_argument("--screening", type=Path, required=True)
    p.add_argument("--crops", type=Path, required=True)
    p.add_argument("--head-crops-dir", type=Path, required=True)
    p.add_argument("--measurement-script", type=Path, required=True)
    p.add_argument("--per-head-output", type=Path, required=True)
    p.add_argument("--per-image-output", type=Path, required=True)
    p.add_argument("--summary-output", type=Path, required=True)
    args = p.parse_args()

    sources = read_csv(args.sources)
    if len(sources) != 13:
        raise ValueError(f"expected 13 frozen direct-image rows, found {len(sources)}")
    by_id = {r["source_id"]: r for r in sources}
    if len(by_id) != len(sources):
        raise ValueError("duplicate source_id in direct-image registry")
    if {r["paper_japan_member_id"] for r in sources} != set(TARGETS):
        raise ValueError("unexpected paper concept in source registry")

    screening = {r["queue_id"]: r for r in read_csv(args.screening)}
    if set(screening) != set(by_id):
        raise ValueError(
            f"screening coverage mismatch missing={sorted(set(by_id)-set(screening))} "
            f"extra={sorted(set(screening)-set(by_id))}"
        )

    crop_rows = read_csv(args.crops) if args.crops.exists() and args.crops.stat().st_size else []
    azami = load_measurement(args.measurement_script)
    per_head = []
    usable_by_image: dict[str, list[dict]] = {sid: [] for sid in by_id}

    for crop in crop_rows:
        sid = crop["queue_id"]
        if sid not in by_id:
            raise ValueError(f"crop references unknown source {sid}")
        crop_path = args.head_crops_dir / Path(crop["crop_path"]).name
        image = cv2.imread(str(crop_path), cv2.IMREAD_COLOR)
        if image is None:
            raise ValueError(f"unreadable detector crop: {crop_path}")
        result = azami.colour_measurement(image)
        confidence = finite(result.get("confidence")) or 0.0
        lightness = finite(result.get("median_lab_lightness"))
        chroma = finite(result.get("median_lab_chroma"))
        usable = bool(
            str(result.get("state", "")) != getattr(azami, "NON_SCOREABLE", "unassessable")
            and confidence >= COLOUR_CONFIDENCE_FLOOR
            and lightness is not None
        )
        source = by_id[sid]
        row = {
            "source_id": sid,
            "paper_japan_member_id": source["paper_japan_member_id"],
            "taxon_name": source["taxon_name"],
            "locality_id": source["locality_id"],
            "observation_date": source["observation_date"],
            "source_type": source["source_type"],
            "det_index": crop.get("det_index", ""),
            "yolo_conf": finite(crop.get("yolo_conf")),
            "measurement_state": result.get("state"),
            "measurement_confidence": confidence,
            "median_lab_lightness": lightness,
            "median_lab_chroma": chroma,
            "floral_pixel_fraction": finite(result.get("floral_pixel_fraction")),
            "mask_quality": finite(result.get("mask_quality")),
            "usable_at_0_55": usable,
        }
        per_head.append(row)
        if usable:
            usable_by_image[sid].append(row)

    per_image = []
    for source in sources:
        sid = source["source_id"]
        screen = screening[sid]
        heads = usable_by_image[sid]
        per_image.append({
            "source_id": sid,
            "paper_japan_member_id": source["paper_japan_member_id"],
            "taxon_name": source["taxon_name"],
            "locality_id": source["locality_id"],
            "observation_date": source["observation_date"],
            "source_type": source["source_type"],
            "screen_status": screen["screen_status"],
            "n_detector_crops": int(screen["n_detections"] or 0),
            "n_colour_usable_crops": len(heads),
            "image_lightness_median": median([r["median_lab_lightness"] for r in heads]),
            "image_chroma_median": median([r["median_lab_chroma"] for r in heads]),
            "image_pass": bool(heads),
        })

    locality_summaries = []
    concept_results = {}
    for mid, taxon in TARGETS.items():
        target_rows = [r for r in per_image if r["paper_japan_member_id"] == mid]
        localities = sorted({r["locality_id"] for r in target_rows})
        passing_lightness = []
        passing_chroma = []
        passing_ids = []
        for locality in localities:
            rows = [r for r in target_rows if r["locality_id"] == locality]
            usable = [r for r in rows if r["image_pass"]]
            lmed = median([r["image_lightness_median"] for r in usable])
            cmed = median([r["image_chroma_median"] for r in usable])
            passed = bool(usable and lmed is not None)
            locality_summaries.append({
                "paper_japan_member_id": mid,
                "taxon_name": taxon,
                "locality_id": locality,
                "registered_images": len(rows),
                "passing_images": len(usable),
                "locality_lightness_median": lmed,
                "locality_chroma_median": cmed,
                "locality_pass": passed,
            })
            if passed:
                passing_ids.append(locality)
                passing_lightness.append(lmed)
                if cmed is not None:
                    passing_chroma.append(cmed)
        promotion = len(passing_ids) >= 2
        concept_results[mid] = {
            "taxon_name": taxon,
            "registered_localities": len(localities),
            "passing_localities": passing_ids,
            "passing_locality_count": len(passing_ids),
            "promotion_rule": ">=2 independent locality/date series with >=1 image containing >=1 detector crop that passes the frozen colour confidence floor 0.55",
            "promotion_pass": promotion,
            "concept_lightness_median_across_locality_medians": median(passing_lightness) if promotion else None,
            "concept_chroma_median_across_locality_medians": median(passing_chroma) if promotion else None,
        }

    result = {
        "contract_version": "eazami_direct_colour_crop_harmonization_v1",
        "status_date": "2026-08-27",
        "detector_contract": {
            "artifact_id": 8076736948,
            "best_weight_sha256": "4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed",
            "production_confidence_threshold": 0.25,
            "crop_pad_ratio": 0.12,
        },
        "measurement_contract": {
            "azami_commit": "4e6f66675004207d331fd620d6bcd894287f316a",
            "colour_confidence_floor": COLOUR_CONFIDENCE_FLOOR,
            "aggregation": "usable detector crops -> image median -> locality median -> concept median; no head/image pseudo-replication",
            "thresholds_changed_after_result": False,
        },
        "registered_direct_images": len(sources),
        "detected_images": sum(r["screen_status"] == "detected" for r in per_image),
        "total_detector_crops": len(crop_rows),
        "total_colour_usable_crops": sum(r["usable_at_0_55"] for r in per_head),
        "concept_results": concept_results,
        "locality_summaries": locality_summaries,
        "claim_boundary": "Detector-crop source harmonization only. A passing concept is a continuous Japan-local colour proxy for sensitivity analysis. This does not define W/C states, transition direction, ancestry, convergence, adaptation, evolutionary rate or molecular reactivation."
    }

    write_csv(args.per_head_output, per_head)
    write_csv(args.per_image_output, per_image)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

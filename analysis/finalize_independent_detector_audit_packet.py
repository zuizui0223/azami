#!/usr/bin/env python3
"""Freeze an exact downloaded detector-audit packet from a larger reserve selection."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-dir", required=True)
    parser.add_argument("--download-results", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-images", type=int, default=1000)
    parser.add_argument("--double-label-fraction", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=20260722)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def score(seed: int, value: Any) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}|{text(value)}".encode()).digest()[:8], "big")


def main() -> None:
    args = parse_args()
    if args.n_images < 1 or not 0 <= args.double_label_fraction <= 1:
        raise ValueError("invalid target size or double-label fraction")
    source = Path(args.selection_dir)
    queue = pd.read_csv(source / "detector_independent_audit_queue.csv", dtype=str, keep_default_na=False)
    blinded = pd.read_csv(source / "detector_independent_audit_blinded_manifest.csv", dtype=str, keep_default_na=False)
    key = pd.read_csv(source / "detector_independent_audit_private_key.csv", dtype=str, keep_default_na=False)
    downloads = pd.read_csv(args.download_results, dtype=str, keep_default_na=False)
    downloads["_order"] = range(len(downloads))
    latest = downloads.sort_values("_order").groupby("queue_id", as_index=False).tail(1)
    successful = set(latest.loc[latest["download_status"].eq("success"), "queue_id"].map(text))
    candidates = queue.loc[queue["audit_id"].map(text).isin(successful)].copy()
    if len(candidates) < args.n_images:
        raise SystemExit(f"Only {len(candidates)} reserve images downloaded successfully; require {args.n_images}")
    final_queue = candidates.head(args.n_images).copy()
    final_ids = list(final_queue["audit_id"].map(text))
    final_id_set = set(final_ids)
    final_blinded = blinded.loc[blinded["audit_id"].map(text).isin(final_id_set)].copy()
    final_key = key.loc[key["audit_id"].map(text).isin(final_id_set)].copy()
    if len(final_queue) != args.n_images or len(final_blinded) != args.n_images or len(final_key) != args.n_images:
        raise AssertionError("filtered packet tables are inconsistent")

    final_blinded["double_label"] = False
    n_double = int(round(args.n_images * args.double_label_fraction))
    selected_double = sorted(final_ids, key=lambda value: score(args.seed + 1, value))[:n_double]
    final_blinded.loc[final_blinded["audit_id"].isin(selected_double), "double_label"] = True
    assignments: list[dict[str, Any]] = []
    for _, row in final_blinded.sort_values("audit_id").iterrows():
        is_double = bool(row["double_label"])
        assignments.append({"audit_id": row["audit_id"], "annotator_id": "annotator_1", "double_label": is_double, "status": "not_started"})
        if is_double:
            assignments.append({"audit_id": row["audit_id"], "annotator_id": "annotator_2", "double_label": True, "status": "not_started"})

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    final_queue.to_csv(out / "detector_independent_audit_queue.csv", index=False)
    final_blinded.to_csv(out / "detector_independent_audit_blinded_manifest.csv", index=False)
    final_key.to_csv(out / "detector_independent_audit_private_key.csv", index=False)
    pd.DataFrame(assignments).to_csv(out / "detector_independent_audit_assignments.csv", index=False)
    pd.DataFrame(columns=[
        "audit_id", "annotator_id", "annotation_status", "assessability",
        "image_quality", "image_width", "image_height", "gt_index", "gt_label",
        "life_stage", "occlusion", "edge_truncated", "gt_x1", "gt_y1", "gt_x2",
        "gt_y2", "notes",
    ]).to_csv(out / "detector_independent_audit_annotation_template.csv", index=False)
    species_qc = final_key.groupby("taxon_name").agg(
        n_selected=("audit_id", "size"),
        n_spatial_blocks=("spatial_block", "nunique"),
        n_latitude_bands=("latitude_band", "nunique"),
    ).reset_index()
    species_qc.to_csv(out / "detector_independent_audit_species_qc.csv", index=False)
    failed = latest.loc[~latest["download_status"].eq("success")].copy()
    unused_success = candidates.iloc[args.n_images:].copy()
    report = {
        "reserve_selected_images": int(len(queue)),
        "reserve_download_success": int(len(successful)),
        "reserve_download_failures": int(len(failed)),
        "final_audit_images": int(len(final_queue)),
        "final_species": int(final_key["taxon_name"].nunique()),
        "final_spatial_blocks": int(final_key["spatial_block"].nunique()),
        "final_double_label_images": int(final_blinded["double_label"].astype(bool).sum()),
        "unused_successful_reserve_images": int(len(unused_success)),
        "selection_rule": "first successful rows in the pre-detector taxon/spatial round-robin order",
        "photo_overlap_with_prior_source": 0,
        "observation_overlap_with_prior_source": 0,
    }
    (out / "detector_independent_audit_report.json").write_text(json.dumps(report, indent=2) + "\n")
    failed.to_csv(out / "detector_independent_audit_download_failures.csv", index=False)
    unused_success.to_csv(out / "detector_independent_audit_unused_reserve.csv", index=False)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

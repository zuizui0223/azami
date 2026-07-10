#!/usr/bin/env python3
"""Assemble the candidate-units table that 47_build_representative_trait_holdout needs.

47 needs one analysable row per head crop and trait, carrying taxon and 10-degree
spatial block so the holdout can be stratified representatively. Those fields are
already produced upstream but in two different tables; this adapter is the small
glue that joins them:

- head-level AI traits (36 output `ai_trait_head_level_with_metadata.csv`):
  annotation_unit_id, trait_id, taxon_name, obs_id, analysis_state_ai_all,
  analysis_state_ai_conservative;
- the environment table (42 output): obs_id -> spatial_block_10deg.

It drops unassessable AI states and rows without a spatial block, then writes
`candidate_units.csv` with exactly the columns 47 requires. The packet manifest
(annotation_unit_id -> image paths) is a separate existing artifact from the head
packet build (28); pass an optional --packet-source here to normalise it into
`packet_manifest.csv` at the same time.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

HEAD_DEFAULTS = {"unit": "annotation_unit_id", "trait": "trait_id", "taxon": "taxon_name", "obs": "obs_id"}
PACKET_IMAGE_COLUMNS = ("source_image", "crop_path", "context_crop_path")
UNASSESSABLE = {"", "unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build candidate_units.csv (and optionally packet_manifest.csv) for 47.")
    parser.add_argument("--head-traits", required=True, help="36 output ai_trait_head_level_with_metadata.csv")
    parser.add_argument("--environment", required=True, help="42 output table carrying obs_id -> spatial block")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative")
    parser.add_argument("--state-column", default=None, help="Override the AI state column (default analysis_state_ai_<mode>)")
    parser.add_argument("--unit-column", default=HEAD_DEFAULTS["unit"])
    parser.add_argument("--trait-column", default=HEAD_DEFAULTS["trait"])
    parser.add_argument("--taxon-column", default=HEAD_DEFAULTS["taxon"])
    parser.add_argument("--obs-column", default=HEAD_DEFAULTS["obs"])
    parser.add_argument("--env-obs-column", default="obs_id")
    parser.add_argument("--block-column", default="spatial_block_10deg")
    parser.add_argument("--traits", nargs="*", default=None, help="Restrict to these trait_ids")
    parser.add_argument("--packet-source", default=None,
                        help="Optional table with annotation_unit_id + image paths to normalise into packet_manifest.csv")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def is_assessable(value: Any) -> bool:
    return text(value) not in UNASSESSABLE


def require(frame: pd.DataFrame, columns: set[str], label: str) -> None:
    if missing := columns.difference(frame.columns):
        raise ValueError(f"{label} missing columns: {sorted(missing)}")


def write_packet_manifest(packet_source: str, unit_column: str, out: Path) -> int:
    source = pd.read_csv(packet_source, dtype=str, keep_default_na=False)
    rename = {unit_column: "annotation_unit_id"} if unit_column != "annotation_unit_id" else {}
    source = source.rename(columns=rename)
    require(source, {"annotation_unit_id", *PACKET_IMAGE_COLUMNS}, "packet source")
    manifest = source[["annotation_unit_id", *PACKET_IMAGE_COLUMNS]].drop_duplicates("annotation_unit_id")
    manifest.to_csv(out / "packet_manifest.csv", index=False, encoding="utf-8-sig")
    return int(len(manifest))


def main() -> None:
    args = parse_args()
    state_column = args.state_column or f"analysis_state_ai_{args.measurement_mode}"

    head = pd.read_csv(args.head_traits, dtype=str, keep_default_na=False)
    env = pd.read_csv(args.environment, dtype=str, keep_default_na=False)
    require(head, {args.unit_column, args.trait_column, args.taxon_column, args.obs_column, state_column}, "head-traits table")
    require(env, {args.env_obs_column, args.block_column}, "environment table")

    head = head.rename(columns={
        args.unit_column: "annotation_unit_id", args.trait_column: "trait_id",
        args.taxon_column: "taxon_name", args.obs_column: "obs_id", state_column: "ai_candidate_state",
    })
    if head.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("head-traits table must have unique annotation_unit_id/trait_id rows")
    if args.traits:
        head = head.loc[head["trait_id"].isin(args.traits)].copy()

    block = (
        env.rename(columns={args.env_obs_column: "obs_id", args.block_column: "spatial_block_10deg"})
        [["obs_id", "spatial_block_10deg"]]
        .drop_duplicates("obs_id")
    )
    n_input = int(len(head))
    head = head.loc[head["ai_candidate_state"].map(is_assessable)].copy()
    n_after_assessable = int(len(head))

    merged = head.merge(block, on="obs_id", how="left", validate="many_to_one")
    has_block = merged["spatial_block_10deg"].map(text).ne("")
    n_no_block = int((~has_block).sum())
    candidate = merged.loc[has_block, [
        "annotation_unit_id", "trait_id", "taxon_name", "spatial_block_10deg", "ai_candidate_state",
    ]].copy()
    if candidate.empty:
        raise ValueError("No candidate units after dropping unassessable states and rows without a spatial block")

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    candidate.to_csv(out / "candidate_units.csv", index=False, encoding="utf-8-sig")
    n_manifest = write_packet_manifest(args.packet_source, args.unit_column, out) if args.packet_source else None

    per_trait = {
        str(trait): {
            "n_candidate_units": int(len(group)),
            "n_taxa": int(group["taxon_name"].nunique()),
            "n_spatial_blocks": int(group["spatial_block_10deg"].nunique()),
        }
        for trait, group in candidate.groupby("trait_id", sort=True)
    }
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "state_column": state_column,
        "n_head_trait_rows_input": n_input,
        "n_after_dropping_unassessable": n_after_assessable,
        "n_dropped_missing_spatial_block": n_no_block,
        "n_candidate_units": int(len(candidate)),
        "n_traits": int(candidate["trait_id"].nunique()),
        "packet_manifest_rows": n_manifest,
        "per_trait": per_trait,
        "next_step": "Feed candidate_units.csv (and packet_manifest.csv) to 47_build_representative_trait_holdout.py.",
    }
    (out / "candidate_units_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

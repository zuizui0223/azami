#!/usr/bin/env python3
"""Attach wide head-plus-peduncle context crops to an existing annotation manifest.

The legacy manifest builder remains usable. This bridge joins context crops by
source image and detection index (or queue ID where available) without exposing
taxonomy in the annotation packet.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Attach context_crop_path to a Chapter 1 annotation manifest.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--context-metadata", required=True, help="head_context_crop_metadata.csv")
    parser.add_argument("--out-path", required=True)
    return parser.parse_args()


def text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    context = pd.read_csv(args.context_metadata, dtype=str, keep_default_na=False)
    required_manifest = {"annotation_unit_id", "source_image", "det_index"}
    required_context = {"source_image", "det_index", "context_crop_path", "context_crop_status"}
    missing = required_manifest.difference(manifest.columns)
    if missing:
        raise ValueError(f"Manifest missing columns: {sorted(missing)}")
    missing = required_context.difference(context.columns)
    if missing:
        raise ValueError(f"Context metadata missing columns: {sorted(missing)}")
    context = context.loc[context["context_crop_status"].map(text).eq("success")].copy()
    context["source_image"] = context["source_image"].map(text)
    context["det_index"] = context["det_index"].map(text)
    if context.duplicated(["source_image", "det_index"]).any():
        raise ValueError("Context metadata has duplicate source_image/det_index combinations")
    output = manifest.copy()
    output["source_image"] = output["source_image"].map(text)
    output["det_index"] = output["det_index"].map(text)
    output = output.merge(
        context[["source_image", "det_index", "context_crop_path"]],
        on=["source_image", "det_index"], how="left", validate="one_to_one",
    )
    output["context_crop_path"] = output["context_crop_path"].fillna("")
    output.to_csv(Path(args.out_path), index=False, encoding="utf-8-sig")
    print(f"[OK] attached context crop paths for {(output['context_crop_path'] != '').sum()} / {len(output)} manifest rows")


if __name__ == "__main__":
    main()

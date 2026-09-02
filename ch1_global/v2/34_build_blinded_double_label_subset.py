#!/usr/bin/env python3
"""Build a public-only subset for predesignated double annotation.

The private key is used only to choose task IDs. The subset intentionally omits
model candidates, margins, taxonomy, locality, selection strata, and the fact
that a task was selected for double annotation.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


PUBLIC_REQUIRED = {"task_id", "annotation_unit_id", "trait_id", "source_image", "crop_path", "context_crop_path"}
PRIVATE_REQUIRED = {"task_id", "trait_id", "double_label"}
PUBLIC_COLUMNS = ["task_id", "annotation_unit_id", "trait_id", "source_image", "crop_path", "context_crop_path"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a blinded second-annotator task subset.")
    parser.add_argument("--public-packet-root", required=True)
    parser.add_argument("--private-audit-key", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--batch-name", required=True)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def copy_once(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists() or target.stat().st_size != source.stat().st_size:
        shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    public_root = Path(args.public_packet_root)
    task_path = public_root / "blinded_trait_audit_tasks.csv"
    key_path = Path(args.private_audit_key)
    if not task_path.is_file():
        raise FileNotFoundError(task_path)
    if not key_path.is_file():
        raise FileNotFoundError(key_path)

    public = pd.read_csv(task_path, dtype=str, keep_default_na=False)
    private = pd.read_csv(key_path, dtype=str, keep_default_na=False)
    if missing := PUBLIC_REQUIRED.difference(public.columns):
        raise ValueError(f"Public audit tasks missing columns: {sorted(missing)}")
    if missing := PRIVATE_REQUIRED.difference(private.columns):
        raise ValueError(f"Private audit key missing columns: {sorted(missing)}")
    if public.empty or public["task_id"].duplicated().any():
        raise ValueError("Public audit task IDs must be nonempty and unique")
    if private["task_id"].duplicated().any():
        raise ValueError("Private audit key task IDs must be unique")

    double_key = private.loc[private["double_label"].map(as_bool), ["task_id", "trait_id"]].copy()
    if double_key.empty:
        raise ValueError("Private audit key has no prespecified double-label tasks")
    selected = public.merge(double_key, on=["task_id", "trait_id"], how="inner", validate="one_to_one")
    if len(selected) != len(double_key):
        unmatched = sorted(set(double_key["task_id"]).difference(selected["task_id"]))
        raise ValueError(f"Designated double-label task IDs are absent from public packet: {unmatched[:10]}")
    selected = selected[PUBLIC_COLUMNS].sort_values("task_id").reset_index(drop=True)

    output = Path(args.out_dir)
    packet_root = output / "blinded_second_annotator_packet"
    packet_root.mkdir(parents=True, exist_ok=True)
    for row in selected.to_dict("records"):
        for field in ("source_image", "crop_path", "context_crop_path"):
            relative = Path(text(row[field]))
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Unsafe packet-relative path in {field}: {relative}")
            source = public_root / relative
            if not source.is_file():
                raise FileNotFoundError(f"Public packet image missing: {source}")
            copy_once(source, packet_root / relative)

    selected.to_csv(packet_root / "blinded_trait_audit_tasks.csv", index=False, encoding="utf-8-sig")
    (packet_root / "README.txt").write_text(
        "Cirsium Ch.1 — independent focused trait audit subset\n\n"
        "This is a blinded subset for an independent annotator. It contains no AI candidate labels, model confidence/margins, taxon names, locality, selection strata, or allocation flags.\n\n"
        "Score only the trait named in each task. Use unassessable whenever the available view cannot support a decision.\n"
        "Use an annotator ID different from the primary annotator and keep responses in a separate responses/<annotator_id>/ folder.\n",
        encoding="utf-8",
    )
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "batch_name": args.batch_name,
        "n_tasks": int(len(selected)),
        "n_unique_annotation_units": int(selected["annotation_unit_id"].nunique()),
        "task_counts_by_trait": selected["trait_id"].value_counts().sort_index().to_dict(),
        "public_columns_only": PUBLIC_COLUMNS,
        "blinding": "No model candidate, margin, taxon, locality, selection stratum, or allocation flag is included in this packet.",
    }
    (output / "blinded_second_annotator_subset_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

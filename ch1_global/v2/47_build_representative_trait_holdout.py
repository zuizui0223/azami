#!/usr/bin/env python3
"""Build a representative, species- and space-stratified blinded trait holdout.

The existing audit chain (30-33) deliberately samples *near the decision margin*
(low-margin uncertainty + high-margin calibration). Its own evaluator states it
"never interprets this non-random audit sample as global population accuracy" and
that "a separate independent, representative holdout is required before any
automatic acceptance policy." This script builds exactly that missing holdout: a
representative random sample, stratified by taxon and 10-degree spatial block, so
that agreement measured on it estimates population-level measurement validity
(and, via the per-species breakdown in 48, species transferability / LOSO).

It selects tasks only; reuse the existing packet-image plumbing (30) and audit
app (31) for the human step. The public task list hides the AI candidate, taxon,
locality, and stratum; those stay in the private key for post-review evaluation
by 48_evaluate_representative_trait_holdout.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

CANDIDATE_REQUIRED = {"annotation_unit_id", "trait_id", "taxon_name", "spatial_block_10deg", "ai_candidate_state"}
UNASSESSABLE = {"", "unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a representative stratified blinded trait holdout selection.")
    parser.add_argument("--candidate-units", required=True, help="CSV with one row per analysable annotation-unit/trait")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--batch-name", required=True)
    parser.add_argument("--n-tasks-per-trait", type=int, default=120)
    parser.add_argument("--double-label-fraction", type=float, default=0.25)
    parser.add_argument("--min-block-taxa", type=int, default=1, help="Minimum candidates required in a taxon x block cell to be sampled")
    parser.add_argument("--seed", default="20260710")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def stable_rank(value: str, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def is_assessable(value: Any) -> bool:
    return text(value) not in UNASSESSABLE


def representative_sample(part: pd.DataFrame, n_target: int, seed: str) -> pd.DataFrame:
    """Round-robin across (taxon, block) strata so no single species/region dominates.

    Within each stratum, order is a stable per-unit hash (deterministic, not tied to
    the AI margin). Round-robin draw approximates proportional-yet-spread coverage.
    """
    if part.empty or n_target <= 0:
        return part.head(0)
    part = part.copy()
    part["_rank"] = part["annotation_unit_id"].map(lambda v: stable_rank(text(v), seed))
    strata = [
        stratum.sort_values("_rank")
        for _, stratum in sorted(part.groupby(["taxon_name", "spatial_block_10deg"], sort=True), key=lambda kv: str(kv[0]))
    ]
    picked: list[pd.Series] = []
    cursor = 0
    while len(picked) < min(n_target, len(part)):
        advanced = False
        for stratum in strata:
            if cursor < len(stratum):
                picked.append(stratum.iloc[cursor])
                advanced = True
                if len(picked) >= min(n_target, len(part)):
                    break
        if not advanced:
            break
        cursor += 1
    return pd.DataFrame(picked).drop(columns="_rank", errors="ignore")


def main() -> None:
    args = parse_args()
    if args.n_tasks_per_trait < 1:
        raise ValueError("--n-tasks-per-trait must be positive")
    if not 0 <= args.double_label_fraction <= 1:
        raise ValueError("--double-label-fraction must be in [0,1]")

    candidates = pd.read_csv(args.candidate_units, dtype=str, keep_default_na=False)
    if missing := CANDIDATE_REQUIRED.difference(candidates.columns):
        raise ValueError(f"Candidate units missing columns: {sorted(missing)}")
    if candidates.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Candidate units must have one row per annotation-unit/trait")

    analysable = candidates.loc[candidates["ai_candidate_state"].map(is_assessable)].copy()
    if analysable.empty:
        raise ValueError("No analysable candidate units after removing unassessable AI states")

    # Drop sparse taxon x block cells so the representative frame is not padded with singletons
    cell_counts = analysable.groupby(["trait_id", "taxon_name", "spatial_block_10deg"])["annotation_unit_id"].transform("size")
    analysable = analysable.loc[cell_counts.ge(args.min_block_taxa)].copy()

    selected_rows: list[pd.DataFrame] = []
    per_trait_report: dict[str, dict[str, int]] = {}
    for trait, part in analysable.groupby("trait_id", sort=True):
        chosen = representative_sample(part, args.n_tasks_per_trait, args.seed)
        chosen = chosen.assign(selection_stratum="representative_random")
        selected_rows.append(chosen)
        per_trait_report[str(trait)] = {
            "n_candidate_units": int(len(part)),
            "n_selected": int(len(chosen)),
            "n_taxa_selected": int(chosen["taxon_name"].nunique()),
            "n_spatial_blocks_selected": int(chosen["spatial_block_10deg"].nunique()),
            "shortfall": int(max(0, args.n_tasks_per_trait - len(chosen))),
        }
    if not selected_rows:
        raise ValueError("No tasks selected")
    selected = pd.concat(selected_rows, ignore_index=True)
    selected["task_id"] = [f"trait_holdout_{index:05d}" for index in range(1, len(selected) + 1)]
    selected["_double_rank"] = selected["task_id"].map(lambda v: stable_rank(v, args.seed))
    n_double = int(round(len(selected) * args.double_label_fraction))
    double_ids = set(selected.sort_values("_double_rank").head(n_double)["task_id"])
    selected["double_label"] = selected["task_id"].isin(double_ids)

    public = selected[["task_id", "annotation_unit_id", "trait_id"]].copy()
    private = selected[[
        "task_id", "annotation_unit_id", "trait_id", "taxon_name",
        "spatial_block_10deg", "selection_stratum", "double_label", "ai_candidate_state",
    ]].copy()
    private["double_label"] = private["double_label"].map(lambda flag: str(bool(flag)).lower())

    out = Path(args.out_dir)
    public_root = out / "representative_trait_holdout_tasks"
    private_root = out / "representative_trait_holdout_key"
    public_root.mkdir(parents=True, exist_ok=True)
    private_root.mkdir(parents=True, exist_ok=True)
    public.to_csv(public_root / "representative_trait_holdout_tasks.csv", index=False, encoding="utf-8-sig")
    (public_root / "README.txt").write_text(
        "One row is one representative trait audit task drawn by taxon x spatial-block strata.\n"
        "AI candidate, taxon, locality and stratum are hidden. Score visible evidence only;\n"
        "choose unassessable when the required view is absent. Reuse the 30/31 packet+app\n"
        "plumbing to attach and display images for these annotation_unit_id values.\n",
        encoding="utf-8",
    )
    private.to_csv(private_root / "representative_trait_holdout_key.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "batch_name": args.batch_name,
        "n_tasks": int(len(selected)),
        "n_double_label_tasks": int(selected["double_label"].sum()),
        "n_traits": int(selected["trait_id"].nunique()),
        "per_trait": per_trait_report,
        "selection_design": "Representative random sample stratified by taxon x 10-degree spatial block. "
        "Selection does not use the model margin, so downstream agreement estimates population-level "
        "measurement validity (48 reports the per-species breakdown for species transferability).",
    }
    (out / "representative_trait_holdout_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

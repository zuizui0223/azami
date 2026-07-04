#!/usr/bin/env python3
"""Build a compact, blinded trait-audit packet from zero-shot CLIP proposals.

The public packet contains only task identity, target trait, and images. It does
not expose CLIP candidates, uncertainty margins, taxon, locality, or selection
stratum. Those remain in a separate private key for post-review evaluation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a compact blinded Cirsium trait-audit packet.")
    parser.add_argument("--packet-manifest", required=True)
    parser.add_argument("--packet-root", required=True)
    parser.add_argument("--proposals-long", required=True)
    parser.add_argument("--prompt-spec", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--batch-name", required=True)
    parser.add_argument("--n-uncertainty-tasks", type=int, default=180)
    parser.add_argument("--n-calibration-tasks", type=int, default=70)
    parser.add_argument("--double-label-fraction", type=float, default=0.25)
    parser.add_argument("--seed", default="20260704")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def stable_hash(value: str, seed: str) -> str:
    return hashlib.sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def weighted_quotas(total: int, traits: list[str], weights: dict[str, float]) -> dict[str, int]:
    if total < 0 or not traits:
        raise ValueError("Invalid task quota request")
    if total == 0:
        return {trait: 0 for trait in traits}
    raw = {trait: total * weights[trait] / sum(weights.values()) for trait in traits}
    quota = {trait: int(math.floor(raw[trait])) for trait in traits}
    remaining = total - sum(quota.values())
    for trait in sorted(traits, key=lambda trait: (raw[trait] - quota[trait], trait), reverse=True)[:remaining]:
        quota[trait] += 1
    return quota


def copy_once(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists() or target.stat().st_size != source.stat().st_size:
        shutil.copy2(source, target)


def main() -> None:
    args = parse_args()
    if args.n_uncertainty_tasks < 1 or args.n_calibration_tasks < 1:
        raise ValueError("Task counts must be positive")
    if not 0 <= args.double_label_fraction <= 1:
        raise ValueError("--double-label-fraction must be in [0,1]")
    packet_root = Path(args.packet_root)
    packet_manifest = pd.read_csv(args.packet_manifest, dtype=str, keep_default_na=False)
    proposals = pd.read_csv(args.proposals_long, dtype=str, keep_default_na=False)
    prompt_spec = json.loads(Path(args.prompt_spec).read_text(encoding="utf-8"))
    if not packet_root.is_dir() or packet_manifest.empty:
        raise ValueError("Trait packet root or manifest is unavailable")
    required_packet = {"annotation_unit_id", "source_image", "crop_path", "context_crop_path"}
    required_proposals = {
        "annotation_unit_id", "trait_id", "proposal_status", "ai_candidate_state",
        "runner_up_state", "similarity_margin", "margin_percentile_within_trait",
    }
    if missing := required_packet.difference(packet_manifest.columns):
        raise ValueError(f"Packet manifest missing columns: {sorted(missing)}")
    if missing := required_proposals.difference(proposals.columns):
        raise ValueError(f"Proposal table missing columns: {sorted(missing)}")
    if packet_manifest["annotation_unit_id"].duplicated().any():
        raise ValueError("Trait packet annotation_unit_id values must be unique")

    candidates = proposals.loc[proposals["proposal_status"].eq("zero_shot_candidate_not_validated")].copy()
    candidates["margin_percentile_within_trait"] = pd.to_numeric(candidates["margin_percentile_within_trait"], errors="raise")
    candidates["similarity_margin"] = pd.to_numeric(candidates["similarity_margin"], errors="raise")
    if candidates.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Proposal candidates must have one row per annotation-unit/trait")
    traits = sorted(candidates["trait_id"].unique().tolist())
    config = prompt_spec.get("traits", {})
    weights: dict[str, float] = {}
    for trait in traits:
        value = config.get(trait, {}).get("review_weight", 1.0)
        try:
            weights[trait] = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError(f"Invalid review weight for {trait!r}") from error
        if weights[trait] <= 0:
            raise ValueError(f"Review weight must be positive for {trait!r}")

    uncertain_quota = weighted_quotas(args.n_uncertainty_tasks, traits, weights)
    calibration_quota = weighted_quotas(args.n_calibration_tasks, traits, {trait: 1.0 for trait in traits})
    selected_rows: list[pd.DataFrame] = []
    for trait in traits:
        part = candidates.loc[candidates["trait_id"].eq(trait)].copy()
        uncertain = part.sort_values(["margin_percentile_within_trait", "annotation_unit_id"], ascending=[True, True]).head(uncertain_quota[trait]).copy()
        uncertain["selection_stratum"] = "low_margin_uncertainty"
        selected_keys = set(zip(uncertain["annotation_unit_id"], uncertain["trait_id"]))
        remaining = part.loc[~part.apply(lambda row: (row["annotation_unit_id"], row["trait_id"]) in selected_keys, axis=1)].copy()
        calibration = remaining.sort_values(["margin_percentile_within_trait", "annotation_unit_id"], ascending=[False, True]).head(calibration_quota[trait]).copy()
        calibration["selection_stratum"] = "high_margin_calibration"
        selected_rows.extend([uncertain, calibration])
    selected = pd.concat(selected_rows, ignore_index=True)
    required_total = args.n_uncertainty_tasks + args.n_calibration_tasks
    if len(selected) != required_total:
        raise ValueError(f"Could not fill requested task count: got={len(selected)}, expected={required_total}")
    selected["task_id"] = [f"trait_audit_{index:04d}" for index in range(1, len(selected) + 1)]
    selected["_double_rank"] = selected["task_id"].map(lambda value: stable_hash(value, args.seed))
    n_double = int(round(len(selected) * args.double_label_fraction))
    double_ids = set(selected.sort_values("_double_rank").head(n_double)["task_id"])
    selected["double_label"] = selected["task_id"].isin(double_ids)

    source_by_unit = packet_manifest.set_index("annotation_unit_id")
    public_rows: list[dict[str, str]] = []
    private_rows: list[dict[str, str]] = []
    output = Path(args.out_dir)
    public_root = output / "blinded_trait_audit_packet"
    for _, row in selected.iterrows():
        unit_id = text(row["annotation_unit_id"])
        packet_row = source_by_unit.loc[unit_id]
        image_values = {key: text(packet_row[key]) for key in ("source_image", "crop_path", "context_crop_path")}
        for relative in image_values.values():
            source = packet_root / relative
            if not source.is_file():
                raise FileNotFoundError(f"Packet image missing: {source}")
            copy_once(source, public_root / relative)
        public_rows.append({
            "task_id": text(row["task_id"]),
            "annotation_unit_id": unit_id,
            "trait_id": text(row["trait_id"]),
            **image_values,
        })
        private_rows.append({
            "task_id": text(row["task_id"]),
            "annotation_unit_id": unit_id,
            "trait_id": text(row["trait_id"]),
            "selection_stratum": text(row["selection_stratum"]),
            "double_label": str(bool(row["double_label"])).lower(),
            "ai_candidate_state": text(row["ai_candidate_state"]),
            "runner_up_state": text(row["runner_up_state"]),
            "similarity_margin": text(row["similarity_margin"]),
            "margin_percentile_within_trait": text(row["margin_percentile_within_trait"]),
            "model_id": text(row.get("model_id", "")),
            "prompt_spec_version": text(row.get("prompt_spec_version", "")),
        })

    public = pd.DataFrame(public_rows)
    private = pd.DataFrame(private_rows)
    if public["task_id"].duplicated().any() or len(public) != required_total:
        raise ValueError("Public task manifest is not uniquely complete")
    public_root.mkdir(parents=True, exist_ok=True)
    public.to_csv(public_root / "blinded_trait_audit_tasks.csv", index=False, encoding="utf-8-sig")
    (public_root / "README.txt").write_text(
        "One row is one trait-specific image audit task.\n"
        "This packet intentionally hides AI candidates, confidence/margins, taxon, locality, and selection stratum.\n"
        "Use the trait-specific audit app and score visible evidence only; choose unassessable when the required view is absent.\n",
        encoding="utf-8",
    )
    private_root = output / "private_audit_key"
    private_root.mkdir(parents=True, exist_ok=True)
    private.to_csv(private_root / "trait_audit_selection_key.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "batch_name": args.batch_name,
        "n_tasks": int(len(public)),
        "n_unique_annotation_units": int(public["annotation_unit_id"].nunique()),
        "n_uncertainty_tasks": int((private["selection_stratum"] == "low_margin_uncertainty").sum()),
        "n_calibration_tasks": int((private["selection_stratum"] == "high_margin_calibration").sum()),
        "n_double_label_tasks": int((private["double_label"] == "true").sum()),
        "task_counts_by_trait": private["trait_id"].value_counts().sort_index().to_dict(),
        "uncertainty_counts_by_trait": private.loc[private["selection_stratum"].eq("low_margin_uncertainty"), "trait_id"].value_counts().sort_index().to_dict(),
        "calibration_counts_by_trait": private.loc[private["selection_stratum"].eq("high_margin_calibration"), "trait_id"].value_counts().sort_index().to_dict(),
        "selection_design": "Low-margin tasks prioritize uncertainty review; high-margin tasks provide a calibration component. AI candidates remain outside the public packet.",
    }
    (output / "trait_audit_packet_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

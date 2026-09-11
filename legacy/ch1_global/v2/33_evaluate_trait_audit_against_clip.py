#!/usr/bin/env python3
"""Evaluate zero-shot CLIP trait candidates against compiled blinded audit evidence.

This reports agreement only within the deliberately selected audit tasks. It
never promotes an automatic trait label to truth, and it never interprets this
non-random audit sample as global population accuracy.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate CLIP trait candidates against blinded human audit responses.")
    parser.add_argument("--canonical-annotations", required=True)
    parser.add_argument("--private-audit-key", required=True)
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def wilson_lower(successes: int, total: int, z: float = 1.96) -> float | None:
    if total <= 0:
        return None
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = proportion + z * z / (2 * total)
    spread = z * math.sqrt((proportion * (1 - proportion) + z * z / (4 * total)) / total)
    return (centre - spread) / denominator


def load_ontology(path: Path) -> dict[str, dict[str, Any]]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"trait_id", "allow_multiple"}
    if missing := required.difference(table.columns):
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    return {
        text(row["trait_id"]): {"allow_multiple": as_bool(row["allow_multiple"])}
        for row in table.to_dict("records")
    }


def is_assessable(value: Any) -> bool:
    return text(value) not in {"", "unassessable"}


def candidate_matches_human(candidate: Any, human_state: Any, allow_multiple: bool) -> bool | None:
    candidate_text, human_text = text(candidate), text(human_state)
    if not candidate_text or not is_assessable(human_text):
        return None
    return candidate_text in set(human_text.split("|")) if allow_multiple else candidate_text == human_text


def summarize_group(frame: pd.DataFrame) -> dict[str, object]:
    assessable = frame.loc[frame["human_assessable"]].copy()
    successes = int(assessable["model_candidate_match"].sum()) if not assessable.empty else 0
    total = int(len(assessable))
    return {
        "n_human_records": int(len(frame)),
        "n_assessable": total,
        "n_unassessable": int(len(frame) - total),
        "n_model_candidate_matches": successes,
        "candidate_match_rate": successes / total if total else None,
        "wilson_lower_95": wilson_lower(successes, total),
    }


def pairwise_agreement(evaluated: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for task_id, group in evaluated.groupby("task_id", sort=False):
        rows = group.sort_values("annotator_id").to_dict("records")
        for left, right in itertools.combinations(rows, 2):
            if left["annotator_id"] == right["annotator_id"]:
                continue
            records.append({
                "task_id": task_id,
                "trait_id": left["trait_id"],
                "double_label": as_bool(left["double_label"]),
                "annotator_a": left["annotator_id"],
                "annotator_b": right["annotator_id"],
                "human_state_a": left["human_state"],
                "human_state_b": right["human_state"],
                "both_assessable": is_assessable(left["human_state"]) and is_assessable(right["human_state"]),
                "exact_agreement": text(left["human_state"]) == text(right["human_state"]),
            })
    return pd.DataFrame(records)


def main() -> None:
    args = parse_args()
    annotations = pd.read_csv(args.canonical_annotations, dtype=str, keep_default_na=False)
    key = pd.read_csv(args.private_audit_key, dtype=str, keep_default_na=False)
    ontology = load_ontology(Path(args.ontology))
    required_annotations = {"task_id", "annotation_unit_id", "trait_id", "human_state", "annotator_id"}
    required_key = {"task_id", "trait_id", "selection_stratum", "double_label", "ai_candidate_state", "runner_up_state", "similarity_margin", "margin_percentile_within_trait"}
    if missing := required_annotations.difference(annotations.columns):
        raise ValueError(f"Canonical annotation file missing columns: {sorted(missing)}")
    if missing := required_key.difference(key.columns):
        raise ValueError(f"Private audit key missing columns: {sorted(missing)}")
    if annotations.empty:
        raise ValueError("No canonical annotations to evaluate")
    if key["task_id"].duplicated().any():
        raise ValueError("Private audit key task IDs must be unique")
    if annotations.duplicated(["task_id", "annotator_id"]).any():
        raise ValueError("Canonical annotations duplicate task_id/annotator_id")
    unknown_traits = sorted(set(annotations["trait_id"]).difference(ontology))
    if unknown_traits:
        raise ValueError(f"Annotations contain traits absent from ontology: {unknown_traits}")

    evaluated = annotations.merge(key, on=["task_id", "trait_id"], how="left", validate="many_to_one", suffixes=("", "_key"))
    if evaluated["selection_stratum"].map(text).eq("").any():
        missing_task_ids = evaluated.loc[evaluated["selection_stratum"].map(text).eq(""), "task_id"].tolist()
        raise ValueError(f"Annotations missing private-key matches: {missing_task_ids[:10]}")
    evaluated["human_assessable"] = evaluated["human_state"].map(is_assessable)
    evaluated["allow_multiple"] = evaluated["trait_id"].map(lambda trait: ontology[trait]["allow_multiple"])
    evaluated["model_candidate_match"] = [
        candidate_matches_human(candidate, human, allow_multiple)
        for candidate, human, allow_multiple in zip(evaluated["ai_candidate_state"], evaluated["human_state"], evaluated["allow_multiple"])
    ]
    evaluated["model_candidate_match"] = evaluated["model_candidate_match"].fillna(False).astype(bool)
    evaluated["similarity_margin"] = pd.to_numeric(evaluated["similarity_margin"], errors="raise")
    evaluated["margin_percentile_within_trait"] = pd.to_numeric(evaluated["margin_percentile_within_trait"], errors="raise")

    summary = pd.DataFrame([
        {"trait_id": trait_id, "selection_stratum": stratum, **summarize_group(group)}
        for (trait_id, stratum), group in evaluated.groupby(["trait_id", "selection_stratum"], sort=True)
    ])
    pairs = pairwise_agreement(evaluated)
    agreement_rows: list[dict[str, object]] = []
    if not pairs.empty:
        for trait_id, group in pairs.groupby("trait_id", sort=True):
            designated = group.loc[group["double_label"]].copy()
            assessable = designated.loc[designated["both_assessable"]]
            successes, total = int(assessable["exact_agreement"].sum()) if not assessable.empty else 0, int(len(assessable))
            agreement_rows.append({
                "trait_id": trait_id,
                "n_pairwise_records_all_tasks": int(len(group)),
                "n_pairwise_records_designated_double_label": int(len(designated)),
                "n_both_assessable_pairs_designated_double_label": total,
                "n_exact_agreements_designated_double_label": successes,
                "exact_agreement_rate_designated_double_label": successes / total if total else None,
                "wilson_lower_95_designated_double_label": wilson_lower(successes, total),
            })
    agreement = pd.DataFrame(agreement_rows)

    policy_rows: list[dict[str, object]] = []
    for trait_id in sorted(set(key["trait_id"])):
        calibration = evaluated.loc[(evaluated["trait_id"] == trait_id) & (evaluated["selection_stratum"] == "high_margin_calibration")]
        calibration_stats = summarize_group(calibration)
        agreement_stats = agreement.loc[agreement["trait_id"].eq(trait_id)]
        pair_count = int(agreement_stats.iloc[0]["n_both_assessable_pairs_designated_double_label"]) if not agreement_stats.empty else 0
        pair_lower = agreement_stats.iloc[0]["wilson_lower_95_designated_double_label"] if not agreement_stats.empty else None
        n_calibration, lower = int(calibration_stats["n_assessable"]), calibration_stats["wilson_lower_95"]
        if n_calibration < 15:
            recommendation = "manual_review_required_insufficient_high_margin_calibration"
        elif lower is None or lower < 0.80:
            recommendation = "manual_review_required_low_model_human_agreement"
        elif pair_count < 10:
            recommendation = "candidate_assist_only_need_more_designated_double_annotation"
        elif pair_lower is None or pair_lower < 0.70:
            recommendation = "candidate_assist_only_low_designated_human_agreement"
        else:
            recommendation = "candidate_assist_only_pending_independent_holdout"
        policy_rows.append({
            "trait_id": trait_id,
            "n_high_margin_assessable": n_calibration,
            "high_margin_candidate_match_rate": calibration_stats["candidate_match_rate"],
            "high_margin_wilson_lower_95": lower,
            "n_designated_double_annotated_assessable_pairs": pair_count,
            "human_human_wilson_lower_95_designated_double_label": pair_lower,
            "recommendation": recommendation,
            "hard_limit": "No automatic trait acceptance is authorised by this selected first-pass audit alone.",
        })
    policy = pd.DataFrame(policy_rows)

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    evaluated.to_csv(output / "model_vs_human_trait_audit_records.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(output / "model_vs_human_summary_by_trait_and_stratum.csv", index=False, encoding="utf-8-sig")
    pairs.to_csv(output / "human_human_pairwise_records.csv", index=False, encoding="utf-8-sig")
    agreement.to_csv(output / "human_human_agreement_by_trait.csv", index=False, encoding="utf-8-sig")
    policy.to_csv(output / "trait_candidate_use_policy.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_human_annotation_records": int(len(annotations)),
        "n_unique_audited_tasks": int(annotations["task_id"].nunique()),
        "n_annotators": int(annotations["annotator_id"].nunique()),
        "n_pairwise_human_records_all_tasks": int(len(pairs)),
        "n_pairwise_human_records_designated_double_label": int(pairs["double_label"].sum()) if not pairs.empty else 0,
        "selection_scope": "Results are conditional on the deliberately selected low-margin uncertainty and high-margin calibration tasks. They are not estimates of global trait-classifier accuracy.",
        "acceptance_boundary": "This report never authorises unreviewed model predictions as final trait observations. A separate independent, representative holdout is required before any automatic acceptance policy.",
    }
    (output / "trait_audit_evaluation_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

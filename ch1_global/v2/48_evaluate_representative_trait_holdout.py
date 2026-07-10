#!/usr/bin/env python3
"""Certify head-level AI traits with one representative human sample.

The production gate uses overall AI-human agreement and agreement between the
two human annotators. Per-species accuracy is still reported descriptively. A
species-specific gate can be enabled explicitly for legacy experiments, but it
is off by default because the original atlas already contains held-out-species
validation.
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

ANNOTATION_REQUIRED = {"task_id", "annotation_unit_id", "trait_id", "human_state", "annotator_id"}
KEY_REQUIRED = {"task_id", "annotation_unit_id", "trait_id", "taxon_name", "double_label"}
UNITS_REQUIRED = {"annotation_unit_id", "trait_id", "taxon_name"}
UNASSESSABLE = {"", "unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Certify AI traits on a representative human holdout")
    parser.add_argument("--canonical-annotations", required=True)
    parser.add_argument("--holdout-key", required=True)
    parser.add_argument("--analysis-units", required=True)
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative")
    parser.add_argument("--primary-annotator", default="")
    parser.add_argument("--min-assessable-per-trait", type=int, default=50)
    parser.add_argument("--accuracy-wilson-floor", type=float, default=0.80)
    parser.add_argument("--human-human-wilson-floor", type=float, default=0.80)
    # Optional legacy species gate. Off in the production design.
    parser.add_argument("--min-species-per-trait", type=int, default=0)
    parser.add_argument("--min-records-per-species", type=int, default=5)
    parser.add_argument("--worst-species-accuracy-floor", type=float, default=0.60)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def is_assessable(value: Any) -> bool:
    return text(value) not in UNASSESSABLE


def wilson_lower(successes: int, total: int, z: float = 1.96) -> float | None:
    if total <= 0:
        return None
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = proportion + z * z / (2 * total)
    spread = z * math.sqrt((proportion * (1 - proportion) + z * z / (4 * total)) / total)
    return (centre - spread) / denominator


def load_ontology(path: Path) -> dict[str, bool]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = {"trait_id", "allow_multiple"}.difference(table.columns)
    if missing:
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    return {text(row["trait_id"]): as_bool(row["allow_multiple"]) for row in table.to_dict("records")}


def state_matches(ai_state: Any, human_state: Any, allow_multiple: bool) -> bool | None:
    ai_text, human_text = text(ai_state), text(human_state)
    if not is_assessable(ai_text) or not is_assessable(human_text):
        return None
    return ai_text in set(human_text.split("|")) if allow_multiple else ai_text == human_text


def choose_state_column(units: pd.DataFrame, mode: str) -> str:
    for column in (f"analysis_state_ai_{mode}", f"observation_ai_{mode}_state"):
        if column in units.columns:
            return column
    raise ValueError(f"Analysis units have no AI state column for mode={mode}")


def choose_primary_annotator(annotations: pd.DataFrame, requested: str) -> str:
    requested = requested.strip()
    available = set(annotations["annotator_id"].map(text))
    if requested:
        if requested not in available:
            raise ValueError(f"Primary annotator {requested!r} is absent")
        return requested
    counts = annotations.groupby("annotator_id")["task_id"].nunique().sort_values(ascending=False)
    if counts.empty:
        raise ValueError("No annotators found")
    maximum = int(counts.iloc[0])
    return sorted(counts.loc[counts.eq(maximum)].index.map(text).tolist())[0]


def human_human_pairs(annotations: pd.DataFrame, double_tasks: set[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for task_id, group in annotations.loc[annotations["task_id"].isin(double_tasks)].groupby("task_id"):
        records = group.drop_duplicates("annotator_id").to_dict("records")
        for left, right in itertools.combinations(records, 2):
            if not (is_assessable(left["human_state"]) and is_assessable(right["human_state"])):
                continue
            rows.append({
                "task_id": task_id,
                "trait_id": left["trait_id"],
                "annotator_left": text(left["annotator_id"]),
                "annotator_right": text(right["annotator_id"]),
                "state_left": text(left["human_state"]),
                "state_right": text(right["human_state"]),
                "exact_agreement": text(left["human_state"]) == text(right["human_state"]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    annotations = pd.read_csv(args.canonical_annotations, dtype=str, keep_default_na=False)
    key = pd.read_csv(args.holdout_key, dtype=str, keep_default_na=False)
    units = pd.read_csv(args.analysis_units, dtype=str, keep_default_na=False)
    ontology = load_ontology(Path(args.ontology))

    for label, frame, required in (
        ("annotations", annotations, ANNOTATION_REQUIRED),
        ("holdout key", key, KEY_REQUIRED),
        ("analysis units", units, UNITS_REQUIRED),
    ):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{label} missing columns: {sorted(missing)}")
    if annotations.empty:
        raise ValueError("No annotations")
    if annotations.duplicated(["task_id", "annotator_id"]).any():
        raise ValueError("Duplicate task/annotator rows")
    if key["task_id"].duplicated().any():
        raise ValueError("Duplicate holdout task IDs")
    if units.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Duplicate analysis-unit rows")

    unknown_traits = sorted(set(annotations["trait_id"]).difference(ontology))
    if unknown_traits:
        raise ValueError(f"Unknown traits: {unknown_traits}")

    state_column = choose_state_column(units, args.measurement_mode)
    primary_annotator = choose_primary_annotator(annotations, args.primary_annotator)
    ai = units[["annotation_unit_id", "trait_id", state_column]].rename(columns={state_column: "ai_state"})
    merged = (
        annotations
        .merge(key[["task_id", "taxon_name", "double_label"]], on="task_id", how="left", validate="many_to_one")
        .merge(ai, on=["annotation_unit_id", "trait_id"], how="left", validate="many_to_one")
    )
    if merged["taxon_name"].map(text).eq("").any() or merged["ai_state"].map(text).eq("").any():
        raise ValueError("Annotations do not match the private key or AI table")

    merged["match"] = [
        state_matches(ai_state, human_state, ontology[trait])
        for ai_state, human_state, trait in zip(merged["ai_state"], merged["human_state"], merged["trait_id"])
    ]
    primary = merged.loc[merged["annotator_id"].map(text).eq(primary_annotator)].copy()
    scoreable = primary.loc[primary["match"].notna()].copy()

    double_tasks = set(key.loc[key["double_label"].map(as_bool), "task_id"])
    human_pairs = human_human_pairs(merged, double_tasks)
    hh_by_trait = {
        trait: (int(group["exact_agreement"].sum()), int(len(group)))
        for trait, group in human_pairs.groupby("trait_id")
    } if not human_pairs.empty else {}

    species_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for trait in sorted(set(key["trait_id"])):
        group = scoreable.loc[scoreable["trait_id"].eq(trait)]
        successes = int(group["match"].sum())
        total = int(len(group))
        accuracy_lower = wilson_lower(successes, total)

        eligible_species = 0
        species_accuracies: list[float] = []
        for taxon, taxon_group in group.groupby("taxon_name"):
            n = int(len(taxon_group))
            hits = int(taxon_group["match"].sum())
            eligible = n >= args.min_records_per_species
            species_rows.append({
                "trait_id": trait,
                "taxon_name": taxon,
                "n_scoreable": n,
                "n_match": hits,
                "accuracy": hits / n if n else None,
                "eligible_for_optional_species_gate": eligible,
            })
            if eligible:
                eligible_species += 1
                species_accuracies.append(hits / n)
        worst_species = min(species_accuracies) if species_accuracies else None

        hh_hits, hh_total = hh_by_trait.get(trait, (0, 0))
        hh_lower = wilson_lower(hh_hits, hh_total)
        reasons: list[str] = []
        if total < args.min_assessable_per_trait:
            reasons.append("insufficient_sample")
        if accuracy_lower is None or accuracy_lower < args.accuracy_wilson_floor:
            reasons.append("low_ai_human_agreement")
        if hh_total == 0:
            reasons.append("no_double_label_evidence")
        elif hh_lower is None or hh_lower < args.human_human_wilson_floor:
            reasons.append("low_human_human_agreement")
        if args.min_species_per_trait > 0:
            if eligible_species < args.min_species_per_trait:
                reasons.append("too_few_species_with_evidence")
            if worst_species is None or worst_species < args.worst_species_accuracy_floor:
                reasons.append("species_generalization_gap")

        summary_rows.append({
            "trait_id": trait,
            "measurement_mode": args.measurement_mode,
            "primary_annotator": primary_annotator,
            "n_scoreable": total,
            "ai_human_accuracy": successes / total if total else None,
            "ai_human_wilson_lower_95": accuracy_lower,
            "human_human_pairs": hh_total,
            "human_human_accuracy": hh_hits / hh_total if hh_total else None,
            "human_human_wilson_lower_95": hh_lower,
            "n_species_described": int(group["taxon_name"].nunique()),
            "n_species_with_min_records": eligible_species,
            "worst_species_accuracy_descriptive": worst_species,
            "gate_decision": "accept_for_analysis" if not reasons else "manual_review_or_exclude",
            "gate_reasons": ";".join(reasons),
        })

    summary = pd.DataFrame(summary_rows)
    per_species = pd.DataFrame(species_rows)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    scoreable.to_csv(out / "representative_holdout_scored_records.csv", index=False, encoding="utf-8-sig")
    human_pairs.to_csv(out / "representative_holdout_human_human_pairs.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(out / "representative_holdout_trait_gate.csv", index=False, encoding="utf-8-sig")
    per_species.to_csv(out / "representative_holdout_per_species_accuracy.csv", index=False, encoding="utf-8-sig")

    accepted = summary.loc[summary["gate_decision"].eq("accept_for_analysis"), "trait_id"].tolist()
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "ai_state_column": state_column,
        "primary_annotator": primary_annotator,
        "n_traits_evaluated": int(len(summary)),
        "traits_accepted_for_analysis": accepted,
        "n_traits_accepted": int(len(accepted)),
        "production_gate": {
            "min_assessable_per_trait": args.min_assessable_per_trait,
            "ai_human_wilson_floor": args.accuracy_wilson_floor,
            "human_human_wilson_floor": args.human_human_wilson_floor,
        },
        "species_gate_enabled": args.min_species_per_trait > 0,
        "semantic_status": "Independent human validation of automated head-level traits. Species transfer is handled primarily by the existing held-out-species analysis.",
    }
    (out / "representative_holdout_certification.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

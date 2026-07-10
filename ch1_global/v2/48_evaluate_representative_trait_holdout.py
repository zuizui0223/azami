#!/usr/bin/env python3
"""Certify AI trait measurement validity on a representative, species-aware holdout.

This is the acceptance instrument the audit chain (30-33) explicitly defers to:
its evaluator (33) can only ever recommend "candidate_assist_only_pending_
independent_holdout" because it audits a margin-stratified, non-representative
sample. Fed the representative holdout from 47 plus human labels, this script
reports population-level agreement AND a per-species breakdown (a leave-one-
species-out diagnostic for whether the model reads the trait rather than species
identity), then applies a pre-declared pass/manual-review gate per trait.

Unlike 33, a passing trait here IS authorised for analysis use, because the
sample is representative. Failing traits are held for manual review or excluded
from the environmental atlas headline. Nothing here establishes causal adaptation.
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
    parser = argparse.ArgumentParser(description="Certify AI trait measurement on a representative holdout.")
    parser.add_argument("--canonical-annotations", required=True, help="Human labels: task_id, annotation_unit_id, trait_id, human_state, annotator_id")
    parser.add_argument("--holdout-key", required=True, help="Private key from 47 (taxon + double_label per task)")
    parser.add_argument("--analysis-units", required=True, help="AI analysis units with observation_ai_*_state columns")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative",
                        help="Which AI state column to certify (observation_ai_<mode>_state)")
    parser.add_argument("--min-assessable-per-trait", type=int, default=50)
    parser.add_argument("--accuracy-wilson-floor", type=float, default=0.80)
    parser.add_argument("--min-species-per-trait", type=int, default=10)
    parser.add_argument("--min-records-per-species", type=int, default=5)
    parser.add_argument("--worst-species-accuracy-floor", type=float, default=0.60)
    parser.add_argument("--human-human-wilson-floor", type=float, default=0.80)
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


def load_ontology(path: Path) -> dict[str, dict[str, Any]]:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    if missing := {"trait_id", "allow_multiple"}.difference(table.columns):
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    return {text(row["trait_id"]): {"allow_multiple": as_bool(row["allow_multiple"])} for row in table.to_dict("records")}


def state_matches(ai_state: Any, human_state: Any, allow_multiple: bool) -> bool | None:
    ai_text, human_text = text(ai_state), text(human_state)
    if not is_assessable(ai_text) or not is_assessable(human_text):
        return None
    return ai_text in set(human_text.split("|")) if allow_multiple else ai_text == human_text


def human_human_agreement(annotations: pd.DataFrame, double_tasks: set[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    designated = annotations.loc[annotations["task_id"].isin(double_tasks)]
    for task_id, group in designated.groupby("task_id", sort=False):
        records = group.drop_duplicates("annotator_id").to_dict("records")
        for left, right in itertools.combinations(records, 2):
            if not (is_assessable(left["human_state"]) and is_assessable(right["human_state"])):
                continue
            rows.append({
                "trait_id": left["trait_id"],
                "task_id": task_id,
                "exact_agreement": text(left["human_state"]) == text(right["human_state"]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    annotations = pd.read_csv(args.canonical_annotations, dtype=str, keep_default_na=False)
    key = pd.read_csv(args.holdout_key, dtype=str, keep_default_na=False)
    units = pd.read_csv(args.analysis_units, dtype=str, keep_default_na=False)
    ontology = load_ontology(Path(args.ontology))

    for name, frame, required in (
        ("canonical annotations", annotations, ANNOTATION_REQUIRED),
        ("holdout key", key, KEY_REQUIRED),
        ("analysis units", units, UNITS_REQUIRED),
    ):
        if missing := required.difference(frame.columns):
            raise ValueError(f"{name} missing columns: {sorted(missing)}")
    state_column = f"observation_ai_{args.measurement_mode}_state"
    if state_column not in units.columns:
        raise ValueError(f"Analysis units missing measurement column: {state_column}")
    if annotations.empty:
        raise ValueError("No canonical annotations to evaluate")
    if annotations.duplicated(["task_id", "annotator_id"]).any():
        raise ValueError("Canonical annotations duplicate task_id/annotator_id")
    if key["task_id"].duplicated().any():
        raise ValueError("Holdout key task_id values must be unique")
    unknown = sorted(set(annotations["trait_id"]).difference(ontology))
    if unknown:
        raise ValueError(f"Annotations reference traits absent from ontology: {unknown}")

    ai_state = units[["annotation_unit_id", "trait_id", state_column]].rename(columns={state_column: "ai_state"})
    merged = (
        annotations
        .merge(key[["task_id", "taxon_name", "double_label"]], on="task_id", how="left", validate="many_to_one")
        .merge(ai_state, on=["annotation_unit_id", "trait_id"], how="left", validate="many_to_one")
    )
    if merged["taxon_name"].map(text).eq("").any():
        raise ValueError("Some annotations have no matching holdout-key task_id")
    merged["allow_multiple"] = merged["trait_id"].map(lambda trait: ontology[trait]["allow_multiple"])
    merged["match"] = [
        state_matches(ai, human, allow_multiple)
        for ai, human, allow_multiple in zip(merged["ai_state"], merged["human_state"], merged["allow_multiple"])
    ]
    # Scoreable = both AI and human assessable (one canonical human row per task assumed;
    # collapse double-labelled tasks to the first annotator for the accuracy denominator).
    merged["scoreable"] = merged["match"].notna()
    scoreable = merged.loc[merged["scoreable"]].drop_duplicates(["task_id", "trait_id"]).copy()

    double_tasks = set(key.loc[key["double_label"].map(as_bool), "task_id"])
    hh = human_human_agreement(merged, double_tasks)
    hh_by_trait = {
        trait: (int(group["exact_agreement"].sum()), int(len(group)))
        for trait, group in hh.groupby("trait_id", sort=True)
    } if not hh.empty else {}

    per_species_rows: list[dict[str, object]] = []
    trait_rows: list[dict[str, object]] = []
    for trait, group in scoreable.groupby("trait_id", sort=True):
        successes, total = int(group["match"].sum()), int(len(group))
        acc_lower = wilson_lower(successes, total)
        species_acc: list[float] = []
        eligible_species = 0
        for taxon, taxon_group in group.groupby("taxon_name", sort=True):
            n = int(len(taxon_group))
            hits = int(taxon_group["match"].sum())
            per_species_rows.append({
                "trait_id": trait, "taxon_name": taxon,
                "n_scoreable": n, "n_match": hits,
                "accuracy": hits / n if n else None,
            })
            if n >= args.min_records_per_species:
                eligible_species += 1
                species_acc.append(hits / n)
        worst_species_acc = min(species_acc) if species_acc else None
        hh_hits, hh_total = hh_by_trait.get(trait, (0, 0))
        hh_lower = wilson_lower(hh_hits, hh_total)

        reasons: list[str] = []
        if total < args.min_assessable_per_trait:
            reasons.append("insufficient_representative_sample")
        if acc_lower is None or acc_lower < args.accuracy_wilson_floor:
            reasons.append("low_population_agreement")
        if eligible_species < args.min_species_per_trait:
            reasons.append("too_few_species_with_evidence")
        if worst_species_acc is not None and worst_species_acc < args.worst_species_accuracy_floor:
            reasons.append("species_generalization_gap")
        if hh_total > 0 and (hh_lower is None or hh_lower < args.human_human_wilson_floor):
            reasons.append("low_human_human_agreement")
        elif hh_total == 0:
            reasons.append("no_double_label_evidence")
        decision = "accept_for_analysis" if not reasons else "manual_review_or_exclude"
        trait_rows.append({
            "trait_id": trait,
            "measurement_mode": args.measurement_mode,
            "n_scoreable": total,
            "population_accuracy": successes / total if total else None,
            "population_accuracy_wilson_lower_95": acc_lower,
            "n_species_eligible": eligible_species,
            "worst_species_accuracy": worst_species_acc,
            "human_human_pairs": hh_total,
            "human_human_wilson_lower_95": hh_lower,
            "gate_decision": decision,
            "gate_reasons": ";".join(reasons),
        })
    trait_summary = pd.DataFrame(trait_rows)
    per_species = pd.DataFrame(per_species_rows)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    scoreable.to_csv(out / "representative_holdout_scored_records.csv", index=False, encoding="utf-8-sig")
    trait_summary.to_csv(out / "representative_holdout_trait_gate.csv", index=False, encoding="utf-8-sig")
    per_species.to_csv(out / "representative_holdout_per_species_accuracy.csv", index=False, encoding="utf-8-sig")
    accepted = trait_summary.loc[trait_summary["gate_decision"].eq("accept_for_analysis"), "trait_id"].tolist() if not trait_summary.empty else []
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "thresholds": {
            "min_assessable_per_trait": args.min_assessable_per_trait,
            "accuracy_wilson_floor": args.accuracy_wilson_floor,
            "min_species_per_trait": args.min_species_per_trait,
            "min_records_per_species": args.min_records_per_species,
            "worst_species_accuracy_floor": args.worst_species_accuracy_floor,
            "human_human_wilson_floor": args.human_human_wilson_floor,
        },
        "n_traits_evaluated": int(len(trait_summary)),
        "traits_accepted_for_analysis": accepted,
        "n_traits_accepted": int(len(accepted)),
        "semantic_status": "Representative-holdout measurement-validity certification for automated image traits. "
        "Accepted traits are authorised as analysis inputs; this establishes measurement validity and species "
        "transferability, not causal adaptation, trait evolution, or pollinator-mediated selection.",
        "acceptance_boundary": "Only traits with gate_decision=accept_for_analysis should headline the "
        "trait-environment atlas. Others are candidate/manual-review only.",
    }
    (out / "representative_holdout_certification.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

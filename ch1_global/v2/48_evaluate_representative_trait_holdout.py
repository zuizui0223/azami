#!/usr/bin/env python3
"""Certify head-level AI trait measurement on an independent human holdout.

Population accuracy is estimated only from the representative-population
stratum. A separately marked species-diagnostic top-up may contribute to the
per-species transferability check, but never to the population-accuracy
estimate. This keeps the two inferential targets distinct.
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
REPRESENTATIVE_STRATUM = "representative_population"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Certify AI trait measurement on a representative holdout.")
    parser.add_argument("--canonical-annotations", required=True)
    parser.add_argument("--holdout-key", required=True)
    parser.add_argument("--analysis-units", required=True, help="Head-level AI table; observation-level legacy fixtures are also accepted")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative")
    parser.add_argument("--primary-annotator", default="", help="Annotator used for AI-human accuracy; default is the annotator with most tasks")
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
    return {
        text(row["trait_id"]): {"allow_multiple": as_bool(row["allow_multiple"])}
        for row in table.to_dict("records")
    }


def state_matches(ai_state: Any, human_state: Any, allow_multiple: bool) -> bool | None:
    ai_text, human_text = text(ai_state), text(human_state)
    if not is_assessable(ai_text) or not is_assessable(human_text):
        return None
    return ai_text in set(human_text.split("|")) if allow_multiple else ai_text == human_text


def choose_state_column(units: pd.DataFrame, mode: str) -> str:
    candidates = [
        f"analysis_state_ai_{mode}",
        f"observation_ai_{mode}_state",
    ]
    for candidate in candidates:
        if candidate in units.columns:
            return candidate
    raise ValueError(f"Analysis units missing AI state column; expected one of {candidates}")


def choose_primary_annotator(annotations: pd.DataFrame, requested: str) -> str:
    requested = requested.strip()
    available = set(annotations["annotator_id"].map(text))
    if requested:
        if requested not in available:
            raise ValueError(f"Primary annotator {requested!r} is absent from canonical annotations")
        return requested
    counts = (
        annotations.assign(annotator_id=annotations["annotator_id"].map(text))
        .groupby("annotator_id")["task_id"]
        .nunique()
        .sort_values(ascending=False)
    )
    if counts.empty or not text(counts.index[0]):
        raise ValueError("Canonical annotations contain no annotator IDs")
    maximum = int(counts.iloc[0])
    return sorted(counts.loc[counts.eq(maximum)].index.tolist())[0]


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

    for name, frame, required in (
        ("canonical annotations", annotations, ANNOTATION_REQUIRED),
        ("holdout key", key, KEY_REQUIRED),
        ("analysis units", units, UNITS_REQUIRED),
    ):
        if missing := required.difference(frame.columns):
            raise ValueError(f"{name} missing columns: {sorted(missing)}")
    if annotations.empty:
        raise ValueError("No canonical annotations to evaluate")
    if annotations.duplicated(["task_id", "annotator_id"]).any():
        raise ValueError("Canonical annotations duplicate task_id/annotator_id")
    if key["task_id"].duplicated().any():
        raise ValueError("Holdout key task_id values must be unique")
    if units.duplicated(["annotation_unit_id", "trait_id"]).any():
        raise ValueError("Analysis units duplicate annotation_unit_id/trait_id")

    unknown = sorted(set(annotations["trait_id"]).difference(ontology))
    if unknown:
        raise ValueError(f"Annotations reference traits absent from ontology: {unknown}")

    state_column = choose_state_column(units, args.measurement_mode)
    primary_annotator = choose_primary_annotator(annotations, args.primary_annotator)
    if "selection_stratum" not in key.columns:
        key["selection_stratum"] = REPRESENTATIVE_STRATUM
    key["selection_stratum"] = key["selection_stratum"].map(text).replace("", REPRESENTATIVE_STRATUM)

    ai_state = units[["annotation_unit_id", "trait_id", state_column]].rename(columns={state_column: "ai_state"})
    key_columns = ["task_id", "taxon_name", "double_label", "selection_stratum"]
    if "ai_candidate_state" in key.columns:
        key_columns.append("ai_candidate_state")
    merged = (
        annotations
        .merge(key[key_columns], on="task_id", how="left", validate="many_to_one")
        .merge(ai_state, on=["annotation_unit_id", "trait_id"], how="left", validate="many_to_one")
    )
    if merged["taxon_name"].map(text).eq("").any():
        raise ValueError("Some annotations have no matching holdout-key task_id")
    if merged["ai_state"].map(text).eq("").any():
        raise ValueError("Some holdout tasks have no matching head-level AI state")
    if "ai_candidate_state" in merged.columns:
        mismatch = merged.loc[
            merged["ai_candidate_state"].map(text).ne("")
            & merged["ai_candidate_state"].map(text).ne(merged["ai_state"].map(text))
        ]
        if not mismatch.empty:
            raise ValueError("Holdout key AI states disagree with the supplied analysis-units table")

    merged["allow_multiple"] = merged["trait_id"].map(lambda trait: ontology[trait]["allow_multiple"])
    merged["match"] = [
        state_matches(ai, human, allow_multiple)
        for ai, human, allow_multiple in zip(merged["ai_state"], merged["human_state"], merged["allow_multiple"])
    ]
    merged["scoreable"] = merged["match"].notna()

    primary = merged.loc[merged["annotator_id"].map(text).eq(primary_annotator)].copy()
    if primary.empty:
        raise ValueError("Primary annotator has no records")
    scoreable_all = primary.loc[primary["scoreable"]].copy()
    population_scoreable = scoreable_all.loc[
        scoreable_all["selection_stratum"].eq(REPRESENTATIVE_STRATUM)
    ].copy()

    double_tasks = set(key.loc[key["double_label"].map(as_bool), "task_id"])
    hh = human_human_agreement(merged, double_tasks)
    hh_by_trait = {
        trait: (int(group["exact_agreement"].sum()), int(len(group)))
        for trait, group in hh.groupby("trait_id", sort=True)
    } if not hh.empty else {}

    per_species_rows: list[dict[str, object]] = []
    trait_rows: list[dict[str, object]] = []
    expected_traits = sorted(set(key["trait_id"]))
    for trait in expected_traits:
        population_group = population_scoreable.loc[population_scoreable["trait_id"].eq(trait)]
        diagnostic_group = scoreable_all.loc[scoreable_all["trait_id"].eq(trait)]

        successes = int(population_group["match"].sum())
        total = int(len(population_group))
        acc_lower = wilson_lower(successes, total)

        species_acc: list[float] = []
        eligible_species = 0
        for taxon, taxon_group in diagnostic_group.groupby("taxon_name", sort=True):
            n = int(len(taxon_group))
            hits = int(taxon_group["match"].sum())
            eligible = n >= args.min_records_per_species
            per_species_rows.append({
                "trait_id": trait,
                "taxon_name": taxon,
                "n_scoreable": n,
                "n_match": hits,
                "accuracy": hits / n if n else None,
                "eligible_for_worst_species_gate": eligible,
            })
            if eligible:
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
        if worst_species_acc is None or worst_species_acc < args.worst_species_accuracy_floor:
            reasons.append("species_generalization_gap")
        if hh_total == 0:
            reasons.append("no_double_label_evidence")
        elif hh_lower is None or hh_lower < args.human_human_wilson_floor:
            reasons.append("low_human_human_agreement")

        trait_rows.append({
            "trait_id": trait,
            "measurement_mode": args.measurement_mode,
            "ai_state_column": state_column,
            "primary_annotator": primary_annotator,
            "n_population_scoreable": total,
            "n_species_diagnostic_scoreable": int(len(diagnostic_group)),
            "population_accuracy": successes / total if total else None,
            "population_accuracy_wilson_lower_95": acc_lower,
            "n_species_eligible": eligible_species,
            "worst_species_accuracy": worst_species_acc,
            "human_human_pairs": hh_total,
            "human_human_accuracy": hh_hits / hh_total if hh_total else None,
            "human_human_wilson_lower_95": hh_lower,
            "gate_decision": "accept_for_analysis" if not reasons else "manual_review_or_exclude",
            "gate_reasons": ";".join(reasons),
        })

    trait_summary = pd.DataFrame(trait_rows)
    per_species = pd.DataFrame(per_species_rows)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    scoreable_all.to_csv(out / "representative_holdout_scored_records.csv", index=False, encoding="utf-8-sig")
    population_scoreable.to_csv(out / "representative_holdout_population_records.csv", index=False, encoding="utf-8-sig")
    hh.to_csv(out / "representative_holdout_human_human_pairs.csv", index=False, encoding="utf-8-sig")
    trait_summary.to_csv(out / "representative_holdout_trait_gate.csv", index=False, encoding="utf-8-sig")
    per_species.to_csv(out / "representative_holdout_per_species_accuracy.csv", index=False, encoding="utf-8-sig")

    accepted = trait_summary.loc[
        trait_summary["gate_decision"].eq("accept_for_analysis"), "trait_id"
    ].tolist()
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "ai_state_column": state_column,
        "primary_annotator": primary_annotator,
        "population_accuracy_stratum": REPRESENTATIVE_STRATUM,
        "species_diagnostic_uses_all_scoreable_strata": True,
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
        "semantic_status": "Representative-population measurement validation plus a separately marked species-transfer diagnostic. Passing authorises a trait as an input to the existing association atlas; it does not establish causal adaptation.",
        "acceptance_boundary": "Only traits with gate_decision=accept_for_analysis may headline the Chapter 1 atlas.",
    }
    (out / "representative_holdout_certification.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

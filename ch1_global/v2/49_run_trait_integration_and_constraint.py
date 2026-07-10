#!/usr/bin/env python3
"""Species-level trait coordination, integration, and constraint (Ch.1 README §3.5).

This is novelty leg 2 of Chapter 1 (global floral trait syndromes / integration /
constraint), which the executed association atlas does not yet cover. It asks not
just "do traits cluster" but "how coordinated and how constrained is capitulum
trait space." It is computed at the SPECIES level, because photo-level
correlations are inflated by within-individual pseudoreplication and shared
measurement error (README §3.5), so the observation-level long table is first
aggregated to one modal state per species per trait.

Outputs (descriptive only; no causal adaptation, evolution, or selection claim):
- species x trait modal-state matrix;
- pairwise trait association (bias-corrected Cramer's V);
- an integration-magnitude permutation test (mean pairwise association vs a null
  that shuffles each trait independently across species);
- a constraint test (number of distinct realised trait-state combinations vs the
  same independence null) plus a table of depleted / forbidden combinations.

These results become explicit hypotheses of correlated evolution for Chapter 2;
they are not evidence of adaptation on their own.
"""

from __future__ import annotations

import argparse
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

LONG_REQUIRED = {"obs_id", "trait_id", "taxon_name"}
UNASSESSABLE = {"", "unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Species-level trait integration and constraint analysis.")
    parser.add_argument("--observation-long", required=True, help="ai_trait_observation_level_long.csv")
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--measurement-mode", choices=["conservative", "all"], default="conservative")
    parser.add_argument("--traits", nargs="*", default=None, help="Restrict to these trait_ids (default: all present)")
    parser.add_argument("--min-obs-per-species-trait", type=int, default=3,
                        help="Minimum assessable observations for a species-trait modal state to be trusted")
    parser.add_argument("--min-species", type=int, default=20, help="Minimum complete-case species to run")
    parser.add_argument("--n-permutations", type=int, default=999)
    parser.add_argument("--min-expected-forbidden", type=float, default=1.0,
                        help="Flag an unobserved combination as forbidden when its independence-expected count exceeds this")
    parser.add_argument("--max-forbidden-enumeration", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260710)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def is_assessable(value: Any) -> bool:
    return text(value) not in UNASSESSABLE


def cramers_v(a: pd.Series, b: pd.Series) -> float | None:
    """Bias-corrected Cramer's V between two categorical series (aligned, no NaN)."""
    table = pd.crosstab(a, b)
    if table.shape[0] < 2 or table.shape[1] < 2:
        return None
    counts = table.to_numpy(dtype=float)
    n = counts.sum()
    if n <= 0:
        return None
    row = counts.sum(axis=1, keepdims=True)
    col = counts.sum(axis=0, keepdims=True)
    expected = row @ col / n
    chi2 = float(((counts - expected) ** 2 / expected).sum())
    phi2 = chi2 / n
    r, k = counts.shape
    phi2corr = max(0.0, phi2 - (k - 1) * (r - 1) / (n - 1))
    rcorr = r - (r - 1) ** 2 / (n - 1)
    kcorr = k - (k - 1) ** 2 / (n - 1)
    denominator = min(kcorr - 1, rcorr - 1)
    if denominator <= 0:
        return None
    return float(np.sqrt(phi2corr / denominator))


def mean_pairwise_association(matrix: pd.DataFrame, traits: list[str]) -> tuple[float, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    values: list[float] = []
    for left, right in itertools.combinations(traits, 2):
        v = cramers_v(matrix[left], matrix[right])
        rows.append({"trait_a": left, "trait_b": right, "cramers_v": v})
        if v is not None:
            values.append(v)
    mean_v = float(np.mean(values)) if values else float("nan")
    return mean_v, pd.DataFrame(rows)


def build_species_matrix(long: pd.DataFrame, state_col: str, traits: list[str], min_obs: int) -> pd.DataFrame:
    work = long.loc[long["trait_id"].isin(traits)].copy()
    work = work.loc[work[state_col].map(is_assessable)]
    modal: dict[str, dict[str, str]] = {}
    for (taxon, trait), group in work.groupby(["taxon_name", "trait_id"], sort=True):
        states = group[state_col].map(text)
        if len(states) < min_obs:
            continue
        modal.setdefault(text(taxon), {})[text(trait)] = states.value_counts().idxmax()
    frame = pd.DataFrame.from_dict(modal, orient="index").reindex(columns=traits)
    frame.index.name = "taxon_name"
    return frame.dropna(axis=0, how="any")


def permutation_null(matrix: pd.DataFrame, traits: list[str], n_perm: int, rng: np.random.Generator) -> np.ndarray:
    values = np.empty(n_perm, dtype=float)
    columns = {trait: matrix[trait].to_numpy() for trait in traits}
    for i in range(n_perm):
        shuffled = pd.DataFrame({trait: rng.permutation(columns[trait]) for trait in traits}, index=matrix.index)
        values[i], _ = mean_pairwise_association(shuffled, traits)
    return values


def distinct_combinations(matrix: pd.DataFrame, traits: list[str]) -> int:
    return int(matrix[traits].drop_duplicates().shape[0])


def combination_null(matrix: pd.DataFrame, traits: list[str], n_perm: int, rng: np.random.Generator) -> np.ndarray:
    values = np.empty(n_perm, dtype=float)
    columns = {trait: matrix[trait].to_numpy() for trait in traits}
    for i in range(n_perm):
        shuffled = pd.DataFrame({trait: rng.permutation(columns[trait]) for trait in traits}, index=matrix.index)
        values[i] = distinct_combinations(shuffled, traits)
    return values


def forbidden_combinations(matrix: pd.DataFrame, traits: list[str], min_expected: float, cap: int) -> pd.DataFrame:
    n = len(matrix)
    marginals = {trait: matrix[trait].value_counts(normalize=True).to_dict() for trait in traits}
    grid_size = int(np.prod([len(marginals[trait]) for trait in traits]))
    if grid_size > cap:
        return pd.DataFrame([{ "note": f"combination grid {grid_size} exceeds cap {cap}; skipped enumeration"}])
    observed = matrix.groupby(traits, sort=False).size().to_dict()
    rows: list[dict[str, object]] = []
    for combo in itertools.product(*[sorted(marginals[trait]) for trait in traits]):
        expected = n * float(np.prod([marginals[trait][state] for trait, state in zip(traits, combo)]))
        key = combo if len(traits) > 1 else combo[0]
        obs = int(observed.get(key, 0))
        rows.append({
            **{trait: state for trait, state in zip(traits, combo)},
            "observed_species": obs,
            "expected_species_if_independent": round(expected, 3),
            "forbidden": bool(obs == 0 and expected >= min_expected),
            "depleted": bool(obs > 0 and expected > 0 and obs / expected < 0.5),
        })
    return pd.DataFrame(rows).sort_values(["forbidden", "expected_species_if_independent"], ascending=[False, False])


def main() -> None:
    args = parse_args()
    long = pd.read_csv(args.observation_long, dtype=str, keep_default_na=False)
    if missing := LONG_REQUIRED.difference(long.columns):
        raise ValueError(f"Observation-long table missing columns: {sorted(missing)}")
    state_col = f"observation_ai_{args.measurement_mode}_state"
    if state_col not in long.columns:
        raise ValueError(f"Observation-long table missing measurement column: {state_col}")
    ontology_traits = set(pd.read_csv(args.ontology, dtype=str, keep_default_na=False)["trait_id"].map(text))

    present = [t for t in sorted(long["trait_id"].map(text).unique()) if t]
    traits = args.traits if args.traits else present
    traits = [t for t in traits if t in ontology_traits and t in present]
    if len(traits) < 2:
        raise ValueError(f"Need >=2 analysable traits; got {traits}")

    rng = np.random.default_rng(args.seed)
    matrix = build_species_matrix(long, state_col, traits, args.min_obs_per_species_trait)
    # keep only traits still varying after aggregation
    traits = [t for t in traits if matrix[t].nunique(dropna=True) >= 2]
    matrix = matrix[traits].dropna(axis=0, how="any") if traits else matrix
    if len(traits) < 2 or len(matrix) < args.min_species:
        raise ValueError(f"Insufficient complete-case data: {len(matrix)} species x {len(traits)} varying traits "
                         f"(need >= {args.min_species} species and >= 2 traits)")

    observed_mean_v, pair_table = mean_pairwise_association(matrix, traits)
    null_v = permutation_null(matrix, traits, args.n_permutations, rng)
    integration_p = float((np.sum(null_v >= observed_mean_v) + 1) / (args.n_permutations + 1))

    observed_combos = distinct_combinations(matrix, traits)
    null_combos = combination_null(matrix, traits, args.n_permutations, rng)
    constraint_p = float((np.sum(null_combos <= observed_combos) + 1) / (args.n_permutations + 1))

    forbidden = forbidden_combinations(matrix, traits, args.min_expected_forbidden, args.max_forbidden_enumeration)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    matrix.reset_index().to_csv(out / "species_trait_modal_matrix.csv", index=False, encoding="utf-8-sig")
    pair_table.to_csv(out / "trait_pairwise_association_cramers_v.csv", index=False, encoding="utf-8-sig")
    forbidden.to_csv(out / "trait_combination_occupancy.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "measurement_mode": args.measurement_mode,
        "n_species_complete_case": int(len(matrix)),
        "traits": traits,
        "integration": {
            "mean_pairwise_cramers_v": observed_mean_v,
            "null_mean": float(np.mean(null_v)),
            "permutation_p_one_sided_high": integration_p,
            "interpretation": "Lower p = traits more coordinated than independent-shuffle null (phenotypic integration).",
        },
        "constraint": {
            "n_distinct_combinations_observed": observed_combos,
            "null_mean_distinct_combinations": float(np.mean(null_combos)),
            "permutation_p_one_sided_low": constraint_p,
            "n_forbidden_combinations": int(forbidden["forbidden"].sum()) if "forbidden" in forbidden.columns else None,
            "interpretation": "Lower p = fewer realised combinations than independence (constraint / forbidden combinations).",
        },
        "n_permutations": args.n_permutations,
        "semantic_status": "Descriptive species-level trait integration and constraint for automated image traits. "
        "Not evidence of adaptation, correlated evolution, or selection; each result is a hypothesis handed to Chapter 2.",
    }
    (out / "trait_integration_constraint_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

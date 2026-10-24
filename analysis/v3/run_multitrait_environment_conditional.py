#!/usr/bin/env python3
"""Exploratory climate-conditioned capitulum integration on the frozen common cohort.

This is an extension, not a refit of the manuscript's frozen trait-environment
tests.  The original nine CHELSA/process predictors are represented by three
outcome-blind environmental PCs; BIO18 and elevation are NOT in the frozen input.

The response is the same 9 constructs x 36 pairs in the 1,734 observation /
42 taxon complete-18 cohort. RV describes non-negative dependence strength and
must never be interpreted as a signed correlation or a causal mechanism.

Two scales are deliberately separate:
  * within taxon: demean traits and environment within each taxon, then
    taxon-balanced weighted least-squares adjustment;
  * among taxa: regress taxon-median traits on environmental PCs calculated
    from taxon-median environments (one row per taxon).

Matched-n resampling draws one residualized within-taxon observation per taxon,
so each compared RV is estimated from 42 rows. These are post-hoc descriptive
sensitivities, NOT individual-pair significance tests.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v3 import run_biological_axis_reanalysis as axis
from analysis.v3 import run_construct_scale_upgrade as upgrade


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traits", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--reference-pairs", type=Path,
                   default=Path("reproducibility/current_reference/upgrade/complete18_construct_pairwise.csv"))
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--environment-pcs", type=int, default=3)
    p.add_argument("--matched-replicates", type=int, default=500)
    p.add_argument("--seed", type=int, default=20261024)
    return p.parse_args()


def weighted_pcs(x: np.ndarray, weights: np.ndarray, n_components: int) -> tuple[np.ndarray, list[float]]:
    x = np.asarray(x, float)
    w = np.asarray(weights, float)
    if x.ndim != 2 or len(x) != len(w) or not np.isfinite(x).all():
        raise ValueError("environment matrix and row weights must be finite")
    if (w <= 0).any():
        raise ValueError("row weights must be positive")
    w = w / w.sum()
    centered = x - np.average(x, weights=w, axis=0)
    std = np.sqrt(np.sum(w[:, None] * centered ** 2, axis=0))
    if (std <= 1e-12).any():
        raise ValueError("constant predictor in common cohort")
    z = centered / std
    _, singular, vt = np.linalg.svd(z * np.sqrt(w[:, None]), full_matrices=False)
    if n_components < 1 or n_components > np.sum(singular > 1e-10):
        raise ValueError("invalid PC rank")
    pcs = z @ vt[:n_components, :].T
    explained = (singular ** 2) / np.sum(singular ** 2)
    return pcs, explained[:n_components].tolist()


def weighted_residual(y: np.ndarray, pcs: np.ndarray, weights: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    if y.ndim == 1:
        y = y[:, None]
    w = np.asarray(weights, float)
    design = np.column_stack([np.ones(len(pcs)), pcs])
    root = np.sqrt(w / w.sum())[:, None]
    beta = np.linalg.lstsq(design * root, y * root, rcond=None)[0]
    return y - design @ beta


def conditioned_data(frame: pd.DataFrame, env_columns: list[str], n_pcs: int
                    ) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[float]]]:
    fmap = frame.attrs["feature_map"]
    feature_columns = list(dict.fromkeys(col for name in upgrade.CORE for col in fmap[name]))
    counts = frame.groupby("taxon_name").size()
    weights = frame.taxon_name.map(1.0 / counts).to_numpy(float)

    taxon_medians = frame.groupby("taxon_name")[[*feature_columns, *env_columns]].median()
    among_pcs, among_var = weighted_pcs(
        taxon_medians[env_columns].to_numpy(float),
        np.ones(len(taxon_medians)), n_pcs)
    among = taxon_medians[feature_columns].copy()
    for j, col in enumerate(feature_columns):
        among[col] = weighted_residual(
            taxon_medians[[col]].to_numpy(float), among_pcs,
            np.ones(len(taxon_medians)))[:, 0]

    # Demeaning before fitting removes taxon intercepts, not the distinct
    # geographic/photographic biases that require other sensitivity analyses.
    within = frame[["obs_id", "taxon_name", *feature_columns]].copy()
    features_c = (frame[feature_columns] -
                  frame.groupby("taxon_name")[feature_columns].transform("mean"))
    env_c = (frame[env_columns] -
             frame.groupby("taxon_name")[env_columns].transform("mean"))
    within_pcs, within_var = weighted_pcs(env_c.to_numpy(float), weights, n_pcs)
    adjusted = weighted_residual(features_c.to_numpy(float), within_pcs, weights)
    within[feature_columns] = adjusted
    return within, among, {
        "within_taxon_pc_explained": within_var,
        "among_taxon_pc_explained": among_var,
    }


def conditioned_pairs(
    within: pd.DataFrame, among: pd.DataFrame, fmap: dict[str, list[str]]
) -> pd.DataFrame:
    taxa = within.taxon_name.to_numpy(dtype=object)
    n_by_taxon = within.groupby("taxon_name").size()
    weights = within.taxon_name.map(1.0 / n_by_taxon).to_numpy(float)
    rows = []
    for i, left in enumerate(upgrade.CORE):
        for right in upgrade.CORE[i + 1:]:
            w = upgrade.rv(within[fmap[left]].to_numpy(float),
                           within[fmap[right]].to_numpy(float), weights)
            a = upgrade.rv(among[fmap[left]].to_numpy(float),
                           among[fmap[right]].to_numpy(float))
            rows.append({
                "construct_left": left, "construct_right": right,
                "conditioned_within_rv": w, "conditioned_among_rv": a,
                "conditioned_among_minus_within": a - w,
            })
    return pd.DataFrame(rows)


def matched_n(
    within: pd.DataFrame, among: pd.DataFrame, fmap: dict[str, list[str]],
    pair_table: pd.DataFrame, reps: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if reps < 50:
        raise ValueError("at least 50 matched replicates required")
    taxa = sorted(within.taxon_name.unique())
    if taxa != list(among.index):
        # Groupby sorts in the canonical analysis. No implicit reordering.
        raise ValueError("within/among taxon cohorts differ")
    positions = {t: np.flatnonzero(within.taxon_name.to_numpy() == t) for t in taxa}
    rng = np.random.default_rng(seed)
    candidates = [(l, r) for i, l in enumerate(upgrade.CORE)
                  for r in upgrade.CORE[i + 1:]]
    fixed = {(r.construct_left, r.construct_right): float(r.conditioned_among_rv)
             for r in pair_table.itertuples()}
    deltas = {key: [] for key in candidates}
    rows = []
    for rep in range(reps):
        selected = np.array([rng.choice(positions[t]) for t in taxa], dtype=int)
        sample = within.iloc[selected]
        within_rvs = []
        for l, r in candidates:
            val = upgrade.rv(sample[fmap[l]].to_numpy(float),
                             sample[fmap[r]].to_numpy(float))
            within_rvs.append(val)
            deltas[(l, r)].append(fixed[(l, r)] - val)
        among_rvs = np.array([fixed[key] for key in candidates])
        within_rvs = np.array(within_rvs)
        rows.append({
            "replicate": rep,
            "median_matched_within_rv": float(np.median(within_rvs)),
            "median_fixed_among_rv": float(np.median(among_rvs)),
            "among_minus_within_median_rv": float(np.median(among_rvs) -
                                                   np.median(within_rvs)),
            "pairs_stronger_among": int(np.sum(among_rvs > within_rvs)),
        })
    summary = []
    for (l, r), values in deltas.items():
        arr = np.asarray(values, float)
        summary.append({
            "construct_left": l, "construct_right": r,
            "matched_delta_median": float(np.median(arr)),
            "matched_delta_low95": float(np.quantile(arr, 0.025)),
            "matched_delta_high95": float(np.quantile(arr, 0.975)),
            "fraction_delta_positive": float(np.mean(arr > 0)),
        })
    return pd.DataFrame(rows), pd.DataFrame(summary)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    t, e = axis.load_data(args.traits, args.environment)
    common, endpoints = upgrade.common_complete_features(t, 5)
    if (len(common), common.taxon_name.nunique(), len(endpoints)) != (1734, 42, 18):
        raise SystemExit("fail closed: wrong frozen complete-18 common cohort")
    fmap = common.attrs["feature_map"]
    _, _, baseline = upgrade.matrices_from_common(common)
    frozen = pd.read_csv(args.reference_pairs)
    keys = ["construct_left", "construct_right"]
    merged = baseline.merge(frozen, on=keys, validate="one_to_one",
                            suffixes=("_rebuilt", "_frozen"))
    if len(merged) != 36:
        raise SystemExit("frozen pair table does not cover 36 pairs")
    for column in ("within_taxon_rv", "among_taxon_rv"):
        if not np.allclose(merged[f"{column}_rebuilt"], merged[f"{column}_frozen"],
                           atol=1e-9, rtol=1e-9):
            raise SystemExit(f"frozen {column} mismatch")

    env_columns = axis.PREDICTORS.copy()
    cols = ["obs_id", "taxon_name", *env_columns]
    e = e[cols].copy()
    if e[["obs_id", "taxon_name"]].duplicated().any():
        raise SystemExit("duplicate environment observation/taxon pairs")
    frame = common.merge(e, on=["obs_id", "taxon_name"], how="left",
                         validate="one_to_one")
    frame.attrs["feature_map"] = fmap
    if len(frame) != len(common) or frame[env_columns].isna().any().any():
        raise SystemExit("missing environmental values in common cohort")
    within, among, pc_info = conditioned_data(frame, env_columns, args.environment_pcs)
    conditional = conditioned_pairs(within, among, fmap)
    matched, edge_intervals = matched_n(
        within, among, fmap, conditional, args.matched_replicates, args.seed)
    combined = baseline.merge(conditional, on=keys, validate="one_to_one")
    combined = combined.merge(edge_intervals, on=keys, validate="one_to_one")

    combined.to_csv(args.out_dir / "multitrait_pairwise_conditional.csv", index=False)
    matched.to_csv(args.out_dir / "multitrait_matched_n_replicates.csv", index=False)
    delta = matched.among_minus_within_median_rv
    report = {
        "analysis_id": "ch1_v3_multitrait_existing_climate_pc_sensitivity_20261024",
        "status": "new_posthoc_exploratory_sensitivity_not_frozen_main_claim",
        "cohort": {"observations": len(common), "taxa": 42, "constructs": 9,
                   "pairs": 36, "endpoints": 18},
        "predictors_adjusted": env_columns,
        "not_available_in_frozen_environment": ["elevation", "BIO18"],
        "adjustment": {
            "method": f"{args.environment_pcs} outcome-blind weighted environmental PCs fitted separately at within and among taxon scales",
            "environment_pcs": args.environment_pcs, **pc_info,
            "geography_phylogeny_measurement_error_adjusted": False,
        },
        "frozen_raw": {
            "median_within_rv": float(baseline.within_taxon_rv.median()),
            "median_among_rv": float(baseline.among_taxon_rv.median()),
            "pairs_stronger_among": int((baseline.among_taxon_rv >
                                         baseline.within_taxon_rv).sum()),
        },
        "climate_conditioned_raw": {
            "median_within_rv": float(conditional.conditioned_within_rv.median()),
            "median_among_rv": float(conditional.conditioned_among_rv.median()),
            "pairs_stronger_among": int((conditional.conditioned_among_rv >
                                         conditional.conditioned_within_rv).sum()),
        },
        "climate_conditioned_equal_n": {
            "replicates": args.matched_replicates,
            "median_among_minus_within": float(delta.median()),
            "difference_low95": float(delta.quantile(0.025)),
            "difference_high95": float(delta.quantile(0.975)),
            "fraction_positive": float((delta > 0).mean()),
            "median_pairs_stronger_among": float(matched.pairs_stronger_among.median()),
        },
        "interpretation_boundary": (
            "Post-hoc, descriptive residual RV sensitivity to three orthogonal "
            "environmental axes. Neither causality, spatial/phylogenetic "
            "independence nor pair-specific significance is established. "
            "RV is unsigned. Raw RV depends on sample size and dimension; "
            "interpret alongside matched-n sensitivity and the existing "
            "permutation-null/measurement audits. No BIO18 or elevation results."
        ),
    }
    (args.out_dir / "multitrait_conditional_report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

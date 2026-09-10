#!/usr/bin/env python3
"""Test whether biological construct integration is conserved across scales.

This analysis is a direct extension of the frozen v2 complete-18 whole-capitulum
synthesis.  The same 18 non-surface endpoints are represented as nine biological
constructs defined without reference to environmental outcomes.  For every pair
of constructs we use the exact same merged observation/taxon cohort at both
scales, requiring at least five paired observations per taxon and at least 20
taxa.  Among-taxon integration is computed from taxon medians.  Within-taxon
integration is computed after taxon centring with equal total weight per taxon.

Association strength between two potentially multivariate constructs is the RV
coefficient.  This avoids forcing hue, involucre form or projection pattern onto
arbitrary single PCs.  Matrix concordance is the Spearman correlation between the
upper triangles of the within- and among-taxon RV matrices.  A QAP-style label
permutation tests whether the named construct-to-construct organization is more
aligned across scales than expected under construct relabelling.

The environmental-signature comparison is descriptive only.  It compares the
relative effect-magnitude profiles across the same nine environmental gradients
for each construct; correlated predictors make predictor-label permutation an
inappropriate confirmatory null.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analysis.v3 import run_biological_axis_reanalysis as axis

CORE_CONSTRUCTS = [
    "presentation_angle",
    "floral_lightness",
    "floral_chroma",
    "floral_hue",
    "head_elongation",
    "head_compactness",
    "involucre_form",
    "projection_prominence",
    "projection_pattern",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--traits", required=True, type=Path)
    parser.add_argument("--environment", required=True, type=Path)
    parser.add_argument("--axis-among", required=True, type=Path)
    parser.add_argument("--axis-within", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--minimum-paired-observations-per-taxon", type=int, default=5)
    parser.add_argument("--minimum-taxa", type=int, default=20)
    parser.add_argument("--qap-permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260910)
    return parser.parse_args()


def weighted_covariance(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    w = w / float(w.sum())
    mean = np.sum(x * w[:, None], axis=0)
    centred = x - mean
    return (centred * w[:, None]).T @ centred


def rv_coefficient(x: np.ndarray, y: np.ndarray, weights: np.ndarray | None = None) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if weights is None:
        weights = np.ones(len(x), dtype=float)
    weights = np.asarray(weights, dtype=float)
    w = weights / float(weights.sum())
    mx = np.sum(x * w[:, None], axis=0)
    my = np.sum(y * w[:, None], axis=0)
    xc = x - mx
    yc = y - my
    sxx = weighted_covariance(x, weights)
    syy = weighted_covariance(y, weights)
    sxy = (xc * w[:, None]).T @ yc
    numerator = float(np.trace(sxy @ sxy.T))
    denominator = float(math.sqrt(np.trace(sxx @ sxx) * np.trace(syy @ syy)))
    return numerator / denominator if denominator > 1e-15 else float("nan")


def construct_observation_features(
    traits: pd.DataFrame,
    construct_id: str,
) -> pd.DataFrame:
    definition = axis.CONSTRUCTS[construct_id]
    members = definition["members"]
    part = traits[traits.endpoint_id.isin(members)][
        ["obs_id", "taxon_name", "endpoint_id", "value"]
    ]
    wide = part.pivot(
        index=["obs_id", "taxon_name"], columns="endpoint_id", values="value"
    ).reset_index()
    wide.columns.name = None
    wide = wide.dropna(subset=members)
    out = wide[["obs_id", "taxon_name"]].copy()
    if definition["kind"] == "scalar":
        means, sds = axis.scaling_from_taxon_medians(traits, members, 5)
        signs = np.asarray(definition["signs"], dtype=float)
        score = ((wide[members] - means) / sds).to_numpy(float) @ signs / len(members)
        out[f"{construct_id}__0"] = score
    else:
        for index, member in enumerate(members):
            out[f"{construct_id}__{index}"] = wide[member].to_numpy(float)
    return out


def pairwise_integration(
    observations: dict[str, pd.DataFrame],
    minimum_per_taxon: int,
    minimum_taxa: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for left_index, left in enumerate(CORE_CONSTRUCTS):
        for right in CORE_CONSTRUCTS[left_index + 1 :]:
            merged = observations[left].merge(
                observations[right], on=["obs_id", "taxon_name"], how="inner"
            )
            left_columns = [c for c in merged.columns if c.startswith(left + "__")]
            right_columns = [c for c in merged.columns if c.startswith(right + "__")]
            counts = merged.groupby("taxon_name").size()
            taxa = counts[counts >= minimum_per_taxon].index
            merged = merged[merged.taxon_name.isin(taxa)].copy()
            n_taxa = int(merged.taxon_name.nunique())
            n_observations = int(len(merged))
            base = {
                "construct_left": left,
                "construct_right": right,
                "n_taxa": n_taxa,
                "n_paired_observations": n_observations,
                "minimum_paired_observations_per_taxon": minimum_per_taxon,
            }
            if n_taxa < minimum_taxa or n_observations < 100:
                rows.append({**base, "status": "insufficient_support"})
                continue

            # Among scale: taxon medians on the exact same merged observation scope.
            medians = merged.groupby("taxon_name")[left_columns + right_columns].median()
            among_rv = rv_coefficient(
                medians[left_columns].to_numpy(float),
                medians[right_columns].to_numpy(float),
            )

            # Within scale: remove taxon means and give every retained taxon equal total weight.
            left_values = merged[left_columns]
            right_values = merged[right_columns]
            left_centered = left_values - left_values.groupby(merged.taxon_name).transform("mean")
            right_centered = right_values - right_values.groupby(merged.taxon_name).transform("mean")
            retained_counts = merged.groupby("taxon_name").size()
            weights = merged.taxon_name.map(1.0 / retained_counts).to_numpy(float)
            within_rv = rv_coefficient(
                left_centered.to_numpy(float),
                right_centered.to_numpy(float),
                weights,
            )
            rows.append(
                {
                    **base,
                    "status": "ok",
                    "within_taxon_rv": within_rv,
                    "among_taxon_rv": among_rv,
                    "delta_among_minus_within": among_rv - within_rv,
                }
            )
    return pd.DataFrame(rows)


def matrix_from_pairs(pairs: pd.DataFrame, value: str) -> pd.DataFrame:
    matrix = pd.DataFrame(
        np.eye(len(CORE_CONSTRUCTS)),
        index=CORE_CONSTRUCTS,
        columns=CORE_CONSTRUCTS,
        dtype=float,
    )
    for _, row in pairs[pairs.status.eq("ok")].iterrows():
        left = row.construct_left
        right = row.construct_right
        matrix.loc[left, right] = matrix.loc[right, left] = float(row[value])
    return matrix


def qap_matrix_alignment(
    within: pd.DataFrame,
    among: pd.DataFrame,
    permutations: int,
    seed: int,
) -> dict[str, float | int]:
    if within.isna().any().any() or among.isna().any().any():
        raise ValueError("QAP requires a complete construct matrix")
    upper = np.triu_indices(len(CORE_CONSTRUCTS), 1)
    w = within.to_numpy(float)
    a = among.to_numpy(float)
    observed = float(spearmanr(w[upper], a[upper]).statistic)
    rng = np.random.default_rng(seed)
    exceed = 0
    simulated: list[float] = []
    for _ in range(permutations):
        order = rng.permutation(len(CORE_CONSTRUCTS))
        permuted = a[np.ix_(order, order)]
        statistic = float(spearmanr(w[upper], permuted[upper]).statistic)
        simulated.append(statistic)
        exceed += int(statistic >= observed - 1e-15)
    return {
        "spearman_upper_triangle": observed,
        "qap_p_value_one_sided": float((exceed + 1) / (permutations + 1)),
        "qap_permutations": permutations,
        "qap_null_mean": float(np.mean(simulated)),
        "qap_null_low_95": float(np.quantile(simulated, 0.025)),
        "qap_null_high_95": float(np.quantile(simulated, 0.975)),
    }


def construct_divergence(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for construct in CORE_CONSTRUCTS:
        part = pairs[
            pairs.status.eq("ok")
            & (
                pairs.construct_left.eq(construct)
                | pairs.construct_right.eq(construct)
            )
        ]
        rows.append(
            {
                "construct_id": construct,
                "n_relations": int(len(part)),
                "mean_within_taxon_rv": float(part.within_taxon_rv.mean()),
                "mean_among_taxon_rv": float(part.among_taxon_rv.mean()),
                "mean_delta_among_minus_within": float(part.delta_among_minus_within.mean()),
                "mean_absolute_scale_difference": float(part.delta_among_minus_within.abs().mean()),
            }
        )
    return pd.DataFrame(rows)


def environmental_signature_alignment(
    among_path: Path,
    within_path: Path,
) -> tuple[pd.DataFrame, dict[str, float | int]]:
    among = pd.read_csv(among_path)
    within = pd.read_csv(within_path)
    rows = []
    among_matrix = []
    within_matrix = []
    for construct in CORE_CONSTRUCTS:
        left = among[among.construct_id.eq(construct)][
            ["predictor", "effect_magnitude"]
        ].rename(columns={"effect_magnitude": "among_effect_magnitude"})
        right = within[within.construct_id.eq(construct)][
            ["predictor", "effect_magnitude"]
        ].rename(columns={"effect_magnitude": "within_effect_magnitude"})
        merged = left.merge(right, on="predictor", validate="one_to_one").sort_values("predictor")
        av = merged.among_effect_magnitude.to_numpy(float)
        wv = merged.within_effect_magnitude.to_numpy(float)
        among_matrix.append(av / np.linalg.norm(av))
        within_matrix.append(wv / np.linalg.norm(wv))
        rows.append(
            {
                "construct_id": construct,
                "n_predictors": int(len(merged)),
                "effect_magnitude_signature_spearman": float(spearmanr(av, wv).statistic),
                "effect_magnitude_signature_cosine": float(
                    np.dot(av, wv) / (np.linalg.norm(av) * np.linalg.norm(wv))
                ),
                "strongest_among_predictor": str(merged.iloc[int(np.argmax(av))].predictor),
                "strongest_within_predictor": str(merged.iloc[int(np.argmax(wv))].predictor),
                "same_strongest_predictor": bool(np.argmax(av) == np.argmax(wv)),
            }
        )
    a = np.asarray(among_matrix, dtype=float)
    w = np.asarray(within_matrix, dtype=float)
    summary = {
        "constructs_compared": len(CORE_CONSTRUCTS),
        "same_strongest_predictor_constructs": int(sum(row["same_strongest_predictor"] for row in rows)),
        "row_normalized_signature_spearman_flattened": float(spearmanr(a.ravel(), w.ravel()).statistic),
        "row_normalized_signature_cosine_flattened": float(
            np.dot(a.ravel(), w.ravel()) / (np.linalg.norm(a) * np.linalg.norm(w))
        ),
        "interpretation": "descriptive_only_correlated_environment_predictors_no_predictor_label_null",
    }
    return pd.DataFrame(rows), summary


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    _, traits, _, _ = axis.load_data(args.traits, args.environment)
    observations = {
        construct: construct_observation_features(traits, construct)
        for construct in CORE_CONSTRUCTS
    }
    pairs = pairwise_integration(
        observations,
        args.minimum_paired_observations_per_taxon,
        args.minimum_taxa,
    )
    expected_pairs = len(CORE_CONSTRUCTS) * (len(CORE_CONSTRUCTS) - 1) // 2
    if int(pairs.status.eq("ok").sum()) != expected_pairs:
        raise SystemExit(
            f"complete-18 construct matrix incomplete: {int(pairs.status.eq('ok').sum())}/{expected_pairs} pairs"
        )
    within_matrix = matrix_from_pairs(pairs, "within_taxon_rv")
    among_matrix = matrix_from_pairs(pairs, "among_taxon_rv")
    qap = qap_matrix_alignment(
        within_matrix,
        among_matrix,
        args.qap_permutations,
        args.seed,
    )
    divergence = construct_divergence(pairs)
    signatures, signature_summary = environmental_signature_alignment(
        args.axis_among,
        args.axis_within,
    )

    pairs.to_csv(args.out_dir / "construct_pairwise_integration.csv", index=False)
    within_matrix.to_csv(args.out_dir / "construct_integration_within_matrix.csv")
    among_matrix.to_csv(args.out_dir / "construct_integration_among_matrix.csv")
    divergence.to_csv(args.out_dir / "construct_scale_divergence.csv", index=False)
    signatures.to_csv(args.out_dir / "construct_environment_signature_alignment.csv", index=False)

    top_among = pairs.sort_values("delta_among_minus_within", ascending=False).head(5)
    top_within = pairs.sort_values("delta_among_minus_within", ascending=True).head(5)
    report = {
        "analysis_id": "ch1_v3_construct_scale_integration_20260910",
        "source_endpoint_scope": "frozen_v2_complete18_non_surface_endpoints",
        "constructs": CORE_CONSTRUCTS,
        "n_constructs": len(CORE_CONSTRUCTS),
        "n_pairwise_relations": expected_pairs,
        "pairwise_cohort_rule": "same paired observations and taxa at both scales; >=5 paired observations per taxon; >=20 taxa",
        "within_rule": "taxon-centred construct components with equal total weight per taxon",
        "among_rule": "taxon medians on the exact same pairwise observation scope",
        "integration_metric": "RV_coefficient",
        "matrix_alignment": qap,
        "median_within_taxon_rv": float(pairs.within_taxon_rv.median()),
        "median_among_taxon_rv": float(pairs.among_taxon_rv.median()),
        "n_relations_stronger_among_than_within": int((pairs.among_taxon_rv > pairs.within_taxon_rv).sum()),
        "top_among_strengthening_relations": top_among[
            ["construct_left", "construct_right", "within_taxon_rv", "among_taxon_rv", "delta_among_minus_within"]
        ].to_dict("records"),
        "top_within_strengthening_relations": top_within[
            ["construct_left", "construct_right", "within_taxon_rv", "among_taxon_rv", "delta_among_minus_within"]
        ].to_dict("records"),
        "environment_signature_alignment": signature_summary,
        "claim_boundary": (
            "construct-level scale-dependent integration synthesis; not functional/genetic modularity, "
            "not plasticity, and not an independent causal environmental analysis"
        ),
    }
    (args.out_dir / "construct_scale_integration_report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

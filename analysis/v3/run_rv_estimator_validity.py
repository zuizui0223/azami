#!/usr/bin/env python3
"""Post-hoc estimator-validity checks for the Chapter 1 RV scale contrast.

The current manuscript compares RV integration among 42 taxon medians with RV
integration within 1,734 observations from the same 42 taxa. Classical RV has a
positive finite-sample baseline that depends on sample size and block dimension.
This script does not redefine the frozen raw-RV result. It adds three bounded
sensitivity diagnostics:

1. equal-n within resampling: one centred observation per taxon (n=42) is drawn
   repeatedly and compared with the fixed among-taxon median matrix;
2. permutation-null centring: chance RV is estimated separately within and
   among taxa for every construct pair and subtracted from observed RV;
3. joint-coordinate standardization: construct coordinates are scaled by their
   across-taxon-median SD before RV is recomputed.

These are estimator diagnostics, not new ecological discovery tests.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v3 import run_biological_axis_reanalysis as axis
from analysis.v3 import run_construct_scale_upgrade as upgrade

CORE = upgrade.CORE


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--traits", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--minimum-complete-observations-per-taxon", type=int, default=5)
    p.add_argument("--equal-n-replicates", type=int, default=1000)
    p.add_argument("--null-permutations", type=int, default=499)
    p.add_argument("--qap-permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260915)
    return p.parse_args()


def stable_seed(seed: int, *parts: object) -> int:
    digest = hashlib.sha256("|".join([str(seed), *map(str, parts)]).encode()).digest()
    return int.from_bytes(digest[:8], "little")


def feature_map(frame: pd.DataFrame) -> dict[str, list[str]]:
    fmap = frame.attrs.get("feature_map")
    if not isinstance(fmap, dict):
        raise ValueError("missing feature_map on common-cohort frame")
    return fmap


def centred_common(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    fmap = feature_map(frame)
    all_features = [c for construct in CORE for c in fmap[construct]]
    for col in all_features:
        out[col] = frame[col] - frame[col].groupby(frame.taxon_name).transform("mean")
    out.attrs["feature_map"] = fmap
    return out


def fixed_among_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    _, _, pairs = upgrade.matrices_from_common(frame)
    return pairs.copy()


def equal_n_within_resampling(
    frame: pd.DataFrame,
    among_pairs: pd.DataFrame,
    replicates: int,
    seed: int,
) -> pd.DataFrame:
    centred = centred_common(frame)
    fmap = feature_map(frame)
    taxa = np.array(sorted(frame.taxon_name.unique()), dtype=object)
    group_indices = {
        taxon: np.flatnonzero(frame.taxon_name.to_numpy(dtype=object) == taxon)
        for taxon in taxa
    }
    among_lookup = {
        (row.construct_left, row.construct_right): float(row.among_taxon_rv)
        for row in among_pairs.itertuples()
    }
    fixed_among_median = float(among_pairs.among_taxon_rv.median())
    rows: list[dict[str, object]] = []
    rng = np.random.default_rng(seed)

    for rep in range(replicates):
        selected = [int(rng.choice(group_indices[taxon])) for taxon in taxa]
        sample = centred.iloc[selected]
        within_values: list[float] = []
        stronger_among = 0
        for i, left in enumerate(CORE):
            lx = fmap[left]
            for right in CORE[i + 1 :]:
                rx = fmap[right]
                within_rv = upgrade.rv(
                    sample[lx].to_numpy(float),
                    sample[rx].to_numpy(float),
                )
                within_values.append(within_rv)
                stronger_among += int(among_lookup[(left, right)] > within_rv)
        within_median = float(np.median(within_values))
        rows.append({
            "replicate": rep,
            "n_rows_within": len(taxa),
            "n_taxa": len(taxa),
            "matched_within_median_rv": within_median,
            "fixed_among_median_rv": fixed_among_median,
            "among_minus_within_median_rv": fixed_among_median - within_median,
            "relations_stronger_among": stronger_among,
        })
    return pd.DataFrame(rows)


def permutation_null_centered_pairs(
    frame: pd.DataFrame,
    observed_pairs: pd.DataFrame,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    fmap = feature_map(frame)
    centred = centred_common(frame)
    taxa = np.array(sorted(frame.taxon_name.unique()), dtype=object)
    medians = frame.groupby("taxon_name")[[c for v in fmap.values() for c in v]].median().reindex(taxa)
    counts = frame.groupby("taxon_name").size()
    weights = frame.taxon_name.map(1.0 / counts).to_numpy(float)
    group_indices = [
        np.flatnonzero(frame.taxon_name.to_numpy(dtype=object) == taxon)
        for taxon in taxa
    ]
    observed_lookup = {
        (row.construct_left, row.construct_right): row
        for row in observed_pairs.itertuples()
    }
    rows: list[dict[str, object]] = []

    for pair_index, left in enumerate(CORE):
        for right in CORE[pair_index + 1 :]:
            lx = fmap[left]
            rx = fmap[right]
            observed = observed_lookup[(left, right)]
            left_within = centred[lx].to_numpy(float)
            right_within = centred[rx].to_numpy(float)
            left_among = medians[lx].to_numpy(float)
            right_among = medians[rx].to_numpy(float)

            rng = np.random.default_rng(stable_seed(seed, left, right))
            within_null = np.empty(permutations, dtype=float)
            among_null = np.empty(permutations, dtype=float)
            base_index = np.arange(len(frame))
            for b in range(permutations):
                perm_index = base_index.copy()
                for idx in group_indices:
                    perm_index[idx] = rng.permutation(idx)
                within_null[b] = upgrade.rv(
                    left_within,
                    right_within[perm_index],
                    weights,
                )
                among_null[b] = upgrade.rv(
                    left_among,
                    right_among[rng.permutation(len(taxa))],
                )

            within_null_mean = float(np.mean(within_null))
            among_null_mean = float(np.mean(among_null))
            within_corrected = float(observed.within_taxon_rv) - within_null_mean
            among_corrected = float(observed.among_taxon_rv) - among_null_mean
            rows.append({
                "construct_left": left,
                "construct_right": right,
                "observed_within_rv": float(observed.within_taxon_rv),
                "observed_among_rv": float(observed.among_taxon_rv),
                "within_null_mean": within_null_mean,
                "among_null_mean": among_null_mean,
                "within_null_low95": float(np.quantile(within_null, 0.025)),
                "within_null_high95": float(np.quantile(within_null, 0.975)),
                "among_null_low95": float(np.quantile(among_null, 0.025)),
                "among_null_high95": float(np.quantile(among_null, 0.975)),
                "null_centered_within_rv": within_corrected,
                "null_centered_among_rv": among_corrected,
                "null_centered_delta_among_minus_within": among_corrected - within_corrected,
            })
    return pd.DataFrame(rows)


def standardize_coordinates_by_taxon_medians(frame: pd.DataFrame) -> pd.DataFrame:
    fmap = feature_map(frame)
    out = frame.copy()
    for construct in CORE:
        for col in fmap[construct]:
            taxon_medians = frame.groupby("taxon_name")[col].median()
            mean = float(taxon_medians.mean())
            sd = float(taxon_medians.std(ddof=0))
            if not np.isfinite(sd) or sd <= 0:
                raise ValueError(f"non-finite or zero taxon-median SD for {col}")
            out[col] = (frame[col] - mean) / sd
    out.attrs["feature_map"] = fmap
    return out


def summarize(
    raw_pairs: pd.DataFrame,
    matched: pd.DataFrame,
    null_pairs: pd.DataFrame,
    standardized_pairs: pd.DataFrame,
    standardized_within: pd.DataFrame,
    standardized_among: pd.DataFrame,
    standardized_qap: dict[str, float | int],
    standardized_module_within: dict[str, float | int],
    standardized_module_among: dict[str, float | int],
    args: argparse.Namespace,
) -> dict[str, object]:
    matched_delta = matched.among_minus_within_median_rv.to_numpy(float)
    matched_counts = matched.relations_stronger_among.to_numpy(int)
    equal_n_gate = bool(np.quantile(matched_delta, 0.025) > 0 and np.mean(matched_delta > 0) >= 0.95)
    pair_count_gate = bool(np.mean(matched_counts > 18) >= 0.95)

    null_within_median = float(null_pairs.null_centered_within_rv.median())
    null_among_median = float(null_pairs.null_centered_among_rv.median())
    null_count = int((null_pairs.null_centered_among_rv > null_pairs.null_centered_within_rv).sum())

    std_within_median = float(standardized_pairs.within_taxon_rv.median())
    std_among_median = float(standardized_pairs.among_taxon_rv.median())
    std_count = int((standardized_pairs.among_taxon_rv > standardized_pairs.within_taxon_rv).sum())

    return {
        "analysis_id": "ch1_v3_rv_estimator_validity_20260915",
        "status": "posthoc_estimator_validity_sensitivity",
        "prospective_status": "not prospective; motivated after the raw-RV result was known",
        "common_cohort": {
            "observations": 1734,
            "taxa": 42,
            "constructs": 9,
            "relations": 36,
        },
        "frozen_raw_result": {
            "median_within_rv": float(raw_pairs.within_taxon_rv.median()),
            "median_among_rv": float(raw_pairs.among_taxon_rv.median()),
            "relations_stronger_among": int((raw_pairs.among_taxon_rv > raw_pairs.within_taxon_rv).sum()),
        },
        "equal_n_within_resampling": {
            "replicates": int(args.equal_n_replicates),
            "within_rows_per_replicate": 42,
            "fixed_among_rows": 42,
            "matched_within_median_rv_median": float(np.median(matched.matched_within_median_rv)),
            "matched_within_median_rv_low95": float(np.quantile(matched.matched_within_median_rv, 0.025)),
            "matched_within_median_rv_high95": float(np.quantile(matched.matched_within_median_rv, 0.975)),
            "among_minus_within_median_rv_median": float(np.median(matched_delta)),
            "among_minus_within_median_rv_low95": float(np.quantile(matched_delta, 0.025)),
            "among_minus_within_median_rv_high95": float(np.quantile(matched_delta, 0.975)),
            "probability_positive_median_difference": float(np.mean(matched_delta > 0)),
            "relations_stronger_among_median": float(np.median(matched_counts)),
            "relations_stronger_among_low95": float(np.quantile(matched_counts, 0.025)),
            "relations_stronger_among_high95": float(np.quantile(matched_counts, 0.975)),
            "probability_majority_relations_stronger_among": float(np.mean(matched_counts > 18)),
            "overall_strength_gate_pass": equal_n_gate,
            "pair_count_majority_gate_pass": pair_count_gate,
        },
        "permutation_null_centering": {
            "permutations_per_pair_per_scale": int(args.null_permutations),
            "median_null_centered_within_rv": null_within_median,
            "median_null_centered_among_rv": null_among_median,
            "median_null_centered_among_minus_within": null_among_median - null_within_median,
            "relations_stronger_among_after_null_centering": null_count,
            "relations_total": 36,
            "role": "diagnostic only; no new significance test",
        },
        "joint_coordinate_standardization": {
            "median_within_rv": std_within_median,
            "median_among_rv": std_among_median,
            "relations_stronger_among": std_count,
            "relations_total": 36,
            "matrix_alignment_rho": float(upgrade.matrix_alignment(standardized_within, standardized_among)),
            "qap_p_one_sided": float(standardized_qap["qap_p_one_sided"]),
            "module_within_p_one_sided": float(standardized_module_within["permutation_p_one_sided"]),
            "module_among_p_one_sided": float(standardized_module_among["permutation_p_one_sided"]),
        },
        "manuscript_decision": {
            "retain_stronger_overall_among_language": equal_n_gate,
            "retain_raw_33_of_36_as_bias_robust_claim": False,
            "recommended_boundary": (
                "If the equal-n gate passes, retain the qualitative statement that integration is stronger "
                "overall among taxa, but label raw 33/36 and raw median-RV differences as descriptive values "
                "from the frozen estimator. Report the equal-n sensitivity in Supporting Information."
            ),
        },
        "claim_boundary": (
            "Estimator robustness supports only a descriptive scale contrast in visible image phenotypes. "
            "It does not establish evolutionary increase in integration, genetic/developmental modularity, "
            "or a causal mechanism."
        ),
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    traits, _ = axis.load_data(args.traits, args.environment)
    common, endpoints = upgrade.common_complete_features(
        traits, args.minimum_complete_observations_per_taxon
    )
    if len(common) != 1734 or common.taxon_name.nunique() != 42 or len(endpoints) != 18:
        raise SystemExit(
            f"unexpected common cohort: {len(common)} observations / "
            f"{common.taxon_name.nunique()} taxa / {len(endpoints)} endpoints"
        )

    raw_within, raw_among, raw_pairs = upgrade.matrices_from_common(common)
    matched = equal_n_within_resampling(
        common, raw_pairs, args.equal_n_replicates, stable_seed(args.seed, "equal_n")
    )
    null_pairs = permutation_null_centered_pairs(
        common, raw_pairs, args.null_permutations, stable_seed(args.seed, "null")
    )

    standardized = standardize_coordinates_by_taxon_medians(common)
    std_within, std_among, std_pairs = upgrade.matrices_from_common(standardized)
    std_qap = upgrade.qap_alignment(
        std_within,
        std_among,
        args.qap_permutations,
        np.random.default_rng(stable_seed(args.seed, "std_qap")),
    )
    std_module_within = upgrade.module_cohesion(
        std_within,
        args.qap_permutations,
        np.random.default_rng(stable_seed(args.seed, "std_module_within")),
    )
    std_module_among = upgrade.module_cohesion(
        std_among,
        args.qap_permutations,
        np.random.default_rng(stable_seed(args.seed, "std_module_among")),
    )

    raw_pairs.to_csv(args.out_dir / "raw_common_cohort_pairs.csv", index=False)
    matched.to_csv(args.out_dir / "equal_n_within_resampling.csv", index=False)
    null_pairs.to_csv(args.out_dir / "permutation_null_centered_pairs.csv", index=False)
    std_pairs.to_csv(args.out_dir / "joint_standardized_pairs.csv", index=False)
    std_within.to_csv(args.out_dir / "joint_standardized_within_matrix.csv")
    std_among.to_csv(args.out_dir / "joint_standardized_among_matrix.csv")

    report = summarize(
        raw_pairs,
        matched,
        null_pairs,
        std_pairs,
        std_within,
        std_among,
        std_qap,
        std_module_within,
        std_module_among,
        args,
    )
    (args.out_dir / "rv_estimator_validity_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

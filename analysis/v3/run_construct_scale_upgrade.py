#!/usr/bin/env python3
"""Strengthen the whole-capitulum synthesis without redefining frozen v2 claims.

This analysis retains only the synthesis layers used by the current manuscript:

1. a complete-18 common-cohort replication in which all nine biological
   constructs, all 36 construct pairs, and both biological scales use the exact
   same observations and taxa;
2. taxon-bootstrap uncertainty for the within-vs-among matrix alignment;
3. module-cohesion tests asking whether biologically related constructs are more
   internally integrated than unrelated constructs at each scale;
4. six predeclared v2 environmental-block signatures.

The frozen v2 endpoint atlas, multiplicity family, and two headline ecological
candidates are not modified by this script.
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

CORE = [
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

MODULES = {
    "presentation_angle": "presentation",
    "floral_lightness": "colour",
    "floral_chroma": "colour",
    "floral_hue": "colour",
    "head_elongation": "head_form",
    "head_compactness": "head_form",
    "involucre_form": "involucre_armature",
    "projection_prominence": "involucre_armature",
    "projection_pattern": "involucre_armature",
}

ENV_BLOCKS = {
    "thermal": ["chelsa_bio01", "chelsa_bio04"],
    "hydric": ["chelsa_bio12", "chelsa_bio15"],
    "radiative_atmospheric": ["chelsa_rsds_mean", "chelsa_vpd_mean"],
    "mechanical": ["chelsa_sfcwind_mean"],
    "growing_season_water": ["chelsa_gsp"],
    "resource_productivity": ["chelsa_npp"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--traits", type=Path, required=True)
    p.add_argument("--environment", type=Path, required=True)
    p.add_argument("--axis-among", type=Path, required=True)
    p.add_argument("--axis-within", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--minimum-complete-observations-per-taxon", type=int, default=5)
    p.add_argument("--bootstrap-replicates", type=int, default=1000)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260910)
    return p.parse_args()


def weighted_covariance(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    m = np.sum(x * w[:, None], axis=0)
    c = x - m
    return (c * w[:, None]).T @ c


def rv(x: np.ndarray, y: np.ndarray, weights: np.ndarray | None = None) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if weights is None:
        weights = np.ones(len(x), dtype=float)
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    xc = x - np.sum(x * w[:, None], axis=0)
    yc = y - np.sum(y * w[:, None], axis=0)
    sxx = weighted_covariance(x, w)
    syy = weighted_covariance(y, w)
    sxy = (xc * w[:, None]).T @ yc
    num = float(np.trace(sxy @ sxy.T))
    den = float(math.sqrt(np.trace(sxx @ sxx) * np.trace(syy @ syy)))
    return num / den if den > 1e-15 else float("nan")


def common_complete_features(t: pd.DataFrame, min_per_taxon: int) -> tuple[pd.DataFrame, list[str]]:
    endpoints = sorted({m for c in CORE for m in axis.CONSTRUCTS[c]["members"]})
    if len(endpoints) != 18:
        raise ValueError(f"expected 18 complete-synthesis endpoints, got {len(endpoints)}")
    part = t[t.endpoint_id.isin(endpoints)][["obs_id", "taxon_name", "endpoint_id", "value"]]
    wide = part.pivot(index=["obs_id", "taxon_name"], columns="endpoint_id", values="value").reset_index()
    wide.columns.name = None
    wide = wide.dropna(subset=endpoints).copy()
    counts = wide.groupby("taxon_name").size()
    keep = counts[counts >= min_per_taxon].index
    wide = wide[wide.taxon_name.isin(keep)].copy()

    out = wide[["obs_id", "taxon_name"]].copy()
    feature_map: dict[str, list[str]] = {}
    for construct in CORE:
        d = axis.CONSTRUCTS[construct]
        members = d["members"]
        cols: list[str] = []
        if d["kind"] == "scalar":
            means, sds = axis.scaling_from_taxon_medians(t, members, 5)
            signs = np.asarray(d["signs"], dtype=float)
            score = ((wide[members] - means) / sds).to_numpy(float) @ signs / len(members)
            col = f"{construct}__0"
            out[col] = score
            cols.append(col)
        else:
            for j, member in enumerate(members):
                col = f"{construct}__{j}"
                out[col] = wide[member].to_numpy(float)
                cols.append(col)
        feature_map[construct] = cols
    out.attrs["feature_map"] = feature_map
    return out, endpoints


def matrices_from_common(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fmap: dict[str, list[str]] = frame.attrs["feature_map"]
    within = pd.DataFrame(np.eye(len(CORE)), index=CORE, columns=CORE, dtype=float)
    among = within.copy()
    rows = []
    counts = frame.groupby("taxon_name").size()
    weights = frame.taxon_name.map(1.0 / counts).to_numpy(float)

    for i, left in enumerate(CORE):
        lx = fmap[left]
        for right in CORE[i + 1 :]:
            rx = fmap[right]
            med = frame.groupby("taxon_name")[lx + rx].median()
            a = rv(med[lx].to_numpy(float), med[rx].to_numpy(float))

            lraw = frame[lx]
            rraw = frame[rx]
            lc = lraw - lraw.groupby(frame.taxon_name).transform("mean")
            rc = rraw - rraw.groupby(frame.taxon_name).transform("mean")
            w = rv(lc.to_numpy(float), rc.to_numpy(float), weights)

            among.loc[left, right] = among.loc[right, left] = a
            within.loc[left, right] = within.loc[right, left] = w
            rows.append({
                "construct_left": left,
                "construct_right": right,
                "within_taxon_rv": w,
                "among_taxon_rv": a,
                "delta_among_minus_within": a - w,
            })
    return within, among, pd.DataFrame(rows)


def matrix_alignment(within: pd.DataFrame, among: pd.DataFrame) -> float:
    upper = np.triu_indices(len(CORE), 1)
    return float(spearmanr(within.to_numpy()[upper], among.to_numpy()[upper]).statistic)


def qap_alignment(within: pd.DataFrame, among: pd.DataFrame, permutations: int, rng: np.random.Generator) -> dict[str, float | int]:
    upper = np.triu_indices(len(CORE), 1)
    w = within.to_numpy(float)
    a = among.to_numpy(float)
    observed = float(spearmanr(w[upper], a[upper]).statistic)
    sims = np.empty(permutations, dtype=float)
    for i in range(permutations):
        order = rng.permutation(len(CORE))
        ap = a[np.ix_(order, order)]
        sims[i] = float(spearmanr(w[upper], ap[upper]).statistic)
    return {
        "rho": observed,
        "qap_p_one_sided": float((np.sum(sims >= observed - 1e-15) + 1) / (permutations + 1)),
        "permutations": permutations,
        "null_mean": float(np.mean(sims)),
        "null_low95": float(np.quantile(sims, 0.025)),
        "null_high95": float(np.quantile(sims, 0.975)),
    }


def bootstrap_alignment(frame: pd.DataFrame, reps: int, seed: int) -> pd.DataFrame:
    taxa = np.array(sorted(frame.taxon_name.unique()), dtype=object)
    rng = np.random.default_rng(seed)
    rows = []
    fmap = frame.attrs["feature_map"]
    for rep in range(reps):
        sampled = rng.choice(taxa, size=len(taxa), replace=True)
        pieces = []
        for j, taxon in enumerate(sampled):
            piece = frame[frame.taxon_name.eq(taxon)].copy()
            piece["taxon_name"] = f"boot_{j:04d}"
            pieces.append(piece)
        boot = pd.concat(pieces, ignore_index=True)
        boot.attrs["feature_map"] = fmap
        w, a, pairs = matrices_from_common(boot)
        rows.append({
            "replicate": rep,
            "matrix_alignment_rho": matrix_alignment(w, a),
            "median_within_rv": float(pairs.within_taxon_rv.median()),
            "median_among_rv": float(pairs.among_taxon_rv.median()),
            "relations_stronger_among": int((pairs.among_taxon_rv > pairs.within_taxon_rv).sum()),
        })
    return pd.DataFrame(rows)


def module_cohesion(matrix: pd.DataFrame, permutations: int, rng: np.random.Generator) -> dict[str, float | int]:
    pairs = []
    for i, left in enumerate(CORE):
        for right in CORE[i + 1 :]:
            pairs.append((left, right, float(matrix.loc[left, right])))
    within_vals = [v for l, r, v in pairs if MODULES[l] == MODULES[r]]
    between_vals = [v for l, r, v in pairs if MODULES[l] != MODULES[r]]
    observed = float(np.mean(within_vals) - np.mean(between_vals))

    labels = np.array([MODULES[c] for c in CORE], dtype=object)
    sims = np.empty(permutations, dtype=float)
    for b in range(permutations):
        perm = rng.permutation(labels)
        mapping = {c: perm[i] for i, c in enumerate(CORE)}
        wi = [v for l, r, v in pairs if mapping[l] == mapping[r]]
        be = [v for l, r, v in pairs if mapping[l] != mapping[r]]
        sims[b] = float(np.mean(wi) - np.mean(be))
    return {
        "mean_within_module_rv": float(np.mean(within_vals)),
        "mean_between_module_rv": float(np.mean(between_vals)),
        "within_minus_between": observed,
        "permutation_p_one_sided": float((np.sum(sims >= observed - 1e-15) + 1) / (permutations + 1)),
        "permutations": permutations,
    }


def block_signatures(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path)
    d = d[d.construct_id.isin(CORE)].copy()
    rows = []
    for construct in CORE:
        part = d[d.construct_id.eq(construct)]
        block_values = {}
        for block, predictors in ENV_BLOCKS.items():
            vals = part[part.predictor.isin(predictors)].effect_magnitude.to_numpy(float)
            if len(vals) != len(predictors):
                raise ValueError(f"missing predictors for {construct} {block}")
            block_values[block] = float(np.sqrt(np.mean(vals ** 2)))
        norm = math.sqrt(sum(v * v for v in block_values.values()))
        row = {"construct_id": construct}
        for block, value in block_values.items():
            row[block] = value
            row[f"norm_{block}"] = value / norm if norm > 0 else np.nan
        row["strongest_block"] = max(block_values, key=block_values.get)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    t, _ = axis.load_data(args.traits, args.environment)
    common, endpoints = common_complete_features(t, args.minimum_complete_observations_per_taxon)
    n_obs = int(len(common))
    n_taxa = int(common.taxon_name.nunique())
    if n_obs < 100 or n_taxa < 20:
        raise SystemExit(f"complete-18 common cohort too small: {n_obs} observations / {n_taxa} taxa")

    within, among, pairs = matrices_from_common(common)
    qap = qap_alignment(within, among, args.permutations, np.random.default_rng(args.seed + 1))
    boot = bootstrap_alignment(common, args.bootstrap_replicates, args.seed + 2)
    boot_ci = {
        "rho_median": float(boot.matrix_alignment_rho.median()),
        "rho_low95": float(boot.matrix_alignment_rho.quantile(0.025)),
        "rho_high95": float(boot.matrix_alignment_rho.quantile(0.975)),
        "probability_rho_positive": float((boot.matrix_alignment_rho > 0).mean()),
    }

    module_within = module_cohesion(within, args.permutations, np.random.default_rng(args.seed + 3))
    module_among = module_cohesion(among, args.permutations, np.random.default_rng(args.seed + 4))

    sig_a = block_signatures(args.axis_among)
    sig_w = block_signatures(args.axis_within)
    sig = sig_a.merge(sig_w, on="construct_id", suffixes=("_among", "_within"))
    same_block = int((sig.strongest_block_among == sig.strongest_block_within).sum())
    block_cols_a = [f"norm_{b}_among" for b in ENV_BLOCKS]
    block_cols_w = [f"norm_{b}_within" for b in ENV_BLOCKS]
    flat_rho = float(spearmanr(sig[block_cols_a].to_numpy().ravel(), sig[block_cols_w].to_numpy().ravel()).statistic)
    flat_cos = float(
        np.dot(sig[block_cols_a].to_numpy().ravel(), sig[block_cols_w].to_numpy().ravel()) /
        (np.linalg.norm(sig[block_cols_a].to_numpy()) * np.linalg.norm(sig[block_cols_w].to_numpy()))
    )

    within.to_csv(args.out_dir / "complete18_construct_integration_within.csv")
    among.to_csv(args.out_dir / "complete18_construct_integration_among.csv")
    pairs.to_csv(args.out_dir / "complete18_construct_pairwise.csv", index=False)
    boot.to_csv(args.out_dir / "complete18_taxon_bootstrap.csv", index=False)
    sig.to_csv(args.out_dir / "six_block_environment_signatures.csv", index=False)

    report = {
        "analysis_id": "ch1_v3_construct_scale_upgrade_20260910",
        "claim_boundary": "construct-level synthesis used by the current manuscript; frozen v2 endpoint conclusions and multiplicity families unchanged",
        "complete18_endpoints": endpoints,
        "common_cohort": {
            "observations": n_obs,
            "taxa": n_taxa,
            "minimum_complete_observations_per_taxon": args.minimum_complete_observations_per_taxon,
        },
        "common_cohort_matrix_alignment": qap,
        "taxon_bootstrap": {"replicates": args.bootstrap_replicates, **boot_ci},
        "median_within_rv": float(pairs.within_taxon_rv.median()),
        "median_among_rv": float(pairs.among_taxon_rv.median()),
        "relations_stronger_among": int((pairs.among_taxon_rv > pairs.within_taxon_rv).sum()),
        "module_cohesion": {"within_taxon": module_within, "among_taxon": module_among},
        "six_environment_block_signature": {
            "blocks": ENV_BLOCKS,
            "same_strongest_block_constructs": same_block,
            "constructs_compared": len(CORE),
            "flattened_normalized_signature_spearman": flat_rho,
            "flattened_normalized_signature_cosine": flat_cos,
            "interpretation": "descriptive block-level comparison; blocks reduce but do not eliminate predictor correlation",
        },
    }
    (args.out_dir / "construct_scale_upgrade_report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare standardized within- vs among-species capitulum hypervolume geometry.

Input is the signed 19-endpoint correlation matrices produced by the v3
whole-capitulum multiscale analysis.  The analysis treats each correlation
matrix as a zero-centred Gaussian ellipsoid in the same standardized 19-D
endpoint space.  Therefore volume is a *correlation-geometry* volume, not raw
trait-space volume in biological units.

Outputs quantify:
  * generalized-variance / ellipsoid volume ratio among vs within;
  * Bhattacharyya overlap and Hellinger separation of the two ellipsoids;
  * generalized eigenaxes solving Sigma_among v = lambda Sigma_within v,
    where lambda > 1 denotes directions relatively expanded among species and
    lambda < 1 denotes directions relatively contracted among species;
  * loading summaries for biological interpretation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--correlations", required=True, type=Path)
    p.add_argument("--contract", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def reconstruct(df: pd.DataFrame, scope: str, scale: str) -> pd.DataFrame:
    x = df[(df["scope"] == scope) & (df["scale"] == scale)]
    names = sorted(set(x["left"]).union(x["right"]))
    out = pd.DataFrame(np.eye(len(names)), index=names, columns=names, dtype=float)
    for row in x.itertuples(index=False):
        out.loc[row.left, row.right] = float(row.value)
        out.loc[row.right, row.left] = float(row.value)
    return out


def generalized_axes(within: np.ndarray, among: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Symmetric whitening avoids a SciPy dependency.
    ew, uw = np.linalg.eigh(within)
    if np.min(ew) <= 0:
        raise ValueError(f"within matrix not positive definite: min eigenvalue={np.min(ew)}")
    invsqrt = uw @ np.diag(1.0 / np.sqrt(ew)) @ uw.T
    whitened = invsqrt @ among @ invsqrt
    vals, q = np.linalg.eigh((whitened + whitened.T) / 2.0)
    vecs = invsqrt @ q
    # Euclidean normalization for readable coefficient summaries.
    vecs = vecs / np.linalg.norm(vecs, axis=0, keepdims=True)
    return vals, vecs


def ellipsoid_metrics(within: np.ndarray, among: np.ndarray) -> dict[str, float]:
    sw = np.linalg.slogdet(within)
    sa = np.linalg.slogdet(among)
    if sw[0] <= 0 or sa[0] <= 0:
        raise ValueError("Correlation matrix determinant must be positive")
    mid = (within + among) / 2.0
    sm = np.linalg.slogdet(mid)
    if sm[0] <= 0:
        raise ValueError("Midpoint covariance determinant must be positive")
    log_volume_ratio = 0.5 * (sa[1] - sw[1])
    db = 0.5 * (sm[1] - 0.5 * (sw[1] + sa[1]))
    bc = float(np.exp(-db))
    return {
        "within_log_generalized_variance": float(sw[1]),
        "among_log_generalized_variance": float(sa[1]),
        "among_vs_within_log_ellipsoid_volume_ratio": float(log_volume_ratio),
        "among_vs_within_ellipsoid_volume_ratio": float(np.exp(log_volume_ratio)),
        "bhattacharyya_distance": float(db),
        "bhattacharyya_coefficient": bc,
        "hellinger_distance": float(np.sqrt(max(0.0, 1.0 - bc))),
    }


def main() -> int:
    args = parse_args()
    corr = pd.read_csv(args.correlations)
    corr = corr[corr["matrix_type"].eq("signed_endpoint_correlation")].copy()
    contract = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    module = contract.set_index("endpoint_id")["module"].to_dict()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    summaries = []
    axes_rows = []
    loading_rows = []
    axis_store: dict[str, tuple[list[str], np.ndarray, np.ndarray]] = {}

    for scope in sorted(corr["scope"].unique()):
        wdf = reconstruct(corr, scope, "within_taxon")
        adf = reconstruct(corr, scope, "among_taxon").loc[wdf.index, wdf.columns]
        names = list(wdf.index)
        w = wdf.to_numpy(float)
        a = adf.to_numpy(float)
        vals, vecs = generalized_axes(w, a)
        met = ellipsoid_metrics(w, a)
        met.update({
            "scope": scope,
            "n_endpoint_dimensions": len(names),
            "within_min_eigenvalue": float(np.linalg.eigvalsh(w).min()),
            "among_min_eigenvalue": float(np.linalg.eigvalsh(a).min()),
            "n_axes_relatively_contracted_among_lambda_lt_1": int(np.sum(vals < 1.0)),
            "n_axes_relatively_expanded_among_lambda_gt_1": int(np.sum(vals > 1.0)),
            "generalized_eigenvalue_geometric_mean": float(np.exp(np.mean(np.log(vals)))),
            "generalized_eigenvalue_median": float(np.median(vals)),
            "strongest_contraction_lambda": float(vals[0]),
            "strongest_expansion_lambda": float(vals[-1]),
        })
        summaries.append(met)
        axis_store[scope] = (names, vals, vecs)
        for i, lam in enumerate(vals):
            axes_rows.append({
                "scope": scope,
                "axis_rank_ascending": i + 1,
                "lambda_among_over_within": float(lam),
                "relative_state": "expanded_among" if lam > 1 else "contracted_among",
            })
            v = vecs[:, i]
            order = np.argsort(np.abs(v))[::-1]
            for rank, j in enumerate(order, start=1):
                loading_rows.append({
                    "scope": scope,
                    "axis_rank_ascending": i + 1,
                    "lambda_among_over_within": float(lam),
                    "loading_rank": rank,
                    "endpoint_id": names[j],
                    "module": module.get(names[j], "unknown"),
                    "coefficient": float(v[j]),
                    "abs_coefficient": float(abs(v[j])),
                })

    # Match min5 axes to min2 axes by absolute Euclidean coefficient cosine.
    matches = []
    scopes = sorted(axis_store)
    if "complete19_min5" in axis_store and "complete19_min2" in axis_store:
        n5, v5, u5 = axis_store["complete19_min5"]
        n2, v2, u2 = axis_store["complete19_min2"]
        if n5 == n2:
            sims = np.abs(u2.T @ u5)
            for i in range(len(v5)):
                j = int(np.argmax(sims[:, i]))
                matches.append({
                    "min5_axis_rank_ascending": i + 1,
                    "min5_lambda": float(v5[i]),
                    "best_min2_axis_rank_ascending": j + 1,
                    "min2_lambda": float(v2[j]),
                    "absolute_coefficient_cosine": float(sims[j, i]),
                })

    pd.DataFrame(summaries).to_csv(out / "capitulum_hypervolume_summary.csv", index=False)
    pd.DataFrame(axes_rows).to_csv(out / "capitulum_hypervolume_generalized_axes.csv", index=False)
    pd.DataFrame(loading_rows).to_csv(out / "capitulum_hypervolume_axis_loadings.csv", index=False)
    pd.DataFrame(matches).to_csv(out / "capitulum_hypervolume_axis_stability.csv", index=False)

    report = {
        "analysis_id": "ch1_capitulum_standardized_hypervolume_geometry_v1",
        "input": str(args.correlations),
        "definition": "zero-centred Gaussian ellipsoid on the standardized 19-endpoint correlation geometry",
        "claim_boundary": (
            "This is not raw-unit trait-space volume and not a kernel-density ecological hypervolume. "
            "It quantifies covariance-geometry compression/expansion and overlap between within- and among-species standardized phenotype organization."
        ),
        "summary": summaries,
    }
    (out / "capitulum_hypervolume_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

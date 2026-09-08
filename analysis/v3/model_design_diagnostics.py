"""Exposure-only diagnostics on an explicitly supplied realized model cohort.

No response values or regressions are read. Pooled within diagnostics do not prove
that each taxon slope is identifiable, so both pooled and taxon-level ranks appear.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .environment_matrix import vif_table


def _summary(matrix: np.ndarray, names: list[str], weights: np.ndarray) -> dict:
    frame = pd.DataFrame(matrix, columns=names)
    frame["weight"] = weights
    vif, result = vif_table(frame, names, "weight")
    result["vif"] = vif.to_dict(orient="records")
    result["constant_variables"] = [name for i, name in enumerate(names) if np.ptp(matrix[:, i]) == 0]
    # Same complete cohort, not pairwise-complete correlations with changing n.
    center = np.average(matrix, axis=0, weights=weights)
    covariance = ((matrix - center) * weights[:, None]).T @ (matrix - center) / weights.sum()
    sd = np.sqrt(np.maximum(np.diag(covariance), 0))
    denom = sd[:, None] * sd[None, :]
    correlation = np.divide(covariance, denom, out=np.full_like(covariance, np.nan), where=denom > 0)
    result["correlation"] = [[float(v) if np.isfinite(v) else None for v in row] for row in correlation]
    result["variables"] = names
    result["full_predictor_rank"] = result["matrix_rank"] == len(names)
    return result


def _weighted_center(matrix: np.ndarray, weights: np.ndarray) -> np.ndarray:
    offset = matrix - matrix[:1]
    return offset - np.average(offset, axis=0, weights=weights)


def _project_out(matrix: np.ndarray, nuisance: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, int]:
    if nuisance.shape[1] == 0:
        return matrix.copy(), 0
    has_intercept = bool(np.any((np.ptp(nuisance, axis=0) == 0) & (nuisance[0] != 0)))
    if has_intercept:
        # Remove the intercept stably first. Alias tolerance is relative to
        # informative variation, not the arbitrary origin of a predictor unit.
        matrix = _weighted_center(matrix, weights)
        nuisance = _weighted_center(nuisance, weights)
    sqrt_w = np.sqrt(weights)
    weighted_n = nuisance * sqrt_w[:, None]
    # Scale columns before rank decisions; preserve zero columns as zero.
    lengths = np.linalg.norm(weighted_n, axis=0)
    retained = lengths > 0
    weighted_n = weighted_n[:, retained] / lengths[retained]
    u, s, _ = np.linalg.svd(weighted_n, full_matrices=False)
    tol = max(weighted_n.shape) * np.finfo(float).eps * (s[0] if len(s) else 0)
    rank = int((s > tol).sum())
    q = u[:, :rank]
    resid = (matrix * sqrt_w[:, None] - q @ (q.T @ (matrix * sqrt_w[:, None]))) / sqrt_w[:, None]
    # Exact aliases can leave floating-point dust; never renormalize it into a slope.
    tiny = np.linalg.norm(resid * sqrt_w[:, None], axis=0) <= 1e-10 * np.maximum(
        np.linalg.norm(matrix * sqrt_w[:, None], axis=0), np.finfo(float).tiny)
    resid[:, tiny] = 0
    return resid, rank + int(has_intercept)


def diagnose(frame: pd.DataFrame, predictors: list[str], nuisance: list[str], *,
             weight_column: str, cohort_id: str) -> dict:
    """Use actual supplied fit weights; report, never silently repair, missing rows.

    The caller supplies finalized calendar/imaging/spatial design columns. Among
    summaries are equal-observation taxon means, including means of each spatial
    basis column; they are not values sampled at a taxon centroid. Among weights
    here are equal-taxon for exposure diagnostics, not fitted variance weights.
    """
    if frame.empty:
        raise ValueError("Empty realized model cohort")
    if not cohort_id or not predictors or len(set(predictors + nuisance)) != len(predictors + nuisance):
        raise ValueError("Named cohort and distinct predictor/nuisance columns required")
    columns = ["obs_id", "accepted_key", weight_column, *predictors, *nuisance]
    if not set(columns) <= set(frame):
        raise ValueError("Required realized-cohort/design columns absent")
    data = frame[columns].copy()
    if data[["obs_id", "accepted_key"]].isna().any().any() or data["obs_id"].astype(str).duplicated().any():
        raise ValueError("Missing identity or duplicated observation")
    numeric = data[[weight_column, *predictors, *nuisance]].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all() or not numeric[weight_column].gt(0).all():
        raise ValueError("Nonfinite design/weights: freeze an explicit missingness cohort first")
    data[numeric.columns] = numeric
    weights = data[weight_column].to_numpy(dtype=float)
    matrix = data[predictors].to_numpy(dtype=float)
    z = data[nuisance].to_numpy(dtype=float)
    within_x = matrix.copy()
    within_z = z.copy()
    taxon_ranks = []
    # Exposure centering uses the exact supplied observation weights.
    for _, indices in data.groupby("accepted_key", sort=True).indices.items():
        local_w = weights[indices]
        within_x[indices] = _weighted_center(matrix[indices], local_w)
        if nuisance:
            within_z[indices] = _weighted_center(z[indices], local_w)
        local, _ = _project_out(within_x[indices], within_z[indices], local_w)
        taxon_ranks.append(_summary(local, predictors, local_w)["matrix_rank"])
    within_resid, nuisance_rank = _project_out(within_x, within_z, weights)
    means = data.groupby("accepted_key", sort=True)[predictors + nuisance].mean()
    among_x = means[predictors].to_numpy(dtype=float)
    among_z = np.column_stack([np.ones(len(means)), means[nuisance].to_numpy(dtype=float)])
    among_resid, among_nuisance_rank = _project_out(among_x, among_z, np.ones(len(means)))
    return {
        "status": "REALIZED_COHORT_ENVIRONMENT_DIAGNOSTICS_REPORTED",
        "cohort_id": cohort_id, "observations": len(data), "taxa": len(means),
        "weight_column": weight_column, "nuisance_columns": nuisance,
        "raw": _summary(matrix, predictors, weights),
        "within": _summary(within_x, predictors, weights),
        "within_nuisance_residualized": _summary(within_resid, predictors, weights),
        "among": _summary(among_x, predictors, np.ones(len(means))),
        "among_nuisance_residualized": _summary(among_resid, predictors, np.ones(len(means))),
        "within_nuisance_rank": nuisance_rank, "among_nuisance_rank": among_nuisance_rank,
        "taxa_with_full_local_predictor_rank": sum(rank == len(predictors) for rank in taxon_ranks),
        "taxa_without_full_local_predictor_rank": sum(rank < len(predictors) for rank in taxon_ranks),
        "trait_values_read": 0, "ecological_models_executed": 0,
        "ecological_fitting_authorized": False,
        "limits": [
            "This diagnoses exposure identification, not confounder sufficiency or a fitted ecological effect.",
            "Local residualized rank is a conservative per-taxon check; it is not the rank of the complete shared-nuisance heterogeneous-slope model.",
            "Final inferential joint-design and covariance checks remain required.",
            "Among exposure diagnostics use equal taxon weights; any fitted variance weights require an additional saved weighted diagnostic.",
        ],
    }

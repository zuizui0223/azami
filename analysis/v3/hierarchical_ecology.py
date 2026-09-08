"""Core estimators for the frozen Chapter 1 v3 hierarchical ecology model.

This module implements only analysis mechanics fixed before v3 trait outcomes:
common taxon+spatial nuisance residualization, supported taxon-specific
standardized conditional slopes, and normal-normal REML partial pooling. It does
not choose predictors, endpoints or cohorts and does not inspect v2 results.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree


SPATIAL_TERMS = ("x", "y", "z", "xy", "xz", "yz", "x2_minus_y2", "three_z2_minus_1")


def spherical_basis(latitude: Iterable[float], longitude: Iterable[float]) -> np.ndarray:
    lat = np.radians(np.asarray(latitude, dtype=float))
    lon = np.radians(np.asarray(longitude, dtype=float))
    if lat.shape != lon.shape or lat.ndim != 1:
        raise ValueError("Latitude/longitude must be matching vectors")
    if not np.all(np.isfinite(lat)) or not np.all(np.isfinite(lon)):
        raise ValueError("Spatial coordinates must be finite")
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack((x, y, z, x*y, x*z, y*z, x*x-y*y, 3*z*z-1.0))


def _demean_by_group(matrix: np.ndarray, groups: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    vector = matrix.ndim == 1
    if vector:
        matrix = matrix[:, None]
    frame = pd.DataFrame(matrix)
    labels = pd.Series(groups, index=frame.index)
    out = (frame - frame.groupby(labels, sort=False).transform("mean")).to_numpy(float)
    return out[:, 0] if vector else out


def common_spatial_residualize(
    response: np.ndarray,
    predictors: np.ndarray,
    taxa: np.ndarray,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """FWL residualization for taxon intercepts plus common spherical basis.

    Taxon demeaning removes taxon intercepts. The demeaned spherical basis is
    then projected out globally from both response and predictor columns. This
    is algebraically the nuisance-removal part of a model containing taxon fixed
    intercepts and common spatial-basis coefficients; it does not estimate eight
    spatial coefficients separately within every taxon.
    """
    y = np.asarray(response, dtype=float)
    x = np.asarray(predictors, dtype=float)
    g = np.asarray(taxa)
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if y.ndim != 1 or x.ndim != 2 or len(y) != len(x) or len(g) != len(y):
        raise ValueError("Response, predictor matrix and taxa have incompatible shapes")
    if len(lat) != len(y) or len(lon) != len(y):
        raise ValueError("Coordinate length mismatch")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(x)):
        raise ValueError("FWL inputs must be finite after cohort filtering")

    basis = spherical_basis(lat, lon)
    y0 = _demean_by_group(y, g)
    x0 = _demean_by_group(x, g)
    b0 = _demean_by_group(basis, g)
    keep = np.std(b0, axis=0, ddof=0) > 1e-12
    b = b0[:, keep]
    if b.shape[1]:
        nuisance_inverse = np.linalg.pinv(b)
        y_resid = y0 - b @ (nuisance_inverse @ y0)
        x_resid = x0 - b @ (nuisance_inverse @ x0)
        rank = int(np.linalg.matrix_rank(b))
    else:
        y_resid, x_resid, rank = y0, x0, 0
    return y_resid, x_resid, {
        "taxa": int(pd.Series(g).nunique()),
        "spatial_basis_columns_retained": int(b.shape[1]),
        "spatial_basis_rank": rank,
        "spatial_terms_retained": [SPATIAL_TERMS[i] for i, flag in enumerate(keep) if flag],
    }


def _standardize(values: np.ndarray) -> tuple[np.ndarray, float]:
    values = np.asarray(values, dtype=float)
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else float("nan")
    if not np.isfinite(sd) or sd <= 1e-12:
        raise ValueError("No estimable within-taxon variation")
    return (values - float(np.mean(values))) / sd, sd


@dataclass(frozen=True)
class TaxonSlope:
    taxon: str
    n_observations: int
    n_cells: int
    predictor_names: tuple[str, ...]
    beta: tuple[float, ...]
    standard_error: tuple[float, ...]
    covariance: tuple[tuple[float, ...], ...]
    residual_df: int


def estimate_taxon_slopes(
    response_residual: np.ndarray,
    predictor_residuals: np.ndarray,
    taxa: np.ndarray,
    cells: np.ndarray,
    predictor_names: list[str],
    minimum_observations: int = 10,
    minimum_cells: int = 4,
    minimum_residual_df: int = 3,
) -> tuple[list[TaxonSlope], pd.DataFrame]:
    """Estimate supported within-taxon standardized conditional slopes.

    All predictors and the response are standardized within taxon after common
    nuisance residualization. Taxa failing source support, predictor rank or
    residual-df rules are returned in the support ledger but not as unstable
    slope estimates.
    """
    y = np.asarray(response_residual, dtype=float)
    x = np.asarray(predictor_residuals, dtype=float)
    g = np.asarray(taxa).astype(str)
    c = np.asarray(cells).astype(str)
    p = len(predictor_names)
    if x.ndim != 2 or x.shape[1] != p or len(y) != len(x) or len(g) != len(y) or len(c) != len(y):
        raise ValueError("Taxon-slope inputs have incompatible dimensions")
    records: list[TaxonSlope] = []
    ledger = []
    for taxon in sorted(set(g)):
        index = np.flatnonzero(g == taxon)
        n = len(index)
        n_cells = len(set(c[index]))
        reasons = []
        if n < minimum_observations:
            reasons.append("below_minimum_observations")
        if n_cells < minimum_cells:
            reasons.append("below_minimum_cells")
        if n < p + minimum_residual_df:
            reasons.append("insufficient_residual_df")
        yy = y[index]
        xx = x[index]
        standardized_x = []
        try:
            y_std, _ = _standardize(yy)
        except ValueError:
            y_std = None
            reasons.append("response_no_variation")
        for column in range(p):
            try:
                z, _ = _standardize(xx[:, column])
                standardized_x.append(z)
            except ValueError:
                standardized_x.append(None)
                reasons.append(f"predictor_no_variation:{predictor_names[column]}")
        if all(z is not None for z in standardized_x):
            design = np.column_stack(standardized_x)
            rank = int(np.linalg.matrix_rank(design))
            if rank < p:
                reasons.append("predictor_rank_deficient")
        else:
            design, rank = None, 0
        residual_df = n - p
        if not reasons and y_std is not None and design is not None:
            xtx_inverse = np.linalg.inv(design.T @ design)
            beta = xtx_inverse @ design.T @ y_std
            residual = y_std - design @ beta
            sigma2 = float(np.dot(residual, residual) / residual_df)
            covariance = xtx_inverse * sigma2
            standard_error = np.sqrt(np.diag(covariance))
            records.append(TaxonSlope(
                taxon=taxon,
                n_observations=n,
                n_cells=n_cells,
                predictor_names=tuple(predictor_names),
                beta=tuple(map(float, beta)),
                standard_error=tuple(map(float, standard_error)),
                covariance=tuple(tuple(map(float, row)) for row in covariance),
                residual_df=int(residual_df),
            ))
        ledger.append({
            "taxon": taxon,
            "n_observations": int(n),
            "n_cells": int(n_cells),
            "n_predictors": int(p),
            "predictor_rank": int(rank),
            "residual_df": int(residual_df),
            "slope_estimable": not reasons,
            "reasons": ";".join(sorted(set(reasons))),
        })
    return records, pd.DataFrame(ledger)


@dataclass(frozen=True)
class RandomEffectsResult:
    n_taxa: int
    hypermean: float
    hypermean_se: float
    confidence_low_95: float
    confidence_high_95: float
    tau2: float
    tau: float
    reml_objective: float
    shrinkage: pd.DataFrame


def random_effects_reml(estimates: Iterable[float], variances: Iterable[float], taxa: Iterable[str]) -> RandomEffectsResult:
    y = np.asarray(list(estimates), dtype=float)
    v = np.asarray(list(variances), dtype=float)
    labels = np.asarray(list(taxa)).astype(str)
    if len(y) != len(v) or len(y) != len(labels) or len(y) < 2:
        raise ValueError("Random-effects meta-analysis needs at least two aligned taxa")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("Taxon slope estimates and variances must be finite with positive variances")

    between = float(np.var(y, ddof=1))
    upper = max(1.0, between * 100.0, float(np.max(v)) * 100.0)

    def objective(tau2: float) -> float:
        vv = v + tau2
        w = 1.0 / vv
        sw = float(w.sum())
        mu = float(np.dot(w, y) / sw)
        q = float(np.dot(w, (y-mu)**2))
        return 0.5 * (float(np.log(vv).sum()) + math.log(sw) + q)

    fit = minimize_scalar(objective, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-12})
    if not fit.success or not np.isfinite(fit.fun):
        raise RuntimeError("REML tau2 optimization failed")
    tau2 = max(0.0, float(fit.x))
    # If the optimum is numerically indistinguishable from the boundary, use zero.
    if objective(0.0) <= objective(tau2) + 1e-10:
        tau2 = 0.0
    w = 1.0 / (v + tau2)
    sw = float(w.sum())
    mu = float(np.dot(w, y) / sw)
    mu_se = math.sqrt(1.0 / sw)

    if tau2 <= 0:
        posterior_mean = np.full(len(y), mu)
        posterior_var = np.zeros(len(y))
    else:
        posterior_var = 1.0 / (1.0/v + 1.0/tau2)
        posterior_mean = posterior_var * (y/v + mu/tau2)
    shrinkage = pd.DataFrame({
        "taxon": labels,
        "raw_slope": y,
        "sampling_variance": v,
        "sampling_se": np.sqrt(v),
        "shrunken_slope": posterior_mean,
        "shrunken_se": np.sqrt(posterior_var),
        "precision_weight": w,
    }).sort_values("taxon", kind="mergesort").reset_index(drop=True)
    return RandomEffectsResult(
        n_taxa=len(y),
        hypermean=mu,
        hypermean_se=mu_se,
        confidence_low_95=mu - 1.959963984540054*mu_se,
        confidence_high_95=mu + 1.959963984540054*mu_se,
        tau2=tau2,
        tau=math.sqrt(tau2),
        reml_objective=float(objective(tau2)),
        shrinkage=shrinkage,
    )


def morans_i(values: np.ndarray, latitude: np.ndarray, longitude: np.ndarray, k: int = 8) -> float:
    values = np.asarray(values, dtype=float)
    if len(values) <= k or not np.all(np.isfinite(values)):
        return float("nan")
    xyz = spherical_basis(latitude, longitude)[:, :3]
    neighbours = cKDTree(xyz).query(xyz, k=k+1)[1][:, 1:]
    centred = values - float(np.mean(values))
    denominator = float(np.dot(centred, centred))
    if denominator <= 0:
        return float("nan")
    numerator = float(np.sum(centred[:, None] * centred[neighbours]))
    return numerator / (k * denominator)

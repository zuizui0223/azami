"""Core estimators for the frozen Chapter 1 v3 hierarchical ecology model.

This module contains pre-outcome numerical building blocks, not an authorized
ecological runner. Common-slope residualization must not be followed by separate
taxon fits: heterogeneous slopes require the joint interaction design below.
Arbitrary shared nuisance columns are supported by the joint solver; their exact
calendar/imaging construction and dependence-aware pooling remain unfinished. It does
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

    This reduction is valid for a COMMON environmental slope. Splitting these
    residuals by taxon afterward is not FWL for heterogeneous slopes: the
    taxon-by-predictor columns must enter the joint model before projection.
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

    This is a separate-group OLS helper, NOT the spatial hierarchical estimator.
    Do not pass global common-spatial residuals to it and interpret the results
    as jointly adjusted taxon slopes. Taxa failing source support, predictor rank or
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
        if n < p + 1 + minimum_residual_df:
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
        # Centering estimates one group intercept even without an explicit column.
        residual_df = n - p - 1
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
class JointSpatialSlopes:
    """Raw-unit joint slopes and full HC3 covariance (taxon-major ordering).

    HC3 assumes independent observational errors. It is a numerical reference,
    not a substitute for shared-image grouping or spatial block uncertainty.
    Off-diagonal taxon covariance must not be discarded when pooling.
    """

    taxa: tuple[str, ...]
    slopes: np.ndarray
    covariance: np.ndarray
    residuals: np.ndarray
    leverage: np.ndarray
    residual_df: int
    spatial_rank: int


def joint_spatial_taxon_slopes(response, predictors, taxa, latitude, longitude) -> JointSpatialSlopes:
    """Historical eight-term spatial wrapper around the general joint solver."""
    result=joint_nuisance_taxon_slopes(response,predictors,taxa,spherical_basis(latitude,longitude))
    return JointSpatialSlopes(result.taxa,result.slopes,result.covariance,result.residuals,
                              result.leverage,result.residual_df,result.nuisance_rank)


@dataclass(frozen=True)
class JointNuisanceSlopes:
    """Joint raw-unit slopes with all shared nuisance terms fitted together.

    HC3 remains an independent-error reference, not spatial/component inference.
    """
    taxa: tuple[str,...]
    slopes: np.ndarray
    covariance: np.ndarray
    residuals: np.ndarray
    leverage: np.ndarray
    residual_df: int
    nuisance_rank: int


def joint_nuisance_taxon_slopes(response,predictors,taxa,nuisance) -> JointNuisanceSlopes:
    """Fit taxon intercepts/slopes and arbitrary shared nuisance effects JOINTLY.

    Use block FWL: first remove each taxon's intercept and predictor columns
    from the common nuisance basis, solve the remaining small-column system,
    then recover all taxon slopes. This avoids an observations-by-all-taxon-
    interactions dense matrix. No response scaling, cohort selection, scientific
    support threshold, timing/imaging coding or pooling is performed here.
    The caller must supply the predeclared aligned observational units.
    """
    y, x = np.asarray(response, float), np.asarray(predictors, float)
    groups = np.asarray(taxa)
    if (y.ndim != 1 or x.ndim != 2 or groups.ndim != 1
            or len(x) != len(y) or len(groups) != len(y) or x.shape[1] == 0):
        raise ValueError("Joint slope inputs have incompatible dimensions")
    if (not np.isfinite(y).all() or not np.isfinite(x).all()
            or pd.isna(groups).any()):
        raise ValueError("Joint slope inputs must be finite with nonmissing taxa")
    groups = groups.astype(str)
    if np.any(groups == ""):
        raise ValueError("Missing taxon label")
    basis=np.asarray(nuisance,float)
    if basis.ndim!=2 or len(basis)!=len(y) or not np.isfinite(basis).all():
        raise ValueError("Nuisance matrix must be finite and aligned")
    labels, p = tuple(sorted(set(groups))), x.shape[1]
    y_perp, b_perp = np.empty_like(y), np.empty_like(basis)
    b_centered = np.empty_like(basis)
    blocks = []
    for label in labels:
        ix = np.flatnonzero(groups == label)
        # Subtract an observed anchor before the mean: exactly constant input
        # then maps to exact zero instead of being admitted as rounding noise.
        xx = x[ix] - x[ix[0]]
        xx -= xx.mean(axis=0)
        yy = y[ix] - y[ix[0]]
        yy -= yy.mean()
        bb = basis[ix] - basis[ix[0]]
        bb -= bb.mean(axis=0)
        # A repeated floating-point coordinate can leave mean-subtraction dust.
        # Exact within-block constants belong to the taxon intercept, not space.
        bb[:, np.ptp(basis[ix], axis=0) == 0] = 0.0
        if len(ix) <= p + 1 or np.linalg.matrix_rank(xx) < p:
            raise ValueError(f"Taxon predictor block not estimable: {label}")
        inverse_x = np.linalg.pinv(xx)
        a = inverse_x @ bb
        y_perp[ix] = yy - xx @ (inverse_x @ yy)
        b_perp[ix] = bb - xx @ a
        b_centered[ix] = bb
        blocks.append((ix, xx, yy, bb, inverse_x, a))

    # Scale the nuisance columns before solving; exact redundant spatial terms
    # are allowed, but confounding of an environmental slope with space is not.
    scales = np.linalg.norm(b_centered, axis=0)
    inactive = scales <= 1e-12
    scales[inactive] = 1.0
    bc, bp = b_centered / scales, b_perp / scales
    bc[:, inactive] = 0.0
    bp[:, inactive] = 0.0
    tolerance = (np.linalg.norm(bc, ord=2) if bc.shape[1] else 0.0) * max(bc.shape) * np.finfo(float).eps
    u, singular, vt = np.linalg.svd(bp, full_matrices=False)
    retained = singular > tolerance
    rank = int(retained.sum())
    centered_rank=int(np.linalg.matrix_rank(bc,tol=tolerance)) if bc.shape[1] else 0
    if rank < centered_rank:
        raise ValueError("Taxon slopes aliased with the common nuisance basis")
    residual_df = len(y) - len(labels) * (p + 1) - rank
    if residual_df <= 0:
        raise ValueError("Joint model has no residual degrees of freedom")
    # Use the SAME threshold in rank accounting and the pseudoinverse.
    inverse_b = (vt[retained].T / singular[retained]) @ u[:, retained].T
    gamma = inverse_b @ y_perp
    slopes, residual, leverage, scaled_a = [], np.empty_like(y), np.empty_like(y), []
    for ix, xx, yy, bb, inverse_x, a in blocks:
        beta = inverse_x @ (yy - (bb / scales) @ gamma)
        slopes.append(beta)
        scaled_a.append(a / scales)
        residual[ix] = yy - xx @ beta - (bb / scales) @ gamma
        leverage[ix] = (1.0 / len(ix) + np.einsum("ij,ji->i", xx, inverse_x)
                        + np.einsum("ij,ji->i", bp[ix], inverse_b[:, ix]))
    if np.any(1.0 - leverage <= 1e-10):
        raise ValueError("HC3 covariance not estimable at unit leverage")
    error2 = (residual / (1.0 - leverage)) ** 2
    gamma_cov = (inverse_b * error2) @ inverse_b.T
    # Full cross-taxon covariance from estimating the shared nuisance effects.
    # Memory is O(N*8 + (taxa*predictors)^2), not O(N*taxa*predictors).
    aa = np.vstack(scaled_a)
    cc = np.vstack([(inverse_x * error2[ix]) @ inverse_b[:, ix].T
                    for ix, _, _, _, inverse_x, _ in blocks])
    covariance = aa @ gamma_cov @ aa.T - cc @ aa.T - aa @ cc.T
    for i, (ix, _, _, _, inverse_x, _) in enumerate(blocks):
        sl = slice(i*p, (i+1)*p)
        covariance[sl, sl] += (inverse_x * error2[ix]) @ inverse_x.T
    covariance = (covariance + covariance.T) / 2.0
    return JointNuisanceSlopes(labels, np.asarray(slopes), covariance, residual,
                              leverage, residual_df, rank)


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
    if values.ndim != 1 or not isinstance(k, int) or k < 1:
        raise ValueError("Moran inputs require a vector and positive integer k")
    if len(latitude) != len(values) or len(longitude) != len(values):
        raise ValueError("Moran coordinate length mismatch")
    if len(values) <= k or not np.all(np.isfinite(values)):
        return float("nan")
    xyz = spherical_basis(latitude, longitude)[:, :3]
    candidates = cKDTree(xyz).query(xyz, k=k+1)[1]
    # With coincident coordinates the query point need not be the first result.
    neighbours = np.array([row[row != i][:k] for i, row in enumerate(candidates)])
    centred = values - float(np.mean(values))
    denominator = float(np.dot(centred, centred))
    if denominator <= 0:
        return float("nan")
    numerator = float(np.sum(centred[:, None] * centred[neighbours]))
    return numerator / (k * denominator)

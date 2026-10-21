"""Spatially cross-fitted original-photo climate novelty control for S2.

A native-range 5th–95th climate envelope evaluated on the SAME native
records used to define that envelope underestimates the true held-out native
novelty rate. This control partitions SOURCE-NATIVE geographic 5-degree cells
into disjoint reference and held-out sets, and compares introduced-range
climate occupancy with the genuinely held-out native occupancy. All
partitioning and thresholds use CLIMATE AND CELLS ONLY, never photo angle,
evolutionary outcomes, original human identities, or posthoc effect signs.

Repeated random source-cell partitions quantify train/test design instability,
NOT independent bootstrap replicates, P values, or causal adaptation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_cached_climate import VARIABLES

SCHEME_VERSION = "v4_S2_geographical_native_cell_crossfit_climate_20261021_v1"
SEED = 20261021
FOLDS = 5
REPEATS = 100
QUANTILES = (.05, .95)
REFERENCE_DEFINITIONS = {
    "annual_precipitation_only": ("chelsa_bio12",),
    "four_without_rsds": tuple(x for x in VARIABLES if x != "chelsa_rsds_mean"),
    "five_original_climate_predictors": tuple(VARIABLES),
    "rsds_only_technical_extrapolation": ("chelsa_rsds_mean",),
}


def geographic_cell_profile(frame: pd.DataFrame) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    if "_cell" not in frame or "obs_id" not in frame:
        raise ValueError("Original site cells and observation identity required")
    if frame.obs_id.isna().any() or frame.obs_id.duplicated().any():
        raise ValueError("Source observation IDs are not one-to-one")
    cells = frame["_cell"].to_numpy()
    if pd.isna(cells).any():
        raise ValueError("Original geographic cells are missing")
    unique, inverse, counts = np.unique(cells, return_inverse=True, return_counts=True)
    if (counts <= 0).any():
        raise ValueError("Invalid observed geographic cell")
    # For each range, the aggregate weight of EVERY distinct source cell is 1.
    weights = 1.0 / counts[inverse]
    return unique, inverse, weights


def _weighted_limits(v: np.ndarray, w: np.ndarray) -> tuple[float,float]:
    if v.ndim != 1 or w.ndim != 1 or len(v) != len(w) or not len(v):
        raise ValueError("Malformed native training climate observations")
    if not np.isfinite(v).all() or not np.isfinite(w).all() or (w <= 0).any():
        raise ValueError("Missing climates or invalid source site weights")
    idx = np.argsort(v, kind="mergesort")
    vs = v[idx]
    ws = w[idx]
    q = (np.cumsum(ws) - ws * 0.5) / ws.sum()
    lo, hi = np.interp(QUANTILES, q, vs, left=vs[0], right=vs[-1])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi < lo:
        raise ValueError("Original source climate envelope is not finite")
    return float(lo), float(hi)


def _outside_matrix(x: np.ndarray, limits: list[tuple[float,float]]) -> np.ndarray:
    if x.ndim != 2 or x.shape[1] != len(limits):
        raise ValueError("Climate dimensions and native reference mismatch")
    mask = np.zeros(len(x), bool)
    for j, (lower, upper) in enumerate(limits):
        mask |= (x[:, j] < lower) | (x[:, j] > upper)
    return mask


def cross_fitted_geographic_novelty(
    native: pd.DataFrame, introduced: pd.DataFrame, *,
    source_name: str, repeats: int = REPEATS, folds: int = FOLDS
) -> dict:
    """Compute climate-only introduced-minus-HELDOUT-native novelty.

    Every native cell is held out once per partition. Introduced pixels are
    evaluated on all folds; native fractions are aggregated across disjoint
    native folds. Nothing is selected using the response or the reported sign.
    """
    if not isinstance(source_name, str) or not source_name.strip():
        raise ValueError("Exact preselected source taxon name required")
    if repeats < 2 or folds < 2 or repeats > 1000 or folds > 10:
        raise ValueError("Registered cross-fitting settings outside bounds")
    required = tuple(VARIABLES)
    for df in (native, introduced):
        if not set(required) <= set(df):
            raise ValueError("The exact five CHELSA climate predictors are missing")
        if df[list(required)].isna().any().any():
            raise ValueError("Missing climate: original pixel cache must not be imputed")
        if not np.isfinite(df[list(required)].to_numpy(float)).all():
            raise ValueError("Non-finite climatic exposure in preselected population")

    n_cells, n_inverse, w_native = geographic_cell_profile(native)
    i_cells, i_inverse, w_intro = geographic_cell_profile(introduced)
    if len(n_cells) < max(25, folds * 5) or len(i_cells) < 15:
        return {
            "status": "STOP_INSUFFICIENT_SOURCE_GEOGRAPHIC_CELLS",
            "source_native_cell_count": len(n_cells),
            "source_introduced_cell_count": len(i_cells),
            "not_equivalent_to_no_niche_shift": True,
        }
    if len(native) < 100 or len(introduced) < 100:
        raise ValueError("Too few original source climate observations to cross-fit")
    n_values = native[list(required)].to_numpy(float)
    i_values = introduced[list(required)].to_numpy(float)
    columns = {name: required.index(name) for name in required}
    seed = int.from_bytes(
        hashlib.sha256(source_name.strip().encode("utf-8")).digest()[:8],
        "little"
    ) ^ SEED
    rng = np.random.default_rng(seed)
    samples = {key:[] for key in REFERENCE_DEFINITIONS}

    for repetition in range(repeats):
        cell_perm = rng.permutation(len(n_cells))
        assignment = np.empty(len(n_cells), dtype=int)
        assignment[cell_perm] = np.arange(len(n_cells)) % folds
        cell_fold = assignment[n_inverse]
        sums = {
            key:{"native_novel_weight":0.0, "native_holdout_weight":0.0,
                 "introduced_novel_rate_sum":0.0}
            for key in REFERENCE_DEFINITIONS
        }
        for k in range(folds):
            train = cell_fold != k
            heldout = ~train
            if int(train.sum()) < 80 or int(heldout.sum()) < 10:
                raise ValueError("Native geographic cross-fit fold has insufficient observations")
            for label, fields in REFERENCE_DEFINITIONS.items():
                ix = [columns[name] for name in fields]
                limits = [
                    _weighted_limits(n_values[train, j], w_native[train])
                    for j in ix
                ]
                n_out = _outside_matrix(n_values[heldout][:, ix], limits)
                i_out = _outside_matrix(i_values[:, ix], limits)
                bucket = sums[label]
                bucket["native_novel_weight"] += float(w_native[heldout] @ n_out.astype(float))
                bucket["native_holdout_weight"] += float(w_native[heldout].sum())
                bucket["introduced_novel_rate_sum"] += float(
                    w_intro @ i_out.astype(float) / w_intro.sum()
                )
        for label, bucket in sums.items():
            native_heldout = (
                bucket["native_novel_weight"] / bucket["native_holdout_weight"]
            )
            introduced_novel = bucket["introduced_novel_rate_sum"] / folds
            samples[label].append(
                (introduced_novel, native_heldout, introduced_novel - native_heldout)
            )

    computed = {}
    for label, rows in samples.items():
        x = np.asarray(rows, float)
        if x.shape != (repeats, 3) or not np.isfinite(x).all():
            raise ValueError("Climate-only cross-fitting failed")
        def quantiles(index):
            return [float(q) for q in np.quantile(x[:, index], (.025, .50, .975))]
        computed[label] = {
            "predictors": list(REFERENCE_DEFINITIONS[label]),
            "introduced_novel_fraction_partition_q025_q50_q975": quantiles(0),
            "genuinely_heldout_native_novel_fraction_partition_q025_q50_q975": quantiles(1),
            "introduced_minus_heldout_native_fraction_partition_q025_q50_q975": quantiles(2),
            "fraction_of_random_partitions_with_positive_excess": float(np.mean(x[:,2] > 0)),
            "head_angle_never_entered_climate_reference_or_fold_selection": True,
            "not_a_sampling_confidence_interval": True,
        }

    return {
        "schema_version": SCHEME_VERSION,
        "status": "EXECUTED_SOURCE_CLIMATE_ONLY_CROSS_FIT_NOT_EVOLUTIONARY_TEST",
        "reference_source_name": source_name,
        "original_native_observations": len(native),
        "original_introduced_observations": len(introduced),
        "original_native_equal_sin_5degree_cells": len(n_cells),
        "original_introduced_equal_sin_5degree_cells": len(i_cells),
        "folds": folds,
        "independent_random_cell_partitions": repeats,
        "native_reference_trained_on_nonheldout_cells": True,
        "introduced_and_native_evaluated_without_angle_values": True,
        "equal_total_weight_per_original_cell": True,
        "climate_reference_variables_from_frozen_model": list(required),
        "results": computed,
        "interpretation_limitations": [
            "Repeated partition quantiles are a design-stability interval, not a 95% ecological sampling confidence interval.",
            "Reference native cells and target introduced cells differ geographically and taxonomically, and introduction origin/source populations are unknown.",
            "Outside any single marginal 5–95% range is not necessarily novel multivariate climatic density.",
            "Training native photo observations are selected by an image-orientation-eligible phenotype; source sampling/recognition bias remains.",
            "All inputs are annual or seasonal climate normals, not actual wetting during anthesis.",
            "Photographed head orientation is camera-relative and no genetic change, plasticity or adaptation is established.",
        ],
    }

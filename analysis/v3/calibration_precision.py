"""Monte Carlo precision over independent outer datasets, valid at repeated looks.

For fixed f in (0,1], 1-f+f*Z/m is nonnegative and has conditional
expectation one when E[Z|past]=m. Its product and a fixed mixture of such
products are test martingales. Apply Ville to both directions with alpha/2;
union across the 16 scenario/grid/metric sequences. No independence among
coefficients, metrics or grids is needed. This is a simple fixed mixture,
not a reimplementation of the optimized betting algorithm in qkad009.
"""
from __future__ import annotations

import json
import math

import numpy as np

from .workflow import ROOT, canonical_digest

BET_FRACTIONS = np.arange(1, 21, dtype=float) / 20


def amendment(root=ROOT):
    spec = json.loads((root / "analysis/v3/calibration_precision_amendment.json").read_text(encoding="utf-8"))
    historical = json.loads((root / spec["historical_execution_contract"]).read_text(encoding="utf-8"))
    if canonical_digest(historical) != spec["historical_contract_canonical_sha256"]:
        raise ValueError("Historical execution contract differs; do not silently amend the experiment")
    if (spec["confidence_sequences"] != len(historical["scenarios"]) * len(historical["grid_degrees"]) * 2
            or spec["alpha_per_sequence"] != spec["simultaneous_error_budget"] / spec["confidence_sequences"]
            or not np.array_equal(spec["fixed_bet_fractions"], BET_FRACTIONS)):
        raise ValueError("Fixed mixture or simultaneous error allocation differs")
    if spec["thresholds"] != {"fwer_half_width_max": .025, "coverage_half_width_max": .015,
                              "fwer_upper_max": .075, "coverage_lower_min": .90}:
        raise ValueError("Do not retune scientific admission thresholds")
    if spec["outer_range"] != {"batch_size_per_scenario": 25, "minimum_per_scenario": 100, "maximum_per_scenario": 400}:
        raise ValueError("Do not retune the outer range")
    return spec


def _bounded(values):
    data = np.asarray(values, dtype=float)
    if data.ndim != 1 or not len(data) or not np.isfinite(data).all() or np.any((data < 0) | (data > 1)):
        raise ValueError("One finite bounded value per independent outer dataset is required")
    return data


def _log_lower_capital(data, mean):
    """Log mixture capital against a too-small mean, for 0 < mean <= 1."""
    factors = (1 - BET_FRACTIONS[:, None]) + BET_FRACTIONS[:, None] * (data / mean)
    with np.errstate(divide="ignore"):
        log_products = np.log(factors).sum(axis=1)
    largest = float(log_products.max())
    return largest + math.log(float(np.exp(log_products - largest).mean()))


def bounded_mean_sequence(values, *, alpha):
    data = _bounded(values)
    if isinstance(alpha, bool) or not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("A fixed error budget in (0,1) is required")
    boundary = math.log(2 / alpha)

    def lower_bound(sample):
        if not np.any(sample > 0):
            return 0.0
        lo, hi = 0.0, float(sample.mean())
        # Capital is decreasing in m and <= 1 at the sample mean by Jensen.
        # Return the outside (lower) bracket; never shrink the confidence set.
        for _ in range(55):
            mid = (lo + hi) / 2
            if _log_lower_capital(sample, mid) >= boundary:
                lo = mid
            else:
                hi = mid
        return max(0.0, lo - 1e-12)

    low = lower_bound(data)
    high = min(1.0, 1 - lower_bound(1 - data))
    return {"mean": float(data.mean()), "low": low, "high": high,
            "half_width": (high - low) / 2, "independent_outer_datasets": len(data),
            "alpha": alpha, "time_uniform": True}


def precision_decision(false_family, coverage_fraction, *, root=ROOT):
    spec = amendment(root)
    false_family, coverage_fraction = _bounded(false_family), _bounded(coverage_fraction)
    if len(false_family) != len(coverage_fraction) or not np.isin(false_family, [0., 1.]).all():
        raise ValueError("Aligned outer-level family indicators and coverage fractions required")
    count = len(false_family)
    if count > spec["outer_range"]["maximum_per_scenario"]:
        raise ValueError("Frozen maximum outer range exceeded")
    fwer = bounded_mean_sequence(false_family, alpha=spec["alpha_per_sequence"])
    coverage = bounded_mean_sequence(coverage_fraction, alpha=spec["alpha_per_sequence"])
    precision = count >= 100 and fwer["half_width"] <= .025 and coverage["half_width"] <= .015
    admission = fwer["high"] <= .075 and coverage["low"] >= .90
    return {"fwer_confidence_sequence": fwer, "coverage_confidence_sequence": coverage,
            "precision_satisfied": precision, "admission_satisfied": admission,
            "independent_unit": "generated_outer_dataset", "coefficient_independence_assumed": False,
            "simultaneous_sequence_count": 16, "simultaneous_error_budget": .05}

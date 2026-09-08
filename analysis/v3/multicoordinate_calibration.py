"""Synthetic final-geometry inputs and precision checks for v3 calibration.

No empirical response or trait-environment association is read. The generator
exercises the frozen 9-predictor order, exact 13-column nuisance constructor and
1/8/4 primary module dimensions before any empirical ecological fit.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .environment_model import definition as environment_definition
from .module_ecology import SCALES
from .nuisance_design import matrix as nuisance_matrix
from .workflow import ROOT

CONTRACT = ROOT / "analysis/v3/final_module_calibration_contract.json"
MODULES = ("orientation", "visible_colour", "gross_shape")


def definition(root: Path = ROOT):
    contract = json.loads((root / "analysis/v3/final_module_calibration_contract.json").read_text(encoding="utf-8"))
    env = environment_definition(root)
    if contract["status"] != "prepared_before_empirical_trait_environment_fitting_not_scheduled_until_tail_mechanics_closes":
        raise ValueError("Final calibration contract is not in the prepared pre-outcome state")
    if contract["predictors"] != len(env["variables"]) or contract["family_slots"] != 36:
        raise ValueError("Calibration predictor/family geometry differs from frozen environment definition")
    if tuple(contract["module_dimensions"]) != MODULES:
        raise ValueError("Calibration module order differs")
    if contract["bootstrap_replicates"] != 999:
        raise ValueError("Tail-resolved calibration must retain 999 planned draws")
    if contract["ecological_models_executed"] != 0 or contract["empirical_trait_environment_values_read"] != 0:
        raise ValueError("Calibration contract no longer represents the pre-outcome boundary")
    return contract, env


def _layout(rng: np.random.Generator, taxa: int):
    counts = np.r_[1, 3, np.clip(np.rint(rng.lognormal(3.7, 0.4, taxa - 2)), 18, 100).astype(int)]
    groups = np.repeat(np.arange(taxa), counts)
    n = len(groups)
    region = rng.integers(0, 16, n)
    centers_lat = np.repeat([27.5, 37.5, 47.5, 57.5], 4)
    centers_lon = np.tile([-112.5, -52.5, 27.5, 117.5], 4)
    side = rng.choice([-1.0, 1.0], n)
    lat = centers_lat[region] + side + rng.uniform(-0.1, 0.1, n)
    lon = centers_lon[region] + side + rng.uniform(-0.1, 0.1, n)
    return counts, groups, region, lat, lon


def _environment(rng: np.random.Generator, groups: np.ndarray, p: int):
    n = len(groups)
    global_factor = rng.normal(size=n)
    moisture = rng.normal(size=n)
    heat = rng.normal(size=n)
    x = np.empty((n, p), dtype=float)
    for j in range(4):
        x[:, j] = 0.25 * global_factor + 0.60 * moisture + math.sqrt(1 - 0.25**2 - 0.60**2) * rng.normal(size=n)
    x[:, 4] = 0.30 * global_factor + 0.30 * heat + math.sqrt(1 - 0.30**2 - 0.30**2) * rng.normal(size=n)
    for j in range(5, 8):
        x[:, j] = 0.25 * global_factor + 0.60 * heat + math.sqrt(1 - 0.25**2 - 0.60**2) * rng.normal(size=n)
    x[:, 8] = 0.15 * global_factor + math.sqrt(1 - 0.15**2) * rng.normal(size=n)
    taxon_means = 0.65 * rng.normal(size=(int(groups.max()) + 1, p))
    x += taxon_means[groups]
    x = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)
    return x


def _raw_nuisance(rng: np.random.Generator, x: np.ndarray, lat: np.ndarray, lon: np.ndarray):
    n = len(x)
    phase = rng.uniform(0, 2 * np.pi, n)
    sin_doy, cos_doy = np.sin(phase), np.cos(phase)
    year = rng.integers(1980, 2027, n).astype(float)
    old = rng.random(n) < 0.02
    year[old] = rng.integers(1970, 1980, int(old.sum()))
    log_size = 4.8 + 0.18 * x[:, 0] + 0.10 * x[:, 4] + rng.normal(0, 0.35, n)
    size = np.exp(log_size)
    log_sharp = 4.0 + 0.15 * x[:, 4] + 0.12 * x[:, 5] + rng.normal(0, 0.4, n)
    sharpness = np.maximum(0.0, np.exp(log_sharp) - 1.0)
    return sin_doy, cos_doy, year, lat, lon, size, sharpness


def _duplicate_components(groups, arrays):
    components = np.arange(len(groups)).astype(str)
    for taxon in range(2, int(groups.max()) + 1, 8):
        positions = np.flatnonzero(groups == taxon)
        if len(positions) < 2:
            continue
        i, j = positions[:2]
        components[j] = components[i]
        for array in arrays:
            array[j] = array[i]
    return components


def _fixed_coefficients(module: str, d: int, p: int, scenario: str):
    coordinate = np.arange(d, dtype=float)
    sign = np.where((coordinate.astype(int) % 2) == 0, 1.0, -1.0)
    scale = 1.0 / np.sqrt(max(1, d))
    beta_w = np.zeros((d, p), dtype=float)
    beta_w[:, 5] = scale * (0.28 + 0.02 * coordinate) * sign
    beta_w[:, 6] = scale * (-0.20 + 0.01 * coordinate)
    beta_w[:, 7] = scale * (0.16 + 0.015 * coordinate) * sign
    beta_w[:, 8] = scale * (0.18 + 0.01 * coordinate)
    beta_a = beta_w.copy()
    if scenario == "scale_difference":
        beta_a[:, 4] = scale * (0.42 + 0.02 * coordinate) * sign
    return beta_w, beta_a


def generate(scenario: str, seed: int, *, root: Path = ROOT):
    contract, env = definition(root)
    if scenario not in contract["scenarios"]:
        raise ValueError("Scenario outside final calibration contract")
    rng = np.random.default_rng(seed)
    taxa, p = int(contract["taxa"]), int(contract["predictors"])
    _, groups, region, lat, lon = _layout(rng, taxa)
    x = _environment(rng, groups, p)
    raw = list(_raw_nuisance(rng, x, lat, lon))
    components = _duplicate_components(groups, [x, *raw])
    sin_doy, cos_doy, year, lat, lon, size, sharpness = raw
    nuisance, nuisance_names = nuisance_matrix(
        sin_doy=sin_doy, cos_doy=cos_doy, observed_year=year,
        latitude=lat, longitude=lon, size=size, sharpness=sharpness, root=root,
    )
    if nuisance.shape[1] != contract["nuisance_columns"]:
        raise ValueError("Generated nuisance geometry differs from final contract")

    n = len(groups)
    global_error = rng.normal(size=n)
    region_error = rng.normal(size=(16, 1)) if scenario != "iid_null" else np.zeros((16, 1))
    responses, truths = {}, {}
    shared_taxon_slopes = rng.normal(size=(taxa, 1, p))
    for module_index, module in enumerate(MODULES):
        d = int(contract["module_dimensions"][module])
        beta_w, beta_a = _fixed_coefficients(module, d, p, scenario)
        tau = 0.75 if scenario == "heterogeneous_slopes" else 0.40
        coordinate_scale = np.linspace(0.8, 1.2, d)[None, :, None]
        deviation = tau * coordinate_scale * (
            0.55 * shared_taxon_slopes + math.sqrt(1 - 0.55**2) * rng.normal(size=(taxa, d, p))
        )
        intercept = 0.45 * rng.normal(size=(taxa, 1)) + 0.55 * rng.normal(size=(taxa, d))
        gamma = np.linspace(0.04, 0.18, nuisance.shape[1])[:, None] * np.linspace(0.8, 1.2, d)[None, :]
        coordinate_factor = np.linspace(0.65, 1.15, d)
        error = 0.35 * global_error[:, None] * coordinate_factor[None, :] + 0.65 * rng.normal(size=(n, d))
        error += 0.40 * region_error[region] * coordinate_factor[None, :]
        y = np.empty((n, d), dtype=float)
        for taxon in range(taxa):
            index = np.flatnonzero(groups == taxon)
            mean_x = x[index].mean(axis=0)
            y[index] = ((x[index] - mean_x) @ (beta_w + deviation[taxon]).T
                        + mean_x @ beta_a.T + nuisance[index] @ gamma
                        + intercept[taxon] + error[index])
        # Exact content duplicates must have matching response values too.
        for taxon in range(2, taxa, 8):
            positions = np.flatnonzero(groups == taxon)
            if len(positions) >= 2 and components[positions[0]] == components[positions[1]]:
                y[positions[1]] = y[positions[0]]
        responses[module] = y
        truths[module] = np.stack([beta_w, beta_a, beta_a - beta_w])

    return {
        "responses": responses,
        "truths": truths,
        "predictors": x,
        "predictor_names": tuple(env["variables"]),
        "taxa": groups,
        "nuisance": nuisance,
        "nuisance_names": tuple(nuisance_names),
        "components": components,
        "latitude": lat,
        "longitude": lon,
        "region": region,
    }


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054):
    if not isinstance(successes, int) or not isinstance(total, int) or total <= 0 or not 0 <= successes <= total:
        raise ValueError("Wilson counts are invalid")
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half), half


def precision_decision(*, outer_with_false_family: int, outer_total: int,
                       covered_coefficients: int, coefficient_total: int, root: Path = ROOT):
    contract, _ = definition(root)
    rule = contract["sequential_outer_rule"]
    fwer = wilson_interval(outer_with_false_family, outer_total)
    coverage = wilson_interval(covered_coefficients, coefficient_total)
    min_n, max_n = rule["minimum_outer_replicates_per_scenario"], rule["maximum_outer_replicates_per_scenario"]
    precision = outer_total >= min_n and fwer[2] <= 0.025 and coverage[2] <= 0.015
    admission = fwer[1] <= 0.075 and coverage[0] >= 0.90
    stop = outer_total >= max_n or (precision and admission)
    return {
        "fwer_wilson95": {"low": fwer[0], "high": fwer[1], "half_width": fwer[2]},
        "coverage_wilson95": {"low": coverage[0], "high": coverage[1], "half_width": coverage[2]},
        "precision_satisfied": precision,
        "admission_satisfied": admission,
        "stop_for_precision_and_admission": stop and admission,
        "stop_at_max_without_admission": outer_total >= max_n and not admission,
        "next_batch_required": not stop,
    }

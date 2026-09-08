"""Frozen calendar, imaging and spherical nuisance construction for Chapter 1 v3.

No file entrypoint and no ecological authorization. Callers must first declare one
realized module cohort with endpoint/module-matched imaging quality and finite source
annotations. This function never drops or imputes rows.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .hierarchical_ecology import SPATIAL_TERMS, spherical_basis
from .workflow import ROOT, canonical_digest

CONTRACT = ROOT / "analysis/v3/nuisance_design_contract.json"


def definition(root: Path = ROOT):
    contract = json.loads((root / "analysis/v3/nuisance_design_contract.json").read_text(encoding="utf-8"))
    expected = [
        "calendar_sin_doy", "calendar_cos_doy", "observation_year_decade",
        "log_head_min_dimension_px", "log1p_head_laplacian_variance",
        *[f"spatial_{name}" for name in SPATIAL_TERMS],
    ]
    if contract["status"] != "fixed_before_empirical_trait_environment_fitting":
        raise ValueError("Nuisance design is not frozen")
    if contract["matrix_order"] != expected:
        raise ValueError("Nuisance matrix order differs from implementation")
    if contract["spatial"]["basis"] != list(SPATIAL_TERMS):
        raise ValueError("Spherical basis differs from hierarchical implementation")
    if contract["ecological_models_executed"] != 0 or contract["empirical_trait_environment_values_read"] != 0:
        raise ValueError("Nuisance contract no longer represents the pre-outcome boundary")
    return {
        "status": "NUISANCE_DESIGN_VERIFIED_NOT_FITTED",
        "contract_canonical_sha256": canonical_digest(contract),
        "columns": tuple(expected),
        "ecological_fitting_authorized": False,
    }


def matrix(*, sin_doy, cos_doy, observed_year, latitude, longitude, size, sharpness, root: Path = ROOT):
    """Construct all predeclared nuisance columns without row mutation.

    `size` and `sharpness` must already be the module-matched observation summaries
    from the same response support. Missing or invalid values are errors so callers
    must report non-eligibility explicitly instead of silently dropping rows here.
    """
    spec = definition(root)
    arrays = [np.asarray(v, dtype=float) for v in
              (sin_doy, cos_doy, observed_year, latitude, longitude, size, sharpness)]
    if any(a.ndim != 1 for a in arrays) or len({len(a) for a in arrays}) != 1 or not len(arrays[0]):
        raise ValueError("Nuisance inputs must be non-empty aligned vectors")
    if not all(np.isfinite(a).all() for a in arrays):
        raise ValueError("Nuisance support must be finite; no silent row deletion or imputation")
    sin_value, cos_value, year, lat, lon, size_value, sharpness_value = arrays
    if np.any(size_value <= 0):
        raise ValueError("Module-matched head size must be positive before log transform")
    if np.any(sharpness_value < 0):
        raise ValueError("Module-matched Laplacian variance cannot be negative")
    if np.any((lat < -90) | (lat > 90) | (lon < -180) | (lon > 180)):
        raise ValueError("Public analysis coordinates are outside valid geographic bounds")
    if np.any(np.abs(sin_value) > 1 + 1e-12) or np.any(np.abs(cos_value) > 1 + 1e-12):
        raise ValueError("Calendar harmonics are outside the annotation definition")

    spatial = spherical_basis(lat, lon)
    out = np.column_stack([
        sin_value,
        cos_value,
        (year - 2010.0) / 10.0,
        np.log(size_value),
        np.log1p(sharpness_value),
        spatial,
    ])
    if out.shape[1] != len(spec["columns"]) or not np.isfinite(out).all():
        raise ValueError("Frozen nuisance matrix construction failed")
    return out, spec["columns"]

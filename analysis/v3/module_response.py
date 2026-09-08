"""Primary v3 response-module transforms fixed before empirical ecology.

The helper is measurement-only. It never reads environment columns, drops rows,
changes endpoint support or authorizes fitting.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .workflow import ROOT, canonical_digest

CONTRACT = ROOT / "analysis/v3/module_response_contract.json"
DECISION = ROOT / "analysis/v3/measurement_qualification_decision_20260908.json"


def definition(root: Path = ROOT):
    contract = json.loads((root / "analysis/v3/module_response_contract.json").read_text(encoding="utf-8"))
    decision = json.loads((root / "analysis/v3/measurement_qualification_decision_20260908.json").read_text(encoding="utf-8"))
    if contract["status"] != "fixed_before_empirical_trait_environment_fitting":
        raise ValueError("Module response geometry is not frozen")
    admitted = {
        row["endpoint_id"] for row in decision["endpoints"]
        if row.get("ecological_route") == "stream_original_required"
    }
    primary = contract["primary_modules"]
    expected = set().union(*(set(primary[module]["endpoints"]) for module in primary))
    if len(expected) != 13 or not expected <= admitted:
        raise ValueError("Primary response geometry is not supported by the frozen measurement decision")
    if set(primary) != {"orientation", "visible_colour", "gross_shape"}:
        raise ValueError("Primary module inventory differs")
    if contract["ecological_models_executed"] != 0 or contract["empirical_trait_environment_values_read"] != 0:
        raise ValueError("Module geometry no longer represents the pre-outcome boundary")
    return {
        "status": "PRIMARY_MODULE_RESPONSE_GEOMETRY_VERIFIED_NOT_FITTED",
        "contract_canonical_sha256": canonical_digest(contract),
        "measurement_decision_canonical_sha256": canonical_digest(decision),
        "modules": primary,
        "ecological_fitting_authorized": False,
    }


def matrix(frame, module: str, *, root: Path = ROOT):
    """Return the fixed response coordinate matrix for one exact joint-support cohort."""
    spec = definition(root)
    if module not in spec["modules"]:
        raise ValueError("Unknown primary module")
    rule = spec["modules"][module]
    missing = [name for name in rule["endpoints"] if name not in frame]
    if missing:
        raise ValueError(f"Missing module endpoints: {missing}")
    values = frame[rule["endpoints"]].to_numpy(dtype=float)
    if values.ndim != 2 or len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Declare finite non-empty exact module support; no silent row deletion")

    if module == "visible_colour":
        # Endpoint order is frozen in the contract: L*, chroma, hue sin/cos, four fractions.
        hue = values[:, 2:4]
        resultant = np.sqrt(np.sum(hue * hue, axis=1))
        if np.any(resultant > 1.0 + 1e-6):
            raise ValueError("Circular hue components exceed the registered unit-disc bound")
        composition = values[:, 4:8]
        if np.any(composition < 0) or np.any(np.abs(composition.sum(axis=1) - 1.0) > 1e-6):
            raise ValueError("Colour composition is not closed; never renormalize a failed row")
        values = np.column_stack([values[:, :4], np.sqrt(composition)])

    if values.shape[1] != len(rule["coordinates"]):
        raise ValueError("Response coordinate geometry differs from contract")
    return values, tuple(rule["coordinates"])

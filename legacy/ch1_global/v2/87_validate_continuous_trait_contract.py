#!/usr/bin/env python3
"""Validate the category-free continuous trait contract for Chapter 1.

The contract separates a biological construct from an image-derived numerical
measurement.  It prevents legacy state labels, duplicated formulas disguised as
different traits, and uncalibrated proxies from silently entering the primary
analysis tier.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "endpoint_id",
    "module",
    "biological_construct",
    "measurement_variable",
    "observation_variable",
    "species_variable",
    "source_status_column",
    "unit",
    "lower_bound",
    "upper_bound",
    "circular_group",
    "compositional_group",
    "formula_id",
    "image_requirement",
    "min_head_dimension_px",
    "analysis_tier",
    "validation_status",
    "category_free",
    "interpretation",
    "not_interpretable_as",
}
ALLOWED_TIERS = {"primary", "candidate", "validation_only", "descriptive_only"}
ALLOWED_VALIDATION = {
    "frozen_technical",
    "requires_botanical_calibration",
    "reference_only",
}
TRUE_VALUES = {"true", "1", "yes"}
FORBIDDEN_CATEGORY_TOKENS = {
    "upright",
    "ascending",
    "horizontal",
    "nodding",
    "pendant",
    "globose",
    "hemispherical",
    "ovoid",
    "cylindrical",
    "urceolate",
    "unarmed",
    "short_spined",
    "long_spined",
    "appressed",
    "spreading",
    "recurved",
    "glabrous",
    "webby",
    "tomentose",
    "glandular",
}


def _text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def _bool(value: Any) -> bool:
    return _text(value).lower() in TRUE_VALUES


def validate_contract(frame: pd.DataFrame) -> dict[str, Any]:
    missing = sorted(REQUIRED_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"Contract is missing required columns: {missing}")
    if frame.empty:
        raise ValueError("Continuous trait contract is empty")

    work = frame.copy()
    for column in REQUIRED_COLUMNS:
        work[column] = work[column].map(_text)

    for column in (
        "endpoint_id",
        "measurement_variable",
        "observation_variable",
        "species_variable",
        "formula_id",
    ):
        if work[column].eq("").any():
            raise ValueError(f"{column} must not be blank")
        duplicated = work.loc[work[column].duplicated(keep=False), column].tolist()
        if duplicated:
            raise ValueError(f"{column} must be unique; duplicates={sorted(set(duplicated))}")

    invalid_tiers = sorted(set(work["analysis_tier"]) - ALLOWED_TIERS)
    if invalid_tiers:
        raise ValueError(f"Unsupported analysis_tier values: {invalid_tiers}")
    invalid_validation = sorted(set(work["validation_status"]) - ALLOWED_VALIDATION)
    if invalid_validation:
        raise ValueError(f"Unsupported validation_status values: {invalid_validation}")
    if not work["category_free"].map(_bool).all():
        bad = work.loc[~work["category_free"].map(_bool), "endpoint_id"].tolist()
        raise ValueError(f"Every endpoint must be category_free=true: {bad}")

    primary = work["analysis_tier"].eq("primary")
    if not work.loc[primary, "validation_status"].eq("frozen_technical").all():
        bad = work.loc[primary & ~work["validation_status"].eq("frozen_technical"), "endpoint_id"].tolist()
        raise ValueError(f"Primary endpoints must have frozen_technical validation: {bad}")

    minimum = pd.to_numeric(work["min_head_dimension_px"], errors="coerce")
    if minimum.isna().any() or minimum.lt(1).any():
        raise ValueError("min_head_dimension_px must be a positive number")
    lower = pd.to_numeric(work["lower_bound"], errors="coerce")
    upper = pd.to_numeric(work["upper_bound"], errors="coerce")
    finite_pair = lower.notna() & upper.notna()
    if (lower[finite_pair] > upper[finite_pair]).any():
        raise ValueError("lower_bound exceeds upper_bound for at least one endpoint")

    for group, part in work.loc[work["circular_group"].ne("")].groupby("circular_group"):
        if len(part) != 2:
            raise ValueError(f"Circular group {group!r} must contain exactly sine and cosine rows")
        formulas = set(part["formula_id"])
        if not any("sine" in value for value in formulas) or not any("cosine" in value for value in formulas):
            raise ValueError(f"Circular group {group!r} lacks explicit sine/cosine formula ids")

    for group, part in work.loc[work["compositional_group"].ne("")].groupby("compositional_group"):
        if len(part) < 2:
            raise ValueError(f"Compositional group {group!r} must contain at least two components")
        if part["analysis_tier"].isin({"primary", "candidate"}).any():
            raise ValueError(
                f"Composition components in {group!r} cannot enter univariate primary/candidate models"
            )

    narrative = (
        work["biological_construct"]
        + " "
        + work["interpretation"]
        + " "
        + work["measurement_variable"]
    ).str.lower()
    for token in FORBIDDEN_CATEGORY_TOKENS:
        # Category words are allowed only in the explicit non-interpretation warning.
        if narrative.str.contains(fr"\b{token}\b", regex=True).any():
            bad = work.loc[narrative.str.contains(fr"\b{token}\b", regex=True), "endpoint_id"].tolist()
            raise ValueError(f"Legacy category token {token!r} leaked into numerical definitions: {bad}")

    if work["interpretation"].eq("").any() or work["not_interpretable_as"].eq("").any():
        raise ValueError("Every endpoint requires positive and negative semantic definitions")

    modules = work.groupby("module")["endpoint_id"].size().sort_index().astype(int).to_dict()
    tiers = work.groupby("analysis_tier")["endpoint_id"].size().sort_index().astype(int).to_dict()
    return {
        "status": "valid",
        "n_endpoints": int(len(work)),
        "n_modules": int(work["module"].nunique()),
        "module_counts": modules,
        "analysis_tier_counts": tiers,
        "n_biological_dimensions": int(
            len(work.loc[work["circular_group"].eq("")])
            + work.loc[work["circular_group"].ne(""), "circular_group"].nunique()
        ),
        "category_columns_allowed_in_analysis": 0,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    frame = pd.read_csv(args.contract, dtype=str, keep_default_na=False)
    report = validate_contract(frame)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

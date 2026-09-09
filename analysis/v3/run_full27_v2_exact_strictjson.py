#!/usr/bin/env python3
"""Strict-JSON launcher for the exact-v2 full-27 lane.

The scientific runner can legitimately return NaN in scalar-only columns of
circular-result rows.  Those missing fields are presentation nulls, not failed
models.  This launcher replaces only non-finite values in the report-row view
with None before the core runner serializes its final JSON.  Numerical CSVs,
permutations, BH correction and all ecological calculations are untouched.
"""
from __future__ import annotations

import math

from analysis.v3 import run_full27_v2_exact_recovered as core

_original_signal_rows = core.signal_rows


def strict_signal_rows(path, scope):
    rows = _original_signal_rows(path, scope)
    clean = []
    for row in rows:
        out = {}
        for key, value in row.items():
            if isinstance(value, float) and not math.isfinite(value):
                out[key] = None
            else:
                out[key] = value
        clean.append(out)
    return clean


def main() -> int:
    core.signal_rows = strict_signal_rows
    return int(core.main())


if __name__ == "__main__":
    raise SystemExit(main())

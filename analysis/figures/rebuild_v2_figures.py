#!/usr/bin/env python3
"""Run the frozen v2 figure builder against tracked reproducibility artifacts.

The historical implementation is kept in ``build_v2_figures.py``. This wrapper
sets repository-root paths explicitly so the implementation can live under
``analysis/figures`` while reading frozen outputs from ``analysis_outputs`` and
writing generated scientific figures to ``reproducibility/figures``.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
IMPLEMENTATION = HERE / "build_v2_figures.py"

spec = importlib.util.spec_from_file_location("azami_v2_figure_builder", IMPLEMENTATION)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load figure builder: {IMPLEMENTATION}")
impl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(impl)

impl.ROOT = ROOT
impl.DEFAULT_INPUT = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"
impl.DEFAULT_OUTPUT = ROOT / "reproducibility" / "figures"
impl.SOURCE = impl.DEFAULT_OUTPUT / "source"


if __name__ == "__main__":
    impl.main()

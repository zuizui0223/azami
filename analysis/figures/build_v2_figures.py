#!/usr/bin/env python3
"""Canonical entry point for rebuilding the frozen v2 result figures.

Generated figures are written to ``reproducibility/figures``. Prepublication
manuscript directories are never created by this entry point. The GitHub Actions
verification workflow rebuilds the same figure set in scratch space and does
not upload another expiring artifact.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
IMPLEMENTATION = HERE / "_build_v2_figures_impl.py"


def _load_impl():
    spec = importlib.util.spec_from_file_location("azami_v2_figure_builder_impl", IMPLEMENTATION)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load figure implementation: {IMPLEMENTATION}")
    impl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(impl)
    impl.ROOT = ROOT
    impl.DEFAULT_INPUT = ROOT / "analysis_outputs" / "v2_full27_environment_atlas_2026-08-27"
    impl.DEFAULT_OUTPUT = ROOT / "reproducibility" / "figures"
    impl.SOURCE = impl.DEFAULT_OUTPUT / "source"
    return impl


def main() -> None:
    _load_impl().main()


if __name__ == "__main__":
    main()

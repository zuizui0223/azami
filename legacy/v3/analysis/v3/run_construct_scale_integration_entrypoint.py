#!/usr/bin/env python3
"""Compatibility launcher for run_construct_scale_integration.

The core integration module was written against an earlier four-value draft of
run_biological_axis_reanalysis.load_data(); the committed axis module returns
(traits, environment).  Keep the scientific implementation unchanged and adapt
only that call contract here.  This wrapper can be removed once the core file is
cleanly edited to unpack two values directly.
"""
from __future__ import annotations

from analysis.v3 import run_construct_scale_integration as core

_original_load_data = core.axis.load_data


def _load_data_compat(*args, **kwargs):
    traits, environment = _original_load_data(*args, **kwargs)
    return None, traits, environment, None


core.axis.load_data = _load_data_compat

if __name__ == "__main__":
    raise SystemExit(core.main())

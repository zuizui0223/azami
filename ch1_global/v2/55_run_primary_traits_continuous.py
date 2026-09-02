#!/usr/bin/env python3
"""OpenCV compatibility layer and production entrypoint.

OpenCV releases differ in whether ``HoughLinesP`` returns an ``(N, 1, 4)`` or
``(N, 4)`` array, especially when only one line is detected. The compatibility
layer normalizes that API boundary. When executed as a script, this file now
routes production runs through the corrected v2 colour and orientation logic.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

import cv2
import numpy as np


_ORIGINAL_HOUGH_LINES_P = cv2.HoughLinesP


def _normalized_hough_lines_p(*args: Any, **kwargs: Any):
    lines = _ORIGINAL_HOUGH_LINES_P(*args, **kwargs)
    if lines is None:
        return None
    return np.asarray(lines).reshape(-1, 1, 4)


cv2.HoughLinesP = _normalized_hough_lines_p

_ENGINE_PATH = Path(__file__).with_name("52_measure_primary_traits_continuous.py")
_SPEC = importlib.util.spec_from_file_location(
    "ch1_continuous_primary_measurement_engine", _ENGINE_PATH
)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"Unable to load continuous measurement engine: {_ENGINE_PATH}")
MEASURE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(MEASURE)


if __name__ == "__main__":
    v2_path = Path(__file__).with_name("56_run_primary_traits_continuous_v2.py")
    v2_spec = importlib.util.spec_from_file_location(
        "ch1_continuous_primary_measurement_v2_entrypoint", v2_path
    )
    if v2_spec is None or v2_spec.loader is None:
        raise RuntimeError(f"Unable to load corrected v2 entrypoint: {v2_path}")
    v2_module = importlib.util.module_from_spec(v2_spec)
    v2_spec.loader.exec_module(v2_module)
    v2_module.main()

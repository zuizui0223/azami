#!/usr/bin/env python3
"""Run continuous primary-trait measurement with OpenCV shape compatibility.

OpenCV releases differ in whether ``HoughLinesP`` returns an ``(N, 1, 4)`` or
``(N, 4)`` array, especially when only one line is detected. The measurement
engine historically expects ``(N, 1, 4)``. This entrypoint normalizes that API
boundary without changing any biological measurement or threshold.
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
    MEASURE.main()

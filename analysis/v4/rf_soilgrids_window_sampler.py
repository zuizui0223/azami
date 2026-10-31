"""Optional equivalent native-pixel SoilGrids sampling with bounded raster windows.

This is not a new dataset or a fitted RF. Source raster must be separately
SHA/source authenticated. Returns exact source pixels, not interpolated values.
Original S1 photo coordinates and any per-pixel outputs must remain private.
"""
from __future__ import annotations
import time
import numpy as np
from rasterio.transform import rowcol
from rasterio.windows import Window


def sample_projected_pixels_windowed(src, xs, ys, *, window_pixels=256,
                                    max_seconds=750.0, clock=time.monotonic):
    """Equivalent to rasterio.sample(masked=True) at in-bounds pixels.

    Read each raster window once and re-use native pixels for all source points
    falling into that window; no resampling, neighbour search or interpolation.
    Time budget is checked before each potentially slow remote read; GDAL's
    own per-HTTP-request timeout remains a separate requirement.
    """
    if src.crs is None or src.count != 1:
        raise ValueError("SOURCE_SINGLE_BAND_CRS_REQUIRED")
    if isinstance(window_pixels, bool) or not isinstance(window_pixels, int) or not 1 <= window_pixels <= 512:
        raise ValueError("INVALID_FIXED_WINDOW_DIMENSION")
    if max_seconds <= 0 or not np.isfinite(max_seconds):
        raise ValueError("INVALID_REMOTE_TIME_BUDGET")
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if xs.ndim != 1 or ys.shape != xs.shape or not np.isfinite(xs).all() or not np.isfinite(ys).all():
        raise ValueError("NONFINITE_OR_MISMATCHED_PROJECTED_POINT_ARRAYS")
    start = clock()
    rr, cc = rowcol(src.transform, xs.tolist(), ys.tolist())
    rr, cc = np.asarray(rr, dtype=np.int64), np.asarray(cc, dtype=np.int64)
    result = np.full(len(xs), np.nan, dtype=float)
    valid = (rr >= 0) & (rr < src.height) & (cc >= 0) & (cc < src.width)
    indices = np.flatnonzero(valid)
    ncols = (src.width + window_pixels - 1) // window_pixels
    block_id = rr[indices] // window_pixels * ncols + cc[indices] // window_pixels
    order = np.argsort(block_id, kind="stable")
    idx_sorted = indices[order]
    ids_sorted = block_id[order]
    cuts = np.r_[0, np.flatnonzero(np.diff(ids_sorted)) + 1, len(idx_sorted)]
    windows_read = 0
    for l, r in zip(cuts[:-1], cuts[1:]):
        if clock() - start > max_seconds:
            raise TimeoutError("SOILGRIDS_NATIVE_PIXEL_WINDOW_BUDGET_EXCEEDED")
        chosen = idx_sorted[l:r]
        row0 = int(rr[chosen[0]] // window_pixels) * window_pixels
        col0 = int(cc[chosen[0]] // window_pixels) * window_pixels
        h = min(window_pixels, src.height - row0)
        w = min(window_pixels, src.width - col0)
        block = src.read(1, window=Window(col0, row0, w, h), masked=True)
        raw = block[rr[chosen] - row0, cc[chosen] - col0]
        mask = np.ma.getmaskarray(raw)
        nums = np.asarray(np.ma.getdata(raw), dtype=float)
        allowed = ~mask & np.isfinite(nums)
        result[chosen[allowed]] = nums[allowed]
        windows_read += 1
    return result, {"n_source_points": int(len(xs)),
                    "n_inside_source_grid": int(valid.sum()),
                    "n_unique_source_windows": int(windows_read),
                    "n_native_pixel_values_finite": int(np.isfinite(result).sum()),
                    "window_pixels": int(window_pixels),
                    "no_interpolation": True,
                    "no_original_coordinates_or_source_ids_exported": True}
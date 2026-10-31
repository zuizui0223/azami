"""Staged SoilGrids34 fallback runner using verified equivalent native pixel windows.

Intentionally NOT enabled in PR #167 while source-real run #41198399573 is
active; this alternative is for a new separate source-verified run if the
current scalar sampler hits its acquisition resource budget.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analysis.v4.rf_soilgrids_34_predictability import run
from analysis.v4.rf_soilgrids_window_sampler import sample_projected_pixels_windowed


def soil_reader_windowed(url: str, latitude, longitude, *, max_seconds: float = 750.) -> np.ndarray:
    """Read exact native pixels only; return NaN for source mask/extent gaps."""
    import rasterio
    from rasterio.warp import transform

    latitude = np.asarray(latitude, dtype=float)
    longitude = np.asarray(longitude, dtype=float)
    if (latitude.ndim != 1 or longitude.shape != latitude.shape or
            not np.isfinite(latitude).all() or not np.isfinite(longitude).all() or
            (np.abs(latitude) > 90).any() or (np.abs(longitude) > 180).any()):
        raise ValueError('INVALID_PRIVATE_SOURCE_WGS84_SITE_COORDINATES')
    path = '/vsicurl?max_retry=1&retry_delay=1&list_dir=no&url=' + url
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
                      GDAL_HTTP_MAX_RETRY='1', GDAL_HTTP_RETRY_DELAY='1',
                      GDAL_HTTP_TIMEOUT='15', GDAL_HTTP_CONNECTTIMEOUT='9',
                      GDAL_CACHEMAX=64):
        with rasterio.open(path) as source:
            if source.count != 1 or source.crs is None:
                raise ValueError('SOIL_SOURCE_RASTER_BAND_OR_CRS_INVALID')
            xs, ys = transform('EPSG:4326', source.crs,
                               longitude.tolist(), latitude.tolist())
            values, _ = sample_projected_pixels_windowed(source, xs, ys,
                                     window_pixels=256, max_seconds=max_seconds)
            return values


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for item in ('native','receipt','archived-seven','archived-seven-sha',
                 'bio-nine','bio-nine-sha','nine-manifest','terrain-table',
                 'terrain-manifest','five-degree','metadata','output'):
        p.add_argument('--'+item,type=str if item.endswith('sha') else Path,required=True)
    p.add_argument('--trees',type=int,default=80)
    p.add_argument('--repeats',type=int,default=2)
    a=p.parse_args()
    report=run(a.native,a.receipt,a.archived_seven,a.archived_seven_sha,
           a.bio_nine,a.bio_nine_sha,a.nine_manifest,a.terrain_table,
           a.terrain_manifest,a.five_degree,a.metadata,a.output,
           trees=a.trees,repeats=a.repeats,reader=soil_reader_windowed)
    report['soil_source_pixel_sampler']='native_256px_GDAL_windows_source_pixel_equivalent'
    report['source_sample_values_not_interpolated']=True
    a.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],
                      'original_s1':report['source']['s1_original_n'],
                      'n_taxa_fit':report['soil_cell_support']['n_taxa_after_soil_requalification'],
                      'models_run':len(report['results']),
                      'source_reader':'native_windowed'}),flush=True)

if __name__=='__main__':
    main()
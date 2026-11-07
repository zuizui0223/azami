"""Synthetic-only native raster semantics, no original S1 or remote sources."""
import numpy as np
import pytest
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from analysis.v4.rf_soilgrids_window_sampler import sample_projected_pixels_windowed


def test_window_reader_matches_masked_source_pixels_across_dispersed_blocks():
    rng = np.random.default_rng(123)
    arr = rng.integers(1, 5000, size=(641, 513), dtype=np.int16)
    arr[30:52, 70:82] = -9999
    transform = from_origin(-180, 80, 0.05, 0.05)
    with MemoryFile() as f:
        with f.open(driver='GTiff', width=513, height=641, count=1,
                    dtype='int16', nodata=-9999, transform=transform,
                    crs='EPSG:4326', tiled=True, blockxsize=256, blockysize=256) as w:
            w.write(arr,1)
        with f.open() as src:
            rows=rng.integers(0,641,size=2500)
            cols=rng.integers(0,513,size=2500)
            # Also check out of bounds / masking and exact pixel centring
            rowcent=rows+.5; colcent=cols+.5
            xs=transform.c + colcent * transform.a
            ys=transform.f + rowcent * transform.e
            xs=np.r_[xs, -999.,999.]
            ys=np.r_[ys, 999.,-999.]
            expected=[]
            for p in src.sample(zip(xs,ys),indexes=1,masked=True):
                q=p[0]
                expected.append(np.nan if np.ma.is_masked(q) else float(q))
            got,diag=sample_projected_pixels_windowed(src,xs,ys,window_pixels=256)
            assert np.allclose(got,np.array(expected),equal_nan=True)
            assert diag['n_source_points']==2502
            assert diag['n_inside_source_grid']==2500
            assert diag['n_unique_source_windows']<=9
            assert diag['no_interpolation']
            assert 'obs_id' not in repr(diag)


def test_window_reader_rejects_unsafe_grain_and_hard_timeout():
    with MemoryFile() as f:
        with f.open(driver='GTiff',width=16,height=16,count=1,dtype='float32',
                    transform=from_origin(-100,60,1,1),crs='EPSG:4326') as w:
            w.write(np.ones((16,16), dtype=np.float32),1)
        with f.open() as src:
            with pytest.raises(ValueError,match='NONFINITE'):
                sample_projected_pixels_windowed(src,[1.,np.nan],[1.,2.])
            with pytest.raises(ValueError,match='INVALID_FIXED'):
                sample_projected_pixels_windowed(src,[1.],[1.],window_pixels=0)
            with pytest.raises(TimeoutError,match='BUDGET_EXCEEDED'):
                sample_projected_pixels_windowed(src,[-99.5],[59.5],max_seconds=0.1,
                                                 clock=iter([0.,0.5]).__next__)
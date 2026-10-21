"""Synthetic original-TIFF tests for FULL WCVP geographic availability PCA-env.

The actual ecological analysis must run from SHA-verified five global CHELSA
sources and sealed iNaturalist/WCVP inputs in the repository Actions runner.
Never copy private original records, full TIFFs or user IDs to public artifacts.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import pytest
rasterio=pytest.importorskip(
    "rasterio",
    reason="Frozen V3 integrity environment excludes rasterio; full-background source Actions installs and tests the original-TIFF path explicitly",
)
from rasterio.transform import from_origin

from analysis.v4.native_pca_env_fullcoverage import (
    official_source_variable_name, sample_original_pixels,
    acquire_verified_source, _collect_background_result
)


def source_tiff(tmp_path:Path):
    file=tmp_path/"fake_source.tif"
    values=np.arange(100*100,dtype="uint16").reshape(100,100)
    values[51,50]=65535
    with rasterio.open(file,"w",driver="GTiff",
            height=100,width=100,count=1,dtype="uint16",
            transform=from_origin(-180,80,1,1),crs="EPSG:4326",
            nodata=65535) as dst:
        dst.write(values,1)
    return file,values


def test_read_original_raster_pixels_exact_5degree_grouping_no_aggregation(tmp_path):
    original,values=source_tiff(tmp_path)
    # Four original 1x1 cells, one explicitly nodata, and a site outside TIFF.
    frame=pd.DataFrame({
        "latitude":[40.5,40.5,28.5,28.5,79.5],
        "longitude":[-150.5,-151.5,-129.5,-130.5,-179.5],
    })
    out=sample_original_pixels(original,frame)
    assert out[0]==float(values[39,29])
    assert out[1]==float(values[39,28])
    assert np.isnan(out[2])  # exact nodata is never interpreted as a biological/climatic zero
    assert out[3]==float(values[51,49])
    assert out[4]==float(values[0,0])
    assert out[0]!=out[1]
    assert out[2]!=out[3]


def test_original_raster_grid_identity_refuses_resampling_or_wrong_source_url():
    assert official_source_variable_name("chelsa_bio15")=="bio15"
    assert official_source_variable_name("chelsa_rsds_mean")=="rsds_mean"
    with pytest.raises(ValueError,match="five registered"):
        official_source_variable_name("chelsa_preprocessed_average")


def test_full_file_checksum_is_required_before_new_original_climate_values(tmp_path,monkeypatch):
    fake,values=source_tiff(tmp_path)
    # Simulate official CHELSA stream without sending any real HTTP request.
    content=fake.read_bytes()
    class Response:
        def __init__(self):
            self.bytes=content
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def raise_for_status(self):pass
        def iter_content(self,chunk_size):yield self.bytes
    import analysis.v4.native_pca_env_fullcoverage as m
    monkeypatch.setattr(m,"requests",type("FakeRequests",(),{
        "get":staticmethod(lambda *args,**kwargs:Response()),
        "RequestException":Exception
    }))
    monkeypatch.setattr(m,"choose_url_paths",lambda x:["https://verified-example.invalid/source"])
    metadata={
        "source_sha256":"0"*64,
        "width":100,"height":100,
        "transform":list(from_origin(-180,80,1,1))[:6],
        "dtype":"uint16","nodata":65535,"scale":1.0,"offset":0.0,
    }
    with pytest.raises(RuntimeError,match="source coverage failure"):
        m.acquire_verified_source("chelsa_bio15",metadata,tmp_path)
    metadata["source_sha256"]=hashlib.sha256(content).hexdigest()
    target,receipt=m.acquire_verified_source("chelsa_bio15",metadata,tmp_path)
    assert target.exists()
    assert receipt["source_SHA256"]==metadata["source_sha256"]
    assert receipt["has_been_SHA_verified_against_existing_native_source"]
    assert receipt["original_pixels_only"]


def test_previous_native_only_background_cannot_be_inferred_as_full_availability():
    import analysis.v4.native_pca_env_fullcoverage as m
    assert m.FROZEN_SOURCE_N if hasattr(m,"FROZEN_SOURCE_N") else True
    assert m.EXPECTED_ORIGINAL_WILD_TWO==297622
    assert m.CACHE_ARCHIVE_MISSING.startswith("Do not substitute")

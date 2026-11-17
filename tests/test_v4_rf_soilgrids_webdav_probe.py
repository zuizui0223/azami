"""All tests mock official network sources; no S1 and no real soil claimed."""
import hashlib
import json
from pathlib import Path
import pytest
from analysis.v4.rf_soilgrids_webdav_probe import (
    SOIL6, CONVERSIONS, official_url, parsed_vrt_metadata, audit,
)


def fake_vrt(code):
    return (f'<VRTDataset rasterXSize="40" rasterYSize="25">'
            '<SRS>EPSG:152160</SRS>'
            '<GeoTransform>0,250,0,0,0,-250</GeoTransform>'
            '<VRTRasterBand dataType="Int16" band="1">'
            '<SimpleSource><SourceFilename relativeToVRT="1">'
            f'{code}_0-5cm_Q0.5/tileSG-001-001.tif'
            '</SourceFilename></SimpleSource></VRTRasterBand></VRTDataset>').encode()


def fake_fetch(url):
    name=url.split('/')[-1].split('_')[0]
    return fake_vrt(name), {"http_status":200, "content_type":"application/xml"}


def test_six_correct_quantile_depth_urls_and_scales():
    assert len(SOIL6)==len(CONVERSIONS)==6
    assert [CONVERSIONS[n][1] for n in SOIL6]==[10,10,10,10,100,10]
    assert all(official_url(n).endswith('_0-5cm_Q0.5.vrt') for n in SOIL6)
    with pytest.raises(ValueError,match="UNREGISTERED"):
        official_url('soil_nitrogen_0_5cm')


def test_valid_vrt_has_real_georeferenced_tiles():
    v=parsed_vrt_metadata(fake_vrt('phh2o'),'phh2o')
    assert v['raster_tile_references']==1 and v['band_count']==1
    assert v['srs_present'] and v['raw_source_tiles_not_downloaded']


def test_invalid_html_and_wrong_property_fail_closed():
    with pytest.raises(ValueError,match="NON_VRT"):
        parsed_vrt_metadata(b'<html><body>200 but not data</body></html>','phh2o')
    with pytest.raises(ValueError,match="WRONG_PROPERTY"):
        parsed_vrt_metadata(fake_vrt('clay'),'phh2o')
    with pytest.raises(ValueError,match="ENTITY"):
        parsed_vrt_metadata(b'<!DOCTYPE root [<!ENTITY x SYSTEM "file:///etc/passwd">]>'+fake_vrt('phh2o'),'phh2o')


def test_complete_six_metadata_source_probe_no_fabricated_pixel_data():
    result=audit(fetcher=fake_fetch)
    assert result['authentic_vrt_layers']==6
    assert result['status'].startswith('all_six_')
    assert result['public_sentinel']['status']=='not_requested'
    assert result['no_private_source_access']
    assert result['not_an_RF_model_or_soil_extraction']
    assert all(len(x['vrt_file_sha256'])==64 for x in result['layers'].values())
    assert 'obs_id' not in json.dumps(result)


def test_failure_one_layer_not_silently_skipped():
    def missing_1(url):
        if '/soc/' in url: raise ConnectionError('network')
        return fake_fetch(url)
    result=audit(fetcher=missing_1)
    assert result['authentic_vrt_layers']==5
    assert result['layers']['soil_soc_0_5cm']['status']=='unavailable_or_invalid_source'
    assert result['status']=='at_least_one_soil_vrt_unavailable'


def test_public_pixel_sampling_explicit_and_not_original_coordinates():
    result=audit(fetcher=fake_fetch, sampler=lambda u,div: {
        'sentinel_sample_status':'one_original_pixel_decoded',
        'sentinel_physical_value':6.2,'public_point_only':True})
    assert result['public_sentinel']['public_point_only']
    assert result['public_sentinel']['sentinel_physical_value']==6.2
    assert 'latitude' not in json.dumps(result)
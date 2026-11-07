import hashlib,json
from pathlib import Path
import pytest
from analysis.v4.rf_soilgrids_28_to34 import (
    SOIL34,TERRAIN28,SOIL6,PAIRED,SOIL_MANIFEST_UNITS,verify_soil_source_manifest
)
from analysis.v4.rf_soilgrids_webdav_probe import official_url, CONVERSIONS
from analysis.v4.rf_soilgrids_private_pilot import VRT_PINNED_SHA256


def valid_manifest(n=100):
    return {'schema':'azami_rf_extra32_source_manifest_v1',
            'source_observations':n,
            'original_S1_point_pixels_uninterpolated':True,
            'all_used_tiff_tile_byte_sha_verified':True,
            'raw_tiff_tile_receipt_private':True,
            'original_pixel_nodata_preserved':True,
            'layers':{
                name:{'dataset':'SoilGrids-v2.0','temporal_support':'static-soil-model',
                      'depth_interval':'0-5cm','statistic':'Q0.50',
                      'source_url':official_url(name),
                      'source_bundle_sha256':VRT_PINNED_SHA256[name],
                      'unit':SOIL_MANIFEST_UNITS[name],
                      'conversion_factor_divisor':CONVERSIONS[name][1],
                      'source_tiff_tiles_hashed':12,
                      'derivation':'Source pixels at unchanged observation locations',
                      'sampling_method':'native original-pixel, no interpolation'}
                for name in SOIL6}}


def test_model_is_28_plus_six_not_future_unmeasured_soil39():
    assert len(TERRAIN28)==28
    assert len(SOIL34)==34
    assert len(set(SOIL34))==34
    assert set(SOIL6).isdisjoint(TERRAIN28)
    assert len(PAIRED)==2
    assert set(PAIRED['terrain28_plus_soilgrids_0to5cm_Q05_34'])==set(SOIL34)


def test_full_original_tiff_provenance_gate(tmp_path):
    f=tmp_path/'manifest.json';f.write_text(json.dumps(valid_manifest()))
    got=verify_soil_source_manifest(f,source_rows=100)
    assert got['source_rows']==100
    assert got['soil6_source_vrt_pinned']
    assert not got['no_photo_ids_coordinates_or_taxon_names_exported'] is False


@pytest.mark.parametrize('tamper',[
    'wrong_depth','wrong_quantile','wrong_vrt_sha','no_tile_sha',
    'unverified_byte_receipts','wrong_unit','different_population','extra_soil',
])
def test_rejects_invalid_soil_manifest_or_surrogate_data(tmp_path,tamper):
    d=valid_manifest()
    item=d['layers']['soil_phh2o_0_5cm']
    if tamper=='wrong_depth':item['depth_interval']='0-30cm'
    if tamper=='wrong_quantile':item['statistic']='mean'
    if tamper=='wrong_vrt_sha':item['source_bundle_sha256']='f'*64
    if tamper=='no_tile_sha':item['source_tiff_tiles_hashed']=0
    if tamper=='unverified_byte_receipts':d['all_used_tiff_tile_byte_sha_verified']=False
    if tamper=='wrong_unit':item['unit']='mg/kg'
    if tamper=='different_population':d['source_observations']=99
    if tamper=='extra_soil':d['layers']['soil_silt_0_5cm']=item.copy()
    f=tmp_path/'source.json';f.write_text(json.dumps(d))
    with pytest.raises(ValueError, match='SOILGRIDS_UNREGISTERED|SOILGRIDS_UNVERIFIED'):
        verify_soil_source_manifest(f,source_rows=100)
import json

from analysis.v3.workflow import ROOT
from analysis.v3.measurement_chunks import contract


def receipts():
    return [json.loads((ROOT/'reproducibility'/name).read_text(encoding='utf-8')) for name in (
        'v3_partial_raw_stream_observation_views_20260908.json',
        'v3_partial_observation_view_preservation_20260908.json')]


def test_partial_view_keeps_full_source_and_all27_without_ecological_admission():
    view,saved=receipts()
    assert view['source_schedule_sha256']==contract()['source_schedule_sha256']
    assert view['source_observations']==319244
    assert view['complete_request_observations']+view['pending_observations']+view['source_only_blocked_observations']==view['source_observations']
    assert view['operational_routes']+view['held_routes']==view['retained_endpoints']==27
    assert view['local_verification']['all27_inventory_rows']==319244*27
    assert view['local_verification']['held_routes_with_nonnull_operational_values']==0
    assert view['ecological_models_executed']==saved['ecological_models_executed']==0
    assert view['ecological_fitting_authorized'] is saved['ecological_fitting_authorized'] is False


def test_protected_view_roundtrip_binds_the_same_database_not_an_image_archive():
    view,saved=receipts()
    assert saved['source_view_database_sha256']==view['database_sha256']
    assert saved['files_restored']==saved['asset']['files']==35
    assert saved['bytes_restored']==saved['asset']['restored_bytes']
    assert saved['sqlite_integrity_check']=='ok'
    assert saved['retained_observation_endpoint_slots']==view['local_verification']['all27_inventory_rows']
    assert saved['draft_verified'] is True and saved['anonymous_release_and_asset_requests']=='404'
    assert saved['source_images_included'] is saved['asset']['source_images_included'] is False
    assert saved['image_requests']==0

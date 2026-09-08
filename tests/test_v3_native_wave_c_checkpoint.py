import json

from analysis.v3.workflow import ROOT, canonical_digest


def read(path):
    return json.loads((ROOT/path).read_text(encoding='utf-8'))


def test_live_wave_checkpoint_keeps_planned_and_verified_denominators_separate():
    receipt = read('reproducibility/v3_native_raw_wave_c_checkpoint_20260908.json')
    batch = read(receipt['batch_contract'])
    assert receipt['batch_contract_canonical_sha256'] == canonical_digest(batch)
    assert receipt['plan_id'] == batch['plan_id']
    assert receipt['planned_wave_chunks'] == len(batch['chunks']) == 128
    assert receipt['checkpoint_chunks'] == [f'c{i:06d}' for i in range(18,31)]
    assert len(receipt['records']) == receipt['locally_restored_chunks'] == 13
    assert receipt['cloud_success_reconciled_chunks'] == 10
    assert receipt['cloud_failed_confirmation_but_local_complete_chunks'] == 3
    assert receipt['ecological_fitting_authorized'] is False
    assert receipt['environment_values_read'] == receipt['ecological_models_executed'] == 0
    assert receipt['aggregate']['selected_observations'] == 598
    assert receipt['aggregate']['verified_photo_units'] == 823
    assert receipt['aggregate']['detected_heads'] == 1390
    for key,value in receipt['aggregate'].items():
        assert value == sum(row['local_restoration'][key] for row in receipt['records'])
    for row in receipt['records']:
        local = row['local_restoration']
        planned = batch['chunks'][row['chunk_id']]
        assert canonical_digest(local) == row['local_report_canonical_sha256']
        assert local['verified_photo_units'] == planned['requests']
        assert local['selected_observations'] == planned['observations']
        assert local['input_packet_sha256'] == planned['packet_sha256']
        assert local['raw_endpoint_slots'] == 27*local['detected_heads']
        assert local['bbox_slots'] == 25*local['detected_heads']
        assert local['local_image_requests'] == local['source_images_persisted'] == 0
        assert len(local['measurement_eligible_heads']) == 27


def test_failed_confirmation_is_not_changed_to_a_successful_cloud_job():
    receipt = read('reproducibility/v3_native_raw_wave_c_checkpoint_20260908.json')
    failed = [r for r in receipt['records'] if r['job_conclusion']=='failure']
    assert [r['chunk_id'] for r in failed] == ['c000022','c000023','c000028']
    for row in failed:
        assert 'cloud' not in row
        assert row['failure']['stage'] == 'checkpoint_return_download'
        assert row['failure']['error_type'] == 'ReadTimeout'
        assert row['failure']['failed_job_recovery_executed'] is False
        assert row['failure']['all_units_committed_before_failure'] == row['local_restoration']['verified_photo_units']
        assert row['failure']['additional_image_requests_by_local_restoration'] == 0
    for row in receipt['records']:
        if row['job_conclusion']=='success':
            assert row['cloud']['output_asset'] == row['local_restoration']['asset']
            assert row['cloud']['commit'] == receipt['execution_commit']

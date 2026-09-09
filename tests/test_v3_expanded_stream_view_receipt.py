import json

from analysis.v3.stream_observation_views import SPECIFICATION
from analysis.v3.workflow import ROOT, canonical_digest, text_digest


def test_expanded_partial_view_preserves_definitions_denominators_and_prior_photos():
    report = json.loads((ROOT / 'reproducibility/v3_partial_raw_stream_observation_views_20260909.json').read_text())
    qa, snapshot = report['local_verification'], report['local_snapshot']
    original = {k: v for k, v in report.items() if k not in {
        'local_verification', 'local_snapshot', 'input_archives', 'reused_pilot_observations', 'prior_public_view', 'execution_boundary'}}
    assert canonical_digest(original) == qa['view_report_canonical_sha256']
    assert report['specification'] == SPECIFICATION
    assert report['implementation_sha256_text_lf'] == text_digest(ROOT / 'analysis/v3/stream_observation_views.py')
    assert report['source_observations'] == 319244 == sum(report[k] for k in
        ('complete_request_observations', 'pending_observations', 'source_only_blocked_observations'))
    assert qa['all27_inventory_rows'] == 319244 * 27
    assert report['retained_endpoints'] == report['operational_routes'] + report['held_routes'] == 27
    assert qa['held_routes_with_nonnull_operational_values'] == qa['pending_slots_with_nonnull_operational_values'] == 0
    assert qa['earlier_photo_units_preserved_unchanged'] == 800
    assert len(report['input_archives']) == len(set(report['input_archives'])) == qa['archives'] == 54
    assert report['verified_photo_units'] == 3559 and report['detected_heads'] == 6145
    assert snapshot['database_sha256'] == qa['database_sha256'] == report['database_sha256']
    assert snapshot['files_restored'] == snapshot['asset']['files'] == 168
    assert snapshot['cloud_roundtrip_verified'] is False
    assert report['ecological_fitting_authorized'] is False


def test_inventory_completion_does_not_erase_cloud_failures_or_claim_full_payload_replay():
    report = json.loads((ROOT / 'reproducibility/v3_wave_c_inventory_completion_20260909.json').read_text())
    assert report['planned_chunks'] == report['protected_final_chunks'] == 128
    assert report['missing_chunks'] == report['resume_chunks'] == []
    assert report['actions_job_conclusion_used_for_completion'] is False
    assert report['protected_numerical_payload_downloaded'] is False
    assert report['execution']['overall_actions_conclusion'] == 'failure'
    assert [r['chunk'] for r in report['execution']['missing_only_resumes']] == ['c000036', 'c000111']
    assert report['execution']['all_other_chunks_reexecuted'] is False
    assert report['ecological_fitting_authorized'] is False


def test_later_cloud_replay_matches_the_exact_expanded_view_without_rewriting_its_local_checkpoint():
    view = json.loads((ROOT / 'reproducibility/v3_partial_raw_stream_observation_views_20260909.json').read_text())
    replay = json.loads((ROOT / 'reproducibility/v3_expanded_stream_view_preservation_20260909.json').read_text())
    assert replay['source_view_report_canonical_sha256'] == view['local_verification']['view_report_canonical_sha256']
    assert replay['source_view_database_sha256'] == view['database_sha256']
    assert replay['verified_photo_units'] == view['verified_photo_units'] == 3559
    assert replay['files_restored'] == view['local_snapshot']['files_restored'] == 168
    assert replay['bytes_restored'] == view['local_snapshot']['bytes_restored'] == 350698483
    for key, value in view['local_snapshot']['asset'].items():
        assert replay['asset'][key] == value
    assert replay['draft_verified'] is True and replay['anonymous_release_and_asset_requests'] == '404'
    assert replay['sqlite_integrity_check'] == 'ok'
    assert replay['retained_observation_endpoint_slots'] == 319244 * 27
    assert replay['source_images_included'] is False and replay['ecological_fitting_authorized'] is False

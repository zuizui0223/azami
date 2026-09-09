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

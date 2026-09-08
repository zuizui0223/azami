import json

from analysis.v3.uniform_colour_profile import RECIPE, METRICS
from analysis.v3.workflow import ROOT, canonical_digest, text_digest


def test_executed_candidate_profile_keeps_complete_grid_and_no_ecology():
    report = json.loads((ROOT/'reproducibility/v3_uniform_colour_technical_profile_20260908.json').read_text())
    recipe = json.loads(RECIPE.read_text())
    assert report['execution_contract']['recipe_canonical_sha256'] == canonical_digest(recipe)
    assert report['execution_contract']['implementation_sha256_text_lf'] == text_digest(ROOT/'analysis/v3/uniform_colour_profile.py')
    assert report['scheduled_heads'] == 2853 and report['condition_rows'] == 2853*14
    assert report['summary_rows'] == len(report['records']) == 14*len(METRICS)*3
    keys = {(r['condition'],r['metric_id'],r['exposure_stratum']) for r in report['records']}
    assert len(keys) == len(report['records'])
    for r in report['records']:
        assert sum(r[k] for k in ('paired_usable_heads','lost_usable_heads','gained_usable_heads','neither_usable_heads')) == r['scheduled_heads']
    assert not report['ecological_fitting_authorized'] and not report['independent_accuracy_estimated']
    assert report['ecological_models_executed'] == report['environment_values_read'] == report['image_requests'] == 0


def test_candidate_profile_protected_copy_binds_exact_numerical_report():
    report = json.loads((ROOT/'reproducibility/v3_uniform_colour_technical_profile_20260908.json').read_text())
    saved = json.loads((ROOT/'reproducibility/v3_uniform_colour_profile_preservation_20260908.json').read_text())
    assert saved['source_profile_report_canonical_sha256'] == canonical_digest(report)
    assert saved['verified_condition_rows'] == report['condition_rows']
    assert saved['files_restored'] == saved['asset']['files'] == 4
    assert saved['bytes_restored'] == saved['asset']['restored_bytes']
    assert saved['draft_verified'] and saved['anonymous_release_and_asset_requests'] == '404'
    assert not saved['source_images_included'] and not saved['ecological_fitting_authorized']

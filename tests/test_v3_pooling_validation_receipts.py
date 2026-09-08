import json
from pathlib import Path

from analysis.v3.workflow import canonical_digest, text_digest

ROOT = Path(__file__).resolve().parents[1]


def read(relative):
    return json.loads((ROOT / relative).read_text(encoding='utf-8'))


def evidence():
    receipt = read('reproducibility/v3_joint_pooling_validation_20260908.json')
    reports = {key: read(value['path']) for key,value in receipt['evidence'].items()}
    return receipt, reports


def test_synthetic_receipts_keep_exact_content_identity():
    receipt, reports = evidence()
    for key,value in receipt['evidence'].items():
        assert canonical_digest(reports[key]) == value['canonical_json_sha256']
    assert receipt['ecological_fitting_authorized'] is False
    assert receipt['ecological_models_executed'] == 0
    assert receipt['empirical_trait_environment_values_read'] == 0


def test_score_correction_preserves_the_entire_simulation_design_and_failures():
    _, reports = evidence()
    before, after = reports['initial_simulation'], reports['scorefix_simulation']
    assert before['execution_contract']['specification'] == after['execution_contract']['specification']
    for report in (before,after):
        assert report['replicate_records'] == 360
        assert report['ecological_fitting_authorized'] is False
        assert [row['zero_rejections'] for row in report['summaries']] == [4,12,8]
        for row in report['summaries']:
            assert row['planned_replicates'] == row['estimable_replicates'] + row['unestimable_replicates'] == 120
            assert row['all_taxa_retained_in_every_estimable_replicate'] is True
        assert report['summaries'][1]['coverage_among_estimable'] == .9


def test_scorefix_historical_report_preserves_its_versioned_implementation():
    receipt, reports = evidence()
    execution = reports['scorefix_simulation']['execution_contract']
    assert execution['implementation_sha256_text_lf'] == '13851f324bfb1d31602bb824d552883cca20e394efbc1724154bd91cdb865fac'
    assert execution['simulation_sha256_text_lf'] == text_digest(ROOT/'analysis/v3/simulate_joint_pooling.py')
    assert execution['pooling_contract_canonical_sha256'] == canonical_digest(read('analysis/v3/joint_partial_pooling_contract.json'))
    assert receipt['fullscale_execution']['implementation_sha256_text_lf'] == execution['implementation_sha256_text_lf']
    generator = receipt['fullscale_generator_reference']
    assert text_digest(ROOT/generator['path']) == generator['sha256_text_lf']


def test_current_contraction_reordering_preserves_model_and_simulation_design():
    _, reports = evidence()
    previous = reports['scorefix_simulation']
    current = read('reproducibility/v3_joint_pooling_synthetic_contractions_20260908.json')
    check = read('reproducibility/v3_pooling_contraction_equivalence_20260908.json')
    assert current['execution_contract']['implementation_sha256_text_lf'] == text_digest(ROOT/'analysis/v3/joint_partial_pooling.py') == check['implementation_sha256_text_lf']
    assert check['reference_sha256_text_lf'] == previous['execution_contract']['implementation_sha256_text_lf']
    assert current['execution_contract']['specification'] == previous['execution_contract']['specification']
    assert current['execution_contract']['pooling_contract_canonical_sha256'] == canonical_digest(read('analysis/v3/joint_partial_pooling_contract.json'))
    assert current['replicate_records'] == 360
    assert [r['zero_rejections'] for r in current['summaries']] == [4,12,8]
    assert current['summaries'][1]['coverage_among_estimable'] == .9
    assert all(r['estimable_replicates']==120 and r['all_taxa_retained_in_every_estimable_replicate'] for r in current['summaries'])
    for row in check['differences']:
        assert row['hypermean_max_absolute_difference'] < 1e-7
        assert row['tau2_max_absolute_difference'] < 1e-6
        assert row['joint_prediction_covariance_max_absolute_difference'] < 1e-6
    assert current['ecological_fitting_authorized'] is check['ecological_fitting_authorized'] is False


def test_paired_crossed_replay_keeps_one_dataset_and_every_planned_draw():
    replay = read('reproducibility/v3_crossed_bootstrap_paired_replay_20260908.json')
    assert replay['independent_generated_datasets']==1 and replay['paired_replays']==2
    assert replay['ecological_fitting_authorized'] is False
    for row in replay['comparisons']:
        assert row['paired_draws']==199
        assert row['candidate_probabilities_and_coverage_decisions_identical'] is True
        assert all(value <= replay['absolute_tolerance'] for name,value in row.items() if name.endswith('max_absolute_difference'))
    execution = replay['evidence'][1]['execution_contract']
    for path,sha in execution['implementation_sha256_text_lf'].items():
        assert sha==text_digest(ROOT/path)
    assert execution['resampling_contract_canonical_sha256']==canonical_digest(read('analysis/v3/dependence_resampling_contract.json'))


def test_fullscale_replay_uses_original_score_tolerance_without_losing_taxa():
    _, reports = evidence()
    report = reports['scorefix_fullscale']
    assert report['observations'] == 319244
    assert report['retained_taxa'] == report['taxa'] == 354
    assert report['prediction_covariance_shape'] == [3186,3186]
    assert report['ecological_fitting_authorized'] is False
    for attempt in report['optimizer_attempts']:
        assert attempt['accepted'] and attempt['stationary']
        assert attempt['projected_gradient_max'] <= 1e-4
        assert not attempt['upper_bound_hit']
        polishing = attempt['stationarity_polishing']
        assert polishing['accepted'] and polishing['minimum_free_curvature'] > 0

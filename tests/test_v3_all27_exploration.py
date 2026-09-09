import csv
import io
import json

import numpy as np
import pandas as pd
import pytest

from analysis.v3 import all27_exploration as a
from analysis.v3 import cloud_all27_exploration as cloud
from analysis.v3.build_observation_measurements import HUE, COMPOSITION


KEYS = [r['endpoint_id'] for r in a.registry()]
HELD = 'bract_projection_roughness'


def numerical_bytes():
    diagnostics, rows = [], []
    for i in range(2):
        diagnostics.append(json.dumps({'head_index': i, 'diagnostics': {
            'head_min_dimension_px': 200 + i * 100, 'head_laplacian_variance': 80 + i * 20}}))
        for key in KEYS:
            rows.append({'head_index': i, 'endpoint_id': key,
                         'value': .25 if key in COMPOSITION else 1 + 8 * i,
                         'status': 'bad_qc' if key == HELD and i else 'usable'})
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
    return ('\n'.join(diagnostics)).encode(), handle.getvalue().encode()


def test_held_coordinates_recover_same_qc_head_quality_without_promotion():
    rows = a.head_summaries(*numerical_bytes())
    assert len(rows) == 27
    held = next(r for r in rows if r['endpoint_id'] == HELD)
    assert held['raw_mean'] == 5 and held['eligible_mean'] == 1
    assert held['head_min_dimension_mean_px'] == 200
    assert held['head_sharpness_mean'] == 80
    assert held['n_eligible_heads'] == 1
    spec, env = a.definition()
    assert spec['measurement_holds_overridden'] is False
    assert spec['planned_coefficient_slots'] == 27 * len(env['variables']) * len(a.SCALES)


def test_joint_hue_and_composition_support_is_not_bypassed():
    diag, raw = numerical_bytes()
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    for row in rows:
        if row['endpoint_id'] in (HUE[0], COMPOSITION[0]):
            row['status'] = 'bad_qc'
    out = io.StringIO(); writer = csv.DictWriter(out, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
    summaries = a.head_summaries(diag, out.getvalue().encode())
    assert all(r['eligible_mean'] is None for r in summaries if r['endpoint_id'] in HUE + COMPOSITION)
    assert len(a.head_summaries(b'', b'head_index,endpoint_id,value,status\n')) == 27


def quality_fixture():
    photos = pd.DataFrame([{'photo_id': p, 'endpoint_id': k, 'value': v, 'size': s, 'sharpness': 100.}
                           for p, v, s in [('p1', 10., 200.), ('p2', 30., 300.)] for k in KEYS])
    images = pd.DataFrame([('1', 'pixel1', 'p1'), ('1', 'pixel2', 'p2'), ('2', 'pixel1', 'p1')],
                          columns=['obs_id', 'pixel_id', 'photo_id'])
    support = pd.DataFrame([('1', 0, 2), ('2', 1, 2)], columns=['obs_id', 'n_pending', 'n_requests'])
    return photos, images, support


def test_equal_images_endpoint_specific_missingness_and_pending():
    photos, images, support = quality_fixture()
    photos.loc[photos.photo_id.eq('p2') & photos.endpoint_id.eq(HELD), 'size'] = np.nan
    out = a.aggregate_quality(photos, images, support)
    row = out[out.obs_id.eq('1') & out.endpoint_id.eq(HELD)].iloc[0]
    assert row.value == 20 and np.isnan(row['size']) and row.sharpness == 100
    assert out[out.obs_id.eq('2')].value.isna().all()
    assert len(out) == 54
    # Missing values are not zeros and do not borrow quality from other heads/photos.
    photos.loc[photos.photo_id.eq('p2') & photos.endpoint_id.eq(HELD), ['value', 'size', 'sharpness']] = np.nan
    out = a.aggregate_quality(photos, images, support)
    row = out[out.obs_id.eq('1') & out.endpoint_id.eq(HELD)].iloc[0]
    assert row.value == 10 and row['size'] == 200 and row.images == 1


@pytest.mark.parametrize('damage', ['duplicate_pixel', 'missing_endpoint', 'unmatched_quality'])
def test_invalid_measurement_joins_stop(damage):
    photos, images, support = quality_fixture()
    if damage == 'duplicate_pixel':
        images = pd.concat([images, images.iloc[:1]])
    elif damage == 'missing_endpoint':
        photos = photos.iloc[1:]
    else:
        photos.loc[0, 'value'] = np.nan
    with pytest.raises(ValueError):
        a.aggregate_quality(photos, images, support)


def fit_frame(key=HELD):
    rng = np.random.default_rng(55931)
    _, env = a.definition()
    n = 300
    f = pd.DataFrame(rng.normal(size=(n, 9)), columns=env['variables'])
    f['obs_id'] = [str(i) for i in range(n)]
    f['accepted_key'] = np.repeat([str(i) for i in range(30)], 10)
    f['endpoint_id'] = key
    f['value'] = rng.uniform(size=n)
    f['size'] = rng.uniform(200, 400, n)
    f['sharpness'] = rng.uniform(90, 300, n)
    f['n_pending'], f['n_requests'] = 0, 1
    f['sin_doy'], f['cos_doy'] = np.sin(np.arange(n)), np.cos(np.arange(n))
    f['observed_year'] = rng.integers(2005, 2027, n)
    f['analysis_latitude'], f['analysis_longitude'] = rng.uniform(20, 60, n), rng.uniform(-160, 160, n)
    return f


def test_all_nine_predictors_fixed_nuisance_and_no_formal_inference():
    f = fit_frame(COMPOSITION[0])
    f.loc[0, 'size'] = np.nan
    seen = []
    def fit(y, x, taxa, b):
        seen.append((y, x, taxa, b))
        return {'coefficients': np.ones((3, 1, 9)), 'diagnostics': []}
    result, data = a.fit_endpoint(f, fit=fit)
    assert len(data) == 299 and result['observations'] == 299
    y, x, taxa, b = seen[0]
    np.testing.assert_allclose(y[:, 0], np.sqrt(data.value))
    np.testing.assert_allclose(x, a.transform(data))
    assert b.shape == (299, 13) and len(set(taxa)) == 30
    assert sum(map(len, result['coefficients'].values())) == 27
    assert result['inferential_status'] == 'exploratory_point_only_no_p_values_or_intervals'
    assert result['design']['observations'] == 299


def test_fit_failure_and_zero_support_keep_every_planned_slot():
    def failure(*args):
        raise RuntimeError('Optimizer not converged; retain failure')
    result, _ = a.fit_endpoint(fit_frame(), fit=failure)
    assert result['status'] == 'not_estimable' and 'not converged' in result['reason']
    assert all(v is None for row in result['coefficients'].values() for v in row.values())
    f = fit_frame(); f['n_pending'] = 1
    result, data = a.fit_endpoint(f, fit=lambda *args: pytest.fail('No data must not be fitted'))
    assert data.empty and result['observations'] == 0 and result['reason']


def test_run_preserves_all729_slots_and_private_membership(tmp_path, monkeypatch):
    prepared, out = tmp_path / 'prepared', tmp_path / 'result'
    prepared.mkdir()
    source = fit_frame().iloc[:2].copy()
    source['native_range_status'] = 'native'
    _, env = a.definition()
    source[['obs_id', *a.SOURCE_FIELDS, *env['variables']]].to_csv(prepared / 'environment_source_private.csv', index=False)
    measures = []
    for key in KEYS:
        f = source[['obs_id', 'value', 'size', 'sharpness', 'n_pending', 'n_requests']].copy()
        f['endpoint_id'], f['component_id'], f['derived_component'] = key, 'private-source', 'private-derived'
        measures.append(f)
    pd.concat(measures).to_csv(prepared / 'all27_observation_qc_private.csv', index=False)
    a.new_json(prepared / 'preparation_report.json', {
        'contract_canonical_sha256': a.canonical_digest(a.definition()[0]),
        'input_sha256': {p.name: a.digest(p) for p in prepared.iterdir()}})
    calls = []
    def fit(frame):
        calls.append(frame.endpoint_id.iloc[0])
        r = {'observations': 2, 'taxa': 1, 'status': 'not_estimable',
             'reason': 'Synthetic support deliberately insufficient',
             'coefficients': {s: dict.fromkeys(env['variables']) for s in a.SCALES}}
        return r, frame
    monkeypatch.setattr(a, 'fit_endpoint', fit)
    result = a.run(prepared, out)
    assert calls == KEYS and result['endpoints_attempted'] == 27
    assert result['endpoints_estimated'] == 0 and result['coefficient_slots_retained'] == 729
    assert result['primary_models_executed'] == result['p_values_reported'] == 0
    coefficients = pd.read_csv(out / 'all729_exploratory_coefficients.csv')
    assert len(coefficients) == 729 and coefficients.estimate.isna().all()
    assert coefficients.measurement_hold_unchanged.sum() == 13 * 27
    assert len(pd.read_csv(out / 'cohort_membership_private.csv')) == 54
    assert 'private-source' not in (out / 'public_report.json').read_text()
    with pytest.raises(ValueError, match='Preserve'):
        a.run(prepared, out)


def test_cloud_batch_rejects_changed_code_or_scope():
    files = ['all27_exploration.py', 'cloud_all27_exploration.py', 'module_ecology.py',
             'joint_partial_pooling.py', 'environment_model.py', 'nuisance_design.py']
    batch = {'status': 'single_partial_checkpoint_all27_point_exploration',
             'exploratory_contract_canonical_sha256': a.canonical_digest(a.definition()[0]),
             'image_requests_authorized': False, 'primary_inference_authorized': False,
             'source_asset': {'source_images_included': False},
             'implementation_sha256_text_lf': {f'analysis/v3/{f}': a.text_digest(a.ROOT / 'analysis/v3' / f) for f in files}}
    assert cloud.validate_batch(batch) is batch
    for key in ('image_requests_authorized', 'primary_inference_authorized'):
        with pytest.raises(ValueError, match='cannot authorize'):
            cloud.validate_batch({**batch, key: True})
    damaged = {**batch, 'implementation_sha256_text_lf': {k: '0' * 64 for k in batch['implementation_sha256_text_lf']}}
    with pytest.raises(ValueError, match='implementation changed'):
        cloud.validate_batch(damaged)


def test_point_design_record_serializes_finite_or_null_values(tmp_path):
    f = fit_frame()
    result, _ = a.fit_endpoint(f, fit=lambda *args: {'coefficients': np.zeros((3, 1, 9)), 'diagnostics': []})
    a.new_json(tmp_path / 'point.json', result)
    assert json.loads((tmp_path / 'point.json').read_text())['observations'] == 300

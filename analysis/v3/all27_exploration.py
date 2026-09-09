"""All-27 basic-QC exploration; never changes primary qualification or tests."""
from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

from .archived_measurement import NumericalSnapshot
from .build_observation_measurements import aggregate_heads, COMPOSITION
from .environment_model import definition as environment_definition, transform
from .image_features import registry
from .model_design_diagnostics import diagnose
from .module_ecology import fit_matched_module, SCALES
from .nuisance_design import matrix as nuisance_matrix
from .protected_artifacts import new_json, require
from .recover_native_source_authority import private_directory
from .stream_observation_views import read_photos
from .workflow import ROOT, canonical_digest, digest, text_digest

CONTRACT = ROOT / 'analysis/v3/all27_exploration_contract.json'
SOURCE_FIELDS = ['accepted_key', 'native_range_status', 'analysis_latitude', 'analysis_longitude',
                 'sin_doy', 'cos_doy', 'observed_year']


def definition():
    spec = json.loads(CONTRACT.read_text(encoding='utf-8'))
    env = environment_definition()
    require(spec['exploratory_fitting_authorized'] is True and spec['measurement_holds_overridden'] is False
            and spec['confirmatory_ecological_fitting_authorized'] is False, 'Exploratory boundary differs')
    require(len(registry()) * len(env['variables']) * len(SCALES) == spec['planned_coefficient_slots'] == 729,
            'Complete all27 coefficient inventory differs')
    return spec, env


def head_summaries(diagnostics: bytes, measurements: bytes):
    """Read numerical bytes only; use the existing QC/quality aggregator unchanged."""
    heads = {}
    for line in diagnostics.decode('utf-8').splitlines():
        row = json.loads(line)
        key = str(row['head_index'])
        require(key not in heads, 'Duplicate head diagnostics')
        heads[key] = {'endpoints': {}, 'diagnostics': row['diagnostics']}
    for row in csv.DictReader(io.StringIO(measurements.decode('utf-8'))):
        head, endpoint = row['head_index'], row['endpoint_id']
        require(head in heads and endpoint not in heads[head]['endpoints'], 'Unmatched or duplicate endpoint')
        value = float(row['value']) if row['value'] else None
        heads[head]['endpoints'][endpoint] = {'value': value, 'status': row['status']}
    return aggregate_heads([heads[k] for k in sorted(heads, key=int)])[0]


def photo_quality(inputs: dict, db):
    """Recover matched QC covariates and reconcile numerical means to the verified view."""
    records, seen = [], {}
    def add(photo, rows):
        identity = canonical_digest(rows)
        if photo in seen:
            require(seen[photo] == identity, 'Repeated photo QC summaries conflict')
            return
        seen[photo] = identity
        old = {r[0]: r[1:] for r in db.execute('SELECT endpoint_id,qc_mean,n_qc FROM photo_endpoints WHERE photo_id=?', (photo,))}
        require(len(old) == len(rows) == 27, 'Photo absent from independently verified view')
        for row in rows:
            key, value = row['endpoint_id'], row['eligible_mean']
            old_value, old_count = old[key]
            require(old_count == row['n_eligible_heads'] and
                    ((old_value is None and value is None) or
                     (old_value is not None and value is not None and np.isclose(old_value, value, rtol=1e-12, atol=1e-12))),
                    'Recovered basic-QC mean differs from verified view')
            records.append({'photo_id': photo, 'endpoint_id': key, 'value': value,
                            'size': row['head_min_dimension_mean_px'], 'sharpness': row['head_sharpness_mean']})
    for photo, row in read_photos(Path(inputs['pilot'])).items():
        add(photo, aggregate_heads(list(row['heads'].values()))[0])
    for entry in inputs['archives']:
        archived = NumericalSnapshot(Path(entry['bundle']), entry['asset'])
        try:
            packet = archived.json('chunk_packet_private.json')
            require(packet['chunk_id'] == entry['chunk_id'], 'Raw numerical archive identity differs')
            for index, item in enumerate(packet['queue']):
                prefix = f'units/u{index:04d}/'
                add(item['photo_id'], head_summaries(archived.read(prefix + 'head_diagnostics_private.jsonl'),
                                                   archived.read(prefix + 'endpoint_measurements_private.csv')))
        finally:
            archived.close()
        print(json.dumps({'qc_archive_recovered': entry['chunk_id'], 'photo_units': len(seen)}), flush=True)
    require(len(seen) == db.execute('SELECT COUNT(*) FROM photo_results').fetchone()[0], 'Not all verified photos recovered')
    return pd.DataFrame(records)


def aggregate_quality(photos, images, support):
    require(not photos.duplicated(['photo_id', 'endpoint_id']).any(), 'Duplicate photo-endpoint')
    require(not images.duplicated(['obs_id', 'pixel_id']).any(), 'Duplicate observation pixel')
    require(not support.obs_id.duplicated().any(), 'Duplicate observation support')
    require(set(photos.endpoint_id) == {r['endpoint_id'] for r in registry()}
            and photos.groupby('photo_id').size().eq(27).all(), 'Incomplete photo endpoint inventory')
    require(not (photos.value.isna() & photos[['size', 'sharpness']].notna().any(axis=1)).any(),
            'Quality exists without a contributing endpoint value')
    require(set(images.obs_id) <= set(support.obs_id), 'Image absent from source support')
    # Intentional image-by-endpoint expansion; one representative per exact pixel.
    x = images.merge(photos, on='photo_id', how='left', validate='many_to_many')
    require(len(x) == len(images) * 27, 'Image-endpoint join changed the registered denominator')
    group = x.groupby(['obs_id', 'endpoint_id'], sort=True)
    result = group.agg(value=('value', 'mean'), size=('size', 'mean'), sharpness=('sharpness', 'mean'),
                       images=('value', 'count'), size_images=('size', 'count'), sharpness_images=('sharpness', 'count')).reset_index()
    for key in ['size', 'sharpness']:
        result.loc[result[key + '_images'].ne(result.images), key] = np.nan
    result = result.merge(support, on='obs_id', validate='many_to_one')
    result.loc[result.n_pending.ne(0) | result.n_requests.eq(0), ['value', 'size', 'sharpness']] = np.nan
    return result.drop(columns=['size_images', 'sharpness_images'])


def prepare(view: Path, inputs: Path, environment: Path, out: Path):
    spec, env = definition()
    out = private_directory(out)
    require(not out.exists(), 'Preserve earlier preparation')
    report = json.loads((view / 'public_report.json').read_text())
    qa = json.loads((view / 'local_verification.json').read_text())
    require(qa['view_report_canonical_sha256'] == canonical_digest(report)
            and qa['input_list_sha256'] == digest(inputs), 'Unverified observation checkpoint inputs')
    database = view / 'observation_views_private.sqlite'
    require(digest(database) == report['database_sha256'] and digest(environment) == env['source_qc_matrix_sha256'],
            'Verified measurement or environment bytes changed')
    out.mkdir(parents=True)
    with sqlite3.connect(database.resolve().as_uri() + '?mode=ro', uri=True) as db:
        photos = photo_quality(json.loads(inputs.read_text()), db)
        images = pd.read_sql_query('SELECT * FROM observation_images', db)
        support = pd.read_sql_query('SELECT obs_id,component_id,n_pending,n_requests FROM observation_support', db)
        groups = pd.read_sql_query('SELECT source_component AS component_id,derived_component FROM derived_components', db)
        support = support.merge(groups, on='component_id', validate='many_to_one')
        quality = aggregate_quality(photos, images, support)
    source = pd.read_csv(environment, usecols=['obs_id', *SOURCE_FIELDS, *env['variables']], dtype={'obs_id': str, 'accepted_key': str})
    require(len(source) == report['source_observations'] and not source.obs_id.duplicated().any()
            and set(source.obs_id) == set(support.obs_id) and source.native_range_status.eq('native').all(), 'Native source identity differs')
    quality.to_csv(out / 'all27_observation_qc_private.csv', index=False)
    source.to_csv(out / 'environment_source_private.csv', index=False)
    receipt = {'status': 'ALL27_EXPLORATORY_INPUT_PREPARED_NO_FITS', 'contract_canonical_sha256': canonical_digest(spec),
               'source_view_database_sha256': report['database_sha256'], 'source_input_list_sha256': digest(inputs),
               'environment_source_matrix_sha256': digest(environment), 'source_observations': len(source),
               'source_taxa': int(source.accepted_key.nunique()), 'verified_photo_units': report['verified_photo_units'],
               'complete_request_observations': report['complete_request_observations'], 'endpoint_slots': 27,
               'input_sha256': {name: digest(out / name) for name in ['all27_observation_qc_private.csv', 'environment_source_private.csv']},
               'image_requests': 0, 'ecological_models_executed': 0, 'confirmatory_ecological_fitting_authorized': False}
    new_json(out / 'preparation_report.json', receipt)
    return receipt


def fit_endpoint(frame, *, fit=fit_matched_module):
    """One simultaneous nine-predictor fit; retain an explicit failed result."""
    _, env = definition()
    require(frame.endpoint_id.nunique() <= 1 and not frame.obs_id.duplicated().any(), 'Mixed endpoint or duplicate observation')
    require(frame.accepted_key.notna().all(), 'Missing taxon identity')
    columns = ['value', 'size', 'sharpness', *env['variables'], 'sin_doy', 'cos_doy', 'observed_year', 'analysis_latitude', 'analysis_longitude']
    numeric = frame[columns].to_numpy(float)
    eligible = np.isfinite(numeric).all(axis=1) & frame['size'].gt(0).to_numpy() & frame.sharpness.ge(0).to_numpy()
    eligible &= frame.n_pending.eq(0).to_numpy() & frame.n_requests.gt(0).to_numpy()
    data = frame.loc[eligible].copy()
    result = {'observations': len(data), 'taxa': int(data.accepted_key.nunique()),
              'excluded_from_supplied_endpoint_rows': len(frame) - len(data), 'status': 'not_estimable',
              'coefficients': {scale: {v: None for v in env['variables']} for scale in SCALES},
              'inferential_status': 'exploratory_point_only_no_p_values_or_intervals'}
    if data.empty:
        result['reason'] = 'No complete endpoint-matched QC and environment observations'
        return result, data
    y = data.value.to_numpy(float)
    if data.endpoint_id.iloc[0] in COMPOSITION:
        require(np.all((y >= 0) & (y <= 1)), 'Invalid QC colour composition')
        y = np.sqrt(y)
    b, names = nuisance_matrix(sin_doy=data.sin_doy, cos_doy=data.cos_doy, observed_year=data.observed_year,
                              latitude=data.analysis_latitude, longitude=data.analysis_longitude,
                              size=data['size'], sharpness=data.sharpness)
    x = transform(data)
    design = data[['obs_id', 'accepted_key']].copy()
    design[env['variables']] = x
    design[list(names)] = b
    design['unit_weight'] = 1.
    result['design'] = diagnose(design, env['variables'], list(names), weight_column='unit_weight', cohort_id=str(data.endpoint_id.iloc[0]))
    try:
        fitted = fit(y[:, None], x, data.accepted_key.to_numpy(), b)
        coefficients = np.asarray(fitted['coefficients'])
        require(coefficients.shape == (3, 1, 9) and np.isfinite(coefficients).all(), 'Incomplete exploratory coefficient vector')
        result.update(status='exploratory_point_estimated_not_qualified',
                      coefficients={scale: dict(zip(env['variables'], coefficients[i, 0].tolist())) for i, scale in enumerate(SCALES)},
                      fit_diagnostics=fitted['diagnostics'])
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        result.update(reason=str(error), error_type=type(error).__name__, optimizer_attempts=getattr(error, 'optimizer_attempts', []))
    return result, data


def run(prepared: Path, out: Path):
    spec, env = definition()
    receipt = json.loads((prepared / 'preparation_report.json').read_text())
    require(receipt['contract_canonical_sha256'] == canonical_digest(spec), 'Prepared exploratory specification differs')
    for name, sha in receipt['input_sha256'].items():
        require(digest(prepared / name) == sha, 'Prepared numerical input changed')
    out = private_directory(out)
    require(not out.exists(), 'Preserve previous exploratory results')
    out.mkdir(parents=True)
    new_json(out / 'execution_contract.json', {'contract': spec, 'contract_canonical_sha256': canonical_digest(spec),
             'preparation_report_sha256': digest(prepared / 'preparation_report.json'),
             'implementation_sha256_text_lf': text_digest(Path(__file__))})
    measurements = pd.read_csv(prepared / 'all27_observation_qc_private.csv', dtype={'obs_id': str})
    source = pd.read_csv(prepared / 'environment_source_private.csv', dtype={'obs_id': str, 'accepted_key': str})
    require(not measurements.duplicated(['obs_id', 'endpoint_id']).any() and not source.obs_id.duplicated().any(), 'Duplicated analysis unit')
    require(source.native_range_status.eq('native').all(), 'Non-native source in exploratory input')
    require(set(measurements.endpoint_id) == {r['endpoint_id'] for r in registry()}, 'All27 scope incomplete')
    joined = measurements.merge(source, on='obs_id', how='left', validate='many_to_one', indicator=True)
    require(joined._merge.eq('both').all(), 'Unmatched environment identity')
    decision = {r['endpoint_id']: r for r in json.loads((ROOT / 'analysis/v3/measurement_qualification_decision_20260908.json').read_text())['endpoints']}
    results, membership = [], []
    for endpoint in registry():
        key = endpoint['endpoint_id']
        result, data = fit_endpoint(joined[joined.endpoint_id.eq(key)])
        result.update(endpoint_id=key, unit=endpoint['unit'], module=decision[key]['module'],
                      original_measurement_route=decision[key]['ecological_route'],
                      measurement_hold_unchanged=decision[key]['ecological_route'] != 'stream_original_required',
                      response_transform='square_root' if key in COMPOSITION else 'identity')
        new_json(out / f'{key}.json', result)
        membership.append(data[['obs_id', 'endpoint_id', 'component_id', 'derived_component']])
        results.append(result)
        print(json.dumps({k: result[k] for k in ['endpoint_id', 'status', 'observations', 'taxa']}), flush=True)
    pd.concat(membership, ignore_index=True).to_csv(out / 'cohort_membership_private.csv', index=False)
    coefficients = [{'endpoint_id': r['endpoint_id'], 'scale': scale, 'predictor': v, 'estimate': r['coefficients'][scale][v],
                     'status': r['status'], 'measurement_hold_unchanged': r['measurement_hold_unchanged']}
                    for r in results for scale in SCALES for v in env['variables']]
    require(len(coefficients) == spec['planned_coefficient_slots'], 'Planned result slots disappeared')
    pd.DataFrame(coefficients).to_csv(out / 'all729_exploratory_coefficients.csv', index=False)
    summary = {'status': 'ALL27_EXPLORATORY_ATTEMPTS_COMPLETE_NOT_CONFIRMATORY',
               'contract_canonical_sha256': canonical_digest(spec), 'endpoints_attempted': len(results),
               'endpoints_estimated': sum(r['status'] == 'exploratory_point_estimated_not_qualified' for r in results),
               'coefficient_slots_retained': len(coefficients), 'measurement_holds_unchanged': 13,
               'image_requests': 0, 'p_values_reported': 0, 'confidence_intervals_reported': 0,
               'primary_models_executed': 0, 'confirmatory_ecological_fitting_authorized': False,
               'input_scope': receipt, 'output_sha256': {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()}}
    new_json(out / 'public_report.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('prepare')
    for arg in ['view', 'inputs', 'environment', 'out']:
        p.add_argument('--' + arg, type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--prepared', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = vars(parser.parse_args())
    action = args.pop('action')
    print(json.dumps(prepare(**args) if action == 'prepare' else run(**args)))


if __name__ == '__main__':
    main()

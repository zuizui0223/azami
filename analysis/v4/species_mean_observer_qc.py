"""Observer-disjoint taxon MEANS, not within-taxon geography or a new RF.

Requires already-authorized original S1 and observer roster; never obtains keys,
recovers encrypted sources, calls a network, or exports taxon/observer records.
Reproducibility of photographic means is not biological measurement validity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

FEATURES = {
    'angle': ['orientation_image_vertical_angle'],
    'lightness': ['corolla_lab_lightness'],
    'chroma': ['corolla_lab_chroma'],
    'hue': ['corolla_hue_sin', 'corolla_hue_cos'],
    'aspect': ['capitulum_outline_aspect_ratio'],
    'circularity': ['capitulum_outline_circularity'],
    'solidity': ['capitulum_outline_solidity'],
    'width_profile': ['capitulum_width_profile_cv'],
}
SEEDS = (20261101, 20261103, 20261105)
METHODS = ('photo_mean', 'observer_equal_mean')
EXPECTED_N, EXPECTED_USERS = 313350, 92839
SOURCE_SHA = '342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed'
ROSTER_SHA = 'a89a586d7aaff7ab4b6c936aa4bfc30c597c152a07bdb49ec7abddb7fee04fef'
MIN_OBSERVATIONS, MIN_OBSERVERS, MIN_TAXA = 5, 3, 5
BOOTSTRAPS = 1000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def side(user: str, seed: int) -> int:
    """Same global assignment as the previous observer-disjoint QC."""
    return hashlib.sha256(f'{seed}\0{user}'.encode('utf-8')).digest()[-1] & 1


def _unit(values: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(values, axis=-1, keepdims=True)
    return np.divide(values, norm, out=np.full_like(values, np.nan), where=norm > 1e-12)


def _divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan),
                     where=denominator > 1e-12)


def _metrics(a: np.ndarray, b: np.ndarray, circular: bool) -> dict:
    """Batch of paired-taxon panels, each shape [taxa, dimensions]."""
    ma, mb = a.mean(axis=1), b.mean(axis=1)
    ac, bc = a - ma[:, None, :], b - mb[:, None, :]
    va = np.mean(np.sum(ac * ac, axis=2), axis=1)
    vb = np.mean(np.sum(bc * bc, axis=2), axis=1)
    covariance = np.mean(np.sum(ac * bc, axis=2), axis=1)
    correlation = np.clip(_divide(covariance, np.sqrt(va * vb)), -1., 1.)
    concordance = np.clip(_divide(2 * covariance, va + vb + np.sum((ma - mb)**2, axis=1)), -1., 1.)
    if circular:
        dot = np.clip(np.sum(a * b, axis=2), -1., 1.)
        cross = b[:, :, 0] * a[:, :, 1] - b[:, :, 1] * a[:, :, 0]
        delta = np.arctan2(cross, dot)
        avg_sin, avg_cos = np.sin(delta).mean(axis=1), np.cos(delta).mean(axis=1)
        bias = np.rad2deg(np.arctan2(avg_sin, avg_cos))
        bias[np.hypot(avg_sin, avg_cos) <= 1e-12] = np.nan
        return {
            'centred_vector_alignment': correlation,
            'joint_vector_concordance': concordance,
            'mean_cosine_agreement': dot.mean(axis=1),
            'angular_bias_degrees': bias,
            'angular_mae_degrees': np.rad2deg(np.abs(delta)).mean(axis=1),
            'angular_rmse_degrees': np.sqrt(np.mean(np.rad2deg(delta)**2, axis=1)),
            'angular_absolute_difference_q95_degrees': np.quantile(np.rad2deg(np.abs(delta)), .95, axis=1),
        }
    x, y = a[:, :, 0], b[:, :, 0]
    rx, ry = rankdata(x, axis=1), rankdata(y, axis=1)
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    rho = np.clip(_divide(np.sum(rx * ry, axis=1),
                          np.sqrt(np.sum(rx * rx, axis=1) * np.sum(ry * ry, axis=1))), -1., 1.)
    diff = y - x
    return {
        'pearson_r': correlation, 'spearman_rho': rho, 'lin_ccc': concordance,
        'bias_b_minus_a': diff.mean(axis=1),
        'mae': np.abs(diff).mean(axis=1),
        'rmse': np.sqrt(np.mean(diff * diff, axis=1)),
        'absolute_difference_q95': np.quantile(np.abs(diff), .95, axis=1),
        'among_taxon_sd_a': np.sqrt(va), 'among_taxon_sd_b': np.sqrt(vb),
    }


def agreement(a: np.ndarray, b: np.ndarray, *, circular: bool = False,
              replicates: int = BOOTSTRAPS, seed: int = SEEDS[0]) -> dict:
    """Equal-taxon agreement, preserving A/B pairs during descriptive bootstrap."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.ndim == 1:
        a, b = a[:, None], b[:, None]
    if (a.shape != b.shape or a.ndim != 2 or len(a) < 2 or
            a.shape[1] != (2 if circular else 1) or
            not np.isfinite(a).all() or not np.isfinite(b).all() or
            not isinstance(replicates, int) or replicates < 1):
        raise ValueError('INVALID_AGREEMENT_INPUT')
    if circular:
        if (np.linalg.norm(a, axis=1) <= 1e-12).any() or (np.linalg.norm(b, axis=1) <= 1e-12).any():
            raise ValueError('UNDEFINED_HUE_DIRECTION')
        a, b = _unit(a), _unit(b)
    points = _metrics(a[None], b[None], circular)
    indices = np.random.default_rng(seed).integers(0, len(a), size=(replicates, len(a)))
    samples = _metrics(a[indices], b[indices], circular)
    intervals = {}
    for name, values in samples.items():
        finite = values[np.isfinite(values)]
        ok = len(finite) >= .95 * replicates
        # A circular bias interval needs circular unwrapping, not scalar quantiles.
        # Report its point estimate only rather than a misleading wraparound CI.
        if name == 'angular_bias_degrees':
            intervals[name] = {'status': 'circular_bias_interval_not_estimated',
                               'finite_replicates': int(len(finite)), 'lower': None, 'upper': None}
        else:
            intervals[name] = {
                'status': 'estimated_conditional_taxon_bootstrap' if ok else 'insufficient_finite_replicates',
                'finite_replicates': int(len(finite)),
                'lower': float(np.quantile(finite, .025)) if ok else None,
                'upper': float(np.quantile(finite, .975)) if ok else None,
            }
    return {
        'metrics': {name: float(v[0]) if np.isfinite(v[0]) else None for name, v in points.items()},
        'bootstrap_95': intervals, 'bootstrap_replicates': replicates,
        'interval_scope': 'paired taxa; conditional on fixed source observations and observer split; not observer, spatial or phylogenetic uncertainty',
    }


def _identity(frame: pd.DataFrame) -> None:
    for column, error in (('obs_id', 'OBSERVATION_IDENTITY'), ('user_id', 'OBSERVER_IDENTITY')):
        if column not in frame or frame[column].isna().any() or frame[column].astype(str).str.strip().eq('').any():
            raise ValueError(error)
    if frame.obs_id.astype(str).duplicated().any():
        raise ValueError('OBSERVATION_IDENTITY')


def _eligible(frame: pd.DataFrame, construct: str) -> pd.DataFrame:
    if construct not in FEATURES:
        raise ValueError('UNKNOWN_CONSTRUCT')
    columns = FEATURES[construct]
    required = {'obs_id', 'user_id', 'accepted_key', *columns, *(c + '__eligible' for c in columns)}
    if not required <= set(frame):
        raise ValueError('MISSING_REQUIRED_COLUMNS')
    _identity(frame)
    values = frame[columns].apply(pd.to_numeric, errors='coerce').to_numpy(float)
    valid = np.isfinite(values).all(axis=1) & frame.accepted_key.notna().to_numpy()
    valid &= frame.accepted_key.astype(str).str.strip().ne('').to_numpy()
    for column in columns:
        valid &= frame[column + '__eligible'].eq(True).to_numpy(bool)
    if construct == 'hue':
        valid &= np.linalg.norm(values, axis=1) >= .2
    good = frame.loc[valid, ['obs_id', 'user_id', 'accepted_key']].copy().astype(str)
    good.rename(columns={'accepted_key': 'taxon'}, inplace=True)
    labels = ['value'] if len(columns) == 1 else ['value0', 'value1']
    selected = _unit(values[valid]) if construct == 'hue' else values[valid]
    for index, label in enumerate(labels):
        good[label] = selected[:, index]
    return good.sort_values('obs_id', kind='mergesort').reset_index(drop=True)


def _pairs(good: pd.DataFrame, construct: str, seed: int) -> dict[str, pd.DataFrame]:
    empty = {name: pd.DataFrame() for name in METHODS}
    if good.empty:
        return empty
    d = good.copy()
    assignment = {u: side(u, seed) for u in sorted(d.user_id.unique())}
    d['half'] = d.user_id.map(assignment)
    if set(d.half.unique()) != {0, 1}:
        return empty
    if set(d.loc[d.half == 0, 'user_id']) & set(d.loc[d.half == 1, 'user_id']):
        raise ValueError('GLOBAL_OBSERVER_OVERLAP')
    values = ['value'] if construct != 'hue' else ['value0', 'value1']
    g = d.groupby(['taxon', 'half'], sort=True, observed=True)
    counts = g.agg(n_observations=('obs_id', 'size'), n_observers=('user_id', 'nunique'))
    means = g[values].mean()
    observer_means = d.groupby(['taxon', 'half', 'user_id'], sort=True, observed=True)[values].mean()
    balanced = observer_means.groupby(level=['taxon', 'half'], sort=True).mean()
    pairs = {}
    for method, table in zip(METHODS, (means, balanced)):
        summary = table.join(counts)
        p = summary.xs(0, level='half').add_suffix('_a').join(
            summary.xs(1, level='half').add_suffix('_b'), how='inner')
        keep = ((p.n_observations_a >= MIN_OBSERVATIONS) & (p.n_observations_b >= MIN_OBSERVATIONS) &
                (p.n_observers_a >= MIN_OBSERVERS) & (p.n_observers_b >= MIN_OBSERVERS))
        if construct == 'hue':
            for suffix in ('a', 'b'):
                keep &= np.linalg.norm(p[[v + '_' + suffix for v in values]], axis=1) >= .2
        p = p.loc[keep].copy()
        if construct == 'hue' and len(p):
            for suffix in ('a', 'b'):
                labels = [v + '_' + suffix for v in values]
                p.loc[:, labels] = _unit(p[labels].to_numpy(float))
        pairs[method] = p.sort_index()
    return pairs


def species_pairs(frame: pd.DataFrame, construct: str, seed: int) -> dict[str, pd.DataFrame]:
    """PRIVATE taxon means: never serialize these tables into a public artifact."""
    return _pairs(_eligible(frame, construct), construct, seed)


def _distribution(values) -> dict:
    v = np.asarray(values, float)
    if not len(v):
        return {'n': 0, 'minimum': None, 'median': None, 'maximum': None}
    return {'n': int(len(v)), 'minimum': float(v.min()), 'median': float(np.median(v)), 'maximum': float(v.max())}


def audit_construct(frame: pd.DataFrame, construct: str, *, seeds: tuple[int, ...] = SEEDS,
                    replicates: int = BOOTSTRAPS) -> dict:
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError('INVALID_OBSERVER_SEEDS')
    good = _eligible(frame, construct)
    prepared = {seed: _pairs(good, construct, seed) for seed in seeds}
    indices = [set(p.index) for methods in prepared.values() for p in methods.values()]
    common = sorted(set.intersection(*indices))
    result = {
        'construct': construct,
        'status': 'estimated_observer_disjoint_taxon_mean_QC' if len(common) >= MIN_TAXA else 'insufficient_supported_taxa',
        'n_trait_qc_observations': int(len(good)),
        'n_trait_qc_taxa_before_half_support': int(good.taxon.nunique()),
        'n_common_taxa_all_splits_and_methods': len(common),
        'eligible_taxon_counts_by_split': [
            {'seed': int(s), 'counts': {m: len(p) for m, p in prepared[s].items()}} for s in seeds],
        'same_taxa_all_splits_and_weighting_methods': True,
        'unit': 'degrees_image_frame_not_gravity' if construct == 'angle' else
                ('circular_hue_direction' if construct == 'hue' else 'original_photographic_trait_units'),
        'splits': [],
    }
    if len(common) < MIN_TAXA:
        return result
    labels = ['value'] if construct != 'hue' else ['value0', 'value1']
    for seed in seeds:
        for method in METHODS:
            p = prepared[seed][method].loc[common]
            a, b = (p[[v + '_' + suffix for v in labels]].to_numpy(float) for suffix in ('a', 'b'))
            result['splits'].append({
                'seed': int(seed), 'method': method, 'global_observer_overlap': 0,
                'support': {
                    'n_taxa': len(p),
                    'n_contributing_observations': int((p.n_observations_a + p.n_observations_b).sum()),
                    'observations_per_half_taxon': _distribution(np.r_[p.n_observations_a, p.n_observations_b]),
                    'observers_per_half_taxon': _distribution(np.r_[p.n_observers_a, p.n_observers_b]),
                },
                **agreement(a, b, circular=construct == 'hue', replicates=replicates, seed=int(seed)),
            })
    return result


def load_source(native: Path, receipt_path: Path, roster: Path, *, expected_n: int = EXPECTED_N,
                expected_users: int = EXPECTED_USERS, source_sha: str = SOURCE_SHA,
                roster_sha: str = ROSTER_SHA, roster_receipt: Path | None = None):
    receipt = json.loads(Path(receipt_path).read_text(encoding='utf-8'))
    digest = sha256(native)
    if (digest != source_sha or digest != receipt.get('prepared_table_sha256') or
            receipt.get('native_observations') != expected_n or
            receipt.get('environment_decoding_status') != 'physical_units_gsp_invalid_encoding_masked_v1'):
        raise ValueError('SOURCE_HASH_OR_RECEIPT_MISMATCH')
    actual_roster_sha = sha256(roster)
    roster_validation = 'fixed_historical_roster_bytes'
    upstream_receipt_sha = None
    if roster_receipt is None:
        if actual_roster_sha != roster_sha:
            raise ValueError('ROSTER_HASH_MISMATCH')
    else:
        # Existing roster writer uses gzip timestamps: regenerated bytes can
        # differ even with identical rows. Require an explicit attestation from
        # the existing authorized upstream run, never silently accept a new hash.
        upstream = json.loads(Path(roster_receipt).read_text(encoding='utf-8'))
        expected = {
            'schema': 'azami_authenticated_observer_roster_receipt_v1',
            'native_decoded_sha256': digest,
            'roster_table_sha256': actual_roster_sha,
            'source_encrypted_sha256': 'dfa98fa162a8029cb6477adf4c26ffc5819cc916a0513d32da14d78d845ab4b0',
            'source_plaintext_zip_sha256': '2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae',
            'n_observations': expected_n, 'n_observers': expected_users,
            'authentication_status': 'original_source_verified_upstream',
        }
        if any(upstream.get(k) != v for k, v in expected.items()):
            raise ValueError('ROSTER_RECEIPT_MISMATCH')
        roster_validation = 'upstream_attestation_checked_not_independently_reauthenticated'
        upstream_receipt_sha = sha256(roster_receipt)
    columns = ['obs_id', 'accepted_key']
    for values in FEATURES.values():
        for value in values:
            columns.extend([value, value + '__eligible'])
    d = pd.read_csv(native, usecols=columns, dtype={'obs_id': str}, low_memory=False)
    u = pd.read_csv(roster, usecols=['obs_id', 'user_id'], dtype=str, low_memory=False)
    if (len(d) != expected_n or len(u) != expected_n or not d.obs_id.is_unique or
            not u.obs_id.is_unique or d.obs_id.isna().any() or u.obs_id.isna().any() or
            u.user_id.isna().any() or u.user_id.nunique() != expected_users or
            set(d.obs_id) != set(u.obs_id)):
        raise ValueError('SOURCE_OR_ROSTER_IDENTITY_DRIFT')
    merged = d.merge(u, on='obs_id', how='left', sort=False, validate='one_to_one')
    _identity(merged)
    return merged, {'source_observations': expected_n, 'source_observers': expected_users,
                    'decoded_source_sha256': digest, 'observer_roster_sha256': actual_roster_sha,
                    'roster_source_validation': roster_validation,
                    'upstream_roster_receipt_sha256': upstream_receipt_sha,
                    'environment_and_coordinates_not_required_or_read': True}


def privacy_check(value):
    forbidden = {'obs_id', 'user_id', 'accepted_key', 'taxon', 'taxon_name', 'latitude', 'longitude',
                 'lat_cell', 'lon_cell', 'photo_id', 'predictions', 'model_weights', 'photographer_id',
                 'species_means', 'taxon_ids', 'observer_ids', 'private_rows'}
    if isinstance(value, dict):
        if forbidden & set(value):
            raise ValueError('PRIVATE_OUTPUT_FIELD')
        for child in value.values():
            privacy_check(child)
    elif isinstance(value, list):
        for child in value:
            privacy_check(child)


def run(native: Path, receipt: Path, roster: Path, output: Path, *,
        roster_receipt: Path | None = None) -> dict:
    frame, source = load_source(native, receipt, roster, roster_receipt=roster_receipt)
    results = [audit_construct(frame, construct) for construct in FEATURES]
    out = {
        'schema': 'azami_observer_disjoint_species_mean_QC_v1',
        'status': 'source_receipts_checked_taxon_mean_QC_not_RF', 'source': source,
        'design': {'constructs': list(FEATURES), 'methods': list(METHODS), 'seeds': list(SEEDS),
                   'min_observations_per_half_taxon': MIN_OBSERVATIONS,
                   'min_observers_per_half_taxon': MIN_OBSERVERS, 'min_common_taxa': MIN_TAXA,
                   'bootstrap_replicates': BOOTSTRAPS,
                   'taxon_panel': 'intersection across all seeds AND both weighting methods per construct',
                   'no_environmental_predictors': True, 'no_geographic_cell_requirement': True},
        'scope': 'repeatability of photographic taxon means, not within-taxon geography, calibration, cause, or RF skill',
        'limits': [
            'Not necessarily the exact S1 environmental-complete cohort; S2 is not included.',
            'Global observer separation does not separate biological individuals or shared regional/seasonal photo conditions.',
            'Paired-taxon intervals condition on fixed observations and observer groups; no observer, spatial or phylogenetic resampling.',
            'Seeds reuse observers and are not independent experiments.',
            'High species-mean agreement does not validate gravity angles, pigment levels, or environmental effects.',
            'No RF-failure mechanism or noise ceiling is inferred from these correlations.',
        ],
        'software': {'python': platform.python_version(), 'numpy': np.__version__,
                     'pandas': pd.__version__, 'scipy': scipy.__version__},
        'private_individual_or_taxon_data_not_exported': True,
        'original_GEB_S1_S2_H3c_and_RF_unchanged': True, 'results': results,
    }
    privacy_check(out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'receipt', 'roster', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--roster-receipt', type=Path,
                        help='Existing authorized upstream attestation when roster gzip bytes were regenerated')
    args = parser.parse_args()
    result = run(args.native, args.receipt, args.roster, args.output, roster_receipt=args.roster_receipt)
    print(json.dumps({'status': result['status'], 'n_constructs': len(result['results'])}), flush=True)


if __name__ == '__main__':
    main()

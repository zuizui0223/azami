"""Original-S1 within-taxon geographical split-half reliability of photo traits.

Exploratory quality-control diagnostic, not an environmental RF fit, ICC,
measurement-error decomposition, adaptation test, or replacement for frozen GEB.
All observation IDs, taxon identities, photographs and coordinates remain private.
Only pooled taxon-equal statistics are emitted outside protected execution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analysis.v4.native_integration import FEATURES

EXPECTED_N = 313350
SEEDS = (20261101, 20261103, 20261105)
GRAINS = (1, 5)
KEYS = ['taxon', 'lat_cell', 'lon_cell']
MIN_HALF = 3
MIN_CELLS = 8
MIN_PHOTOS = 40
N_BOOT = 500
N_NULL = 199


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _distribution(values: np.ndarray) -> dict:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return {'n': 0, 'q05': None, 'q50': None, 'q95': None}
    q = np.quantile(x, [.05, .5, .95])
    return {'n': int(len(x)), 'q05': float(q[0]),
            'q50': float(q[1]), 'q95': float(q[2])}


def _unit(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    d = np.linalg.norm(x, axis=1, keepdims=True)
    return np.divide(x, d, out=np.full_like(x, np.nan), where=d > 1e-12)


def _hash_obs_ids(frame: pd.DataFrame, seeds: tuple[int, ...]) -> pd.DataFrame:
    """Hash original obs_id, not location or phenotype; never export hashes."""
    if not frame.obs_id.is_unique or frame.obs_id.isna().any():
        raise ValueError('DUPLICATE_OR_MISSING_ORIGINAL_OBS_ID')
    ids = frame.obs_id.astype(str).to_numpy()
    for seed in seeds:
        frame[f'_split_hash_{seed}'] = np.fromiter(
            (int.from_bytes(hashlib.sha256(f'{seed}\0{pid}'.encode('utf-8')).digest()[:8], 'big')
             for pid in ids), dtype=np.uint64, count=len(ids))
    return frame


def _pair_cells(good: pd.DataFrame, cols: list[str], seed: int, *,
                min_half: int, min_cells: int, min_photos: int,
                hue: bool) -> pd.DataFrame:
    hkey = f'_split_hash_{seed}'
    # Within-cell hash ordering, alternating halves, guarantees balance when n is odd.
    ranked = good.sort_values(KEYS + [hkey, 'obs_id'], kind='mergesort').copy()
    ranked['half'] = (ranked.groupby(KEYS, sort=False, observed=True)
                      .cumcount().to_numpy() % 2).astype(np.int8)
    grouped = ranked.groupby(KEYS + ['half'], sort=False, observed=True)
    med = grouped[cols].median()
    med['n'] = grouped.size()
    a = med.xs(0, level='half').add_suffix('_a')
    b = med.xs(1, level='half').add_suffix('_b')
    paired = a.join(b, how='inner').reset_index()
    paired = paired.loc[(paired.n_a >= min_half) &
                        (paired.n_b >= min_half)].copy()
    if hue and len(paired):
        norm_a = np.linalg.norm(paired[[c+'_a' for c in cols]].to_numpy(float), axis=1)
        norm_b = np.linalg.norm(paired[[c+'_b' for c in cols]].to_numpy(float), axis=1)
        paired = paired.loc[(norm_a >= .20) & (norm_b >= .20)].copy()
    if paired.empty:
        return paired
    tax = paired.groupby('taxon', observed=True).agg(
        n_cells=('n_a', 'size'), n_a=('n_a', 'sum'), n_b=('n_b', 'sum'))
    keep = tax.index[(tax.n_cells >= min_cells) &
                     ((tax.n_a + tax.n_b) >= min_photos)]
    return paired.loc[paired.taxon.isin(keep)].copy()


def _evaluate_pairs(paired: pd.DataFrame, cols: list[str], *,
                    hue: bool, seed: int) -> dict:
    if paired.empty:
        return {'status': 'insufficient_independently_supported_cells',
                'n_taxa': 0, 'n_cells': 0}
    r_values, spearmans, diff_rmses, spread_rmses, ratios = [], [], [], [], []
    arrays = []
    for _, block in paired.groupby('taxon', observed=True, sort=False):
        a = block[[c+'_a' for c in cols]].to_numpy(float)
        b = block[[c+'_b' for c in cols]].to_numpy(float)
        if hue:
            a, b = _unit(a), _unit(b)
        a0 = a - a.mean(axis=0)
        b0 = b - b.mean(axis=0)
        denom = np.sqrt(np.sum(a0*a0) * np.sum(b0*b0))
        if not np.isfinite(denom) or denom <= 1e-12:
            continue
        r = float(np.sum(a0*b0) / denom)
        diff = float(np.sqrt(np.mean(np.sum((a-b)**2, axis=1))))
        spread = float(np.sqrt((np.mean(np.sum(a0*a0, axis=1)) +
                                np.mean(np.sum(b0*b0, axis=1))) / 2))
        r_values.append(np.clip(r, -1., 1.))
        diff_rmses.append(diff)
        spread_rmses.append(spread)
        ratios.append(diff/spread if spread > 1e-12 else np.nan)
        if not hue:
            spearmans.append(float(spearmanr(a[:, 0], b[:, 0]).statistic))
        arrays.append((a0, b0, denom))
    r_values = np.asarray(r_values, dtype=float)
    if len(r_values) < 5:
        return {'status': 'insufficient_finite_taxon_correlations',
                'n_taxa': int(len(r_values)), 'n_cells': int(len(paired))}
    # Species-level resampling is descriptive, not phylogenetically independent.
    rng = np.random.default_rng(8191 + int(seed))
    boot = np.median(r_values[rng.integers(0, len(r_values),
                        size=(N_BOOT, len(r_values)))], axis=1)
    null = []
    for _ in range(N_NULL):
        fake = [np.sum(a0 * b0[rng.permutation(len(b0))]) / denom
                for a0, b0, denom in arrays]
        null.append(float(np.median(fake)))
    return {
        'status': 'estimated_exploratory_reliability',
        'n_taxa': int(len(r_values)),
        'n_qualified_cells': int(len(paired)),
        'n_contributing_photos': int((paired.n_a+paired.n_b).sum()),
        'photos_per_half_cell': _distribution(np.r_[paired.n_a.to_numpy(), paired.n_b.to_numpy()]),
        'cells_per_taxon': _distribution(paired.groupby('taxon', observed=True).size().to_numpy()),
        'taxon_equal_spatial_alignment': _distribution(r_values),
        'median_alignment_taxon_bootstrap_95': [float(np.quantile(boot, .025)),
                                                 float(np.quantile(boot, .975))],
        'fraction_taxa_alignment_positive': float(np.mean(r_values > 0)),
        'fraction_taxa_alignment_at_least_0_5': float(np.mean(r_values >= .5)),
        'taxon_equal_spearman': (None if hue else _distribution(spearmans)),
        'taxon_equal_half_median_difference_RMSE': _distribution(diff_rmses),
        'taxon_equal_among_cell_half_median_RMSE': _distribution(spread_rmses),
        'half_difference_over_among_cell_RMSE': _distribution(ratios),
        'within_taxon_cell_labels_permuted_null_median_alignment': _distribution(null),
        'null_permutations': N_NULL,
        'taxon_bootstrap_replicates': N_BOOT,
        'not_an_ICC_measurement_error_or_causal_effect': True,
    }


def audit_construct(frame: pd.DataFrame, construct: str, *, degrees: int,
                    seeds: tuple[int, ...] = SEEDS,
                    min_half: int = MIN_HALF,
                    min_cells: int = MIN_CELLS,
                    min_photos: int = MIN_PHOTOS) -> dict:
    if construct not in FEATURES or degrees not in GRAINS or min_half < 2:
        raise ValueError('INVALID_REGISTERED_SPLIT_HALF_DESIGN')
    cols = list(FEATURES[construct]); hue = construct == 'hue'
    necessary = {'accepted_key', 'latitude', 'longitude', 'obs_id', *cols,
                 *(c+'__eligible' for c in cols),
                 *(f'_split_hash_{s}' for s in seeds)}
    if not necessary.issubset(frame.columns):
        raise ValueError('ORIGINAL_PHOTO_SPLIT_SOURCE_COLUMNS_MISSING')
    lat = pd.to_numeric(frame.latitude, errors='coerce').to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors='coerce').to_numpy(float)
    y = frame[cols].apply(pd.to_numeric, errors='coerce').to_numpy(float)
    good_mask = (np.isfinite(lat) & np.isfinite(lon) &
                 (np.abs(lat) <= 90) & (np.abs(lon) <= 180) &
                 np.isfinite(y).all(axis=1) & frame.accepted_key.notna().to_numpy())
    for c in cols:
        good_mask &= frame[c+'__eligible'].eq(True).to_numpy()
    if hue:
        good_mask &= (np.linalg.norm(y, axis=1) >= .20)
    subset = frame.loc[good_mask, ['obs_id', 'accepted_key', *(f'_split_hash_{s}' for s in seeds)]].copy()
    subset.rename(columns={'accepted_key': 'taxon'}, inplace=True)
    subset['taxon'] = subset.taxon.astype(str)
    subset['lat_cell'] = np.floor((lat[good_mask]+90)/degrees).astype(int)
    subset['lon_cell'] = np.floor((lon[good_mask]+180)/degrees).astype(int)
    if hue:
        y_good = _unit(y[good_mask])
    else:
        y_good = y[good_mask]
    for i, c in enumerate(cols):
        subset[c] = y_good[:, i]
    splits = []
    for seed in seeds:
        paired = _pair_cells(subset, cols, seed, min_half=min_half,
                             min_cells=min_cells, min_photos=min_photos, hue=hue)
        outcome = _evaluate_pairs(paired, cols, hue=hue, seed=seed)
        splits.append({'split_seed': int(seed), **outcome})
    return {'construct': construct, 'degrees': degrees,
            'n_original_photo_trait_QC_eligible': int(good_mask.sum()),
            'splits': splits}


def run(native: Path, receipt_path: Path, output: Path, *,
        expected_n: int = EXPECTED_N,
        constructs: tuple[str, ...] | None = None,
        grains: tuple[int, ...] = GRAINS,
        seeds: tuple[int, ...] = SEEDS) -> dict:
    receipt = json.loads(Path(receipt_path).read_text(encoding='utf-8'))
    input_sha = sha256(native)
    if (input_sha != receipt.get('prepared_table_sha256') or
            receipt.get('native_observations') != expected_n or
            receipt.get('environment_decoding_status') !=
            'physical_units_gsp_invalid_encoding_masked_v1'):
        raise ValueError('ORIGINAL_S1_DECODE_AND_SHA_NOT_AUTHENTICATED')
    chosen = tuple(FEATURES) if constructs is None else tuple(constructs)
    if (not chosen or len(chosen) != len(set(chosen)) or
            not set(chosen).issubset(FEATURES) or
            not grains or len(grains) != len(set(grains)) or
            not set(grains).issubset(GRAINS) or
            not seeds or len(seeds) != len(set(seeds))):
        raise ValueError('INVALID_SPLIT_HALF_SELECTION')
    required = ['obs_id', 'accepted_key', 'latitude', 'longitude']
    for c in chosen:
        for trait in FEATURES[c]:
            required.extend((trait, trait+'__eligible'))
    df = pd.read_csv(native, usecols=list(dict.fromkeys(required)),
                     dtype={'obs_id': str}, low_memory=False)
    if len(df) != expected_n:
        raise ValueError('ORIGINAL_S1_ROW_COUNT_DRIFT')
    df = _hash_obs_ids(df, seeds)
    results = [audit_construct(df, c, degrees=d, seeds=seeds)
               for c in chosen for d in grains]
    out = {
        'schema': 'azami_v4_original_s1_independent_photo_split_half_geo_reliability_v1',
        'status': 'source_authenticated_exploratory_QC_not_RF_or_GEB_reanalysis',
        'source_original_S1_observations': expected_n,
        'source_decoded_sha256': input_sha,
        'design': {'constructs': list(chosen), 'grains_degrees': list(grains),
                   'fixed_hash_seeds': list(seeds),
                   'minimum_photos_per_half_cell': MIN_HALF,
                   'minimum_geographical_cells_per_taxon': MIN_CELLS,
                   'minimum_contributing_photos_per_taxon': MIN_PHOTOS,
                   'unit_of_partition': 'unique_original_observation_id',
                   'split_algorithm': 'SHA256(seed+obs_id), balanced alternating within taxon-cell'},
        'spatial_alignment_definition': 'taxon-wise centred cross-half cell-vector cosine; scalar Pearson analogue, circular hue joint unit-sin/cos',
        'null_definition': 'within taxon permute half-B cell labels; exploratory calibration only',
        'interpretation_limit': 'reproducible photo-derived geographic pattern, not genetic or adaptive signal, independent individuals, direct camera measurement accuracy or RF transfer skill',
        'photo_QC_not_RF28_or_Soil34_environment_complete_case_cohort': True,
        'frozen_GEB_S1_S2_H3c_untouched': True,
        'source_taxon_cell_photo_location_data_not_exported': True,
        'results': results,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2, ensure_ascii=False, allow_nan=False)+'\n',
                      encoding='utf-8')
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native', type=Path, required=True)
    p.add_argument('--receipt', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    r = run(a.native, a.receipt, a.output)
    print(json.dumps({'status': r['status'],
                      'n_construct_grain_summaries': len(r['results']),
                      'n_partition_replicates': len(SEEDS)}), flush=True)


if __name__ == '__main__':
    main()

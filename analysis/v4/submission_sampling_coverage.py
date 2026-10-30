"""Coarse coverage of the exact V4 analytical cohorts, never source points.

This module reuses the frozen cohort selectors and verifies their coefficients.
It makes no new ecological test. Public rows contain pooled five-degree cell
indices and counts only; cells with fewer than five records are omitted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

CELL_DEGREES = 5
MIN_CELL_RECORDS = 5
EXPECTED_PANEL_COUNTS = {
    's1_angle': 168907,
    'arvense_native': 22936,
    'arvense_introduced': 26171,
    'vulgare_native': 17335,
    'vulgare_introduced': 27200,
}
PANEL_FIELDS = {
    'panel', 'total_records', 'total_taxa', 'total_cells', 'shown_records',
    'suppressed_records', 'shown_cells', 'suppressed_cells', 'cells',
}
CELL_FIELDS = {'lat_bin', 'lon_bin', 'records', 'taxa'}


def aggregate_panel(frame: pd.DataFrame, name: str, *, taxon_column: str):
    coordinates = frame[['latitude', 'longitude']].to_numpy(float)
    if (not np.isfinite(coordinates).all() or
            np.any(np.abs(coordinates[:, 0]) > 90) or
            np.any(np.abs(coordinates[:, 1]) > 180)):
        raise ValueError('INVALID_SOURCE_COORDINATE')
    if frame[taxon_column].isna().any():
        raise ValueError('MISSING_SOURCE_TAXON')
    work = pd.DataFrame({
        'lat_bin': np.minimum(np.floor((coordinates[:, 0] + 90) / CELL_DEGREES), 35).astype(int),
        'lon_bin': np.minimum(np.floor((coordinates[:, 1] + 180) / CELL_DEGREES), 71).astype(int),
        'taxon': frame[taxon_column].astype(str).to_numpy(),
    })
    counts = work.groupby(['lat_bin', 'lon_bin'], sort=True).agg(
        records=('taxon', 'size'), taxa=('taxon', 'nunique')).reset_index()
    shown = counts.loc[counts.records >= MIN_CELL_RECORDS]
    result = {
        'panel': name,
        'total_records': int(len(frame)),
        'total_taxa': int(frame[taxon_column].nunique()),
        'total_cells': int(len(counts)),
        'shown_records': int(shown.records.sum()),
        'suppressed_records': int(counts.loc[counts.records < MIN_CELL_RECORDS, 'records'].sum()),
        'shown_cells': int(len(shown)),
        'suppressed_cells': int(len(counts) - len(shown)),
        'cells': [{key: int(row[key]) for key in ('lat_bin', 'lon_bin', 'records', 'taxa')}
                  for row in shown.to_dict('records')],
    }
    validate_panel(result)
    return result


def validate_panel(panel):
    if set(panel) != PANEL_FIELDS or not isinstance(panel['panel'], str):
        raise ValueError('UNEXPECTED_PUBLIC_PANEL_FIELD')
    for key in PANEL_FIELDS - {'panel', 'cells'}:
        if type(panel[key]) is not int or panel[key] < 0:
            raise ValueError('INVALID_PUBLIC_PANEL_COUNT')
    cells = panel['cells']
    if not isinstance(cells, list):
        raise ValueError('INVALID_PUBLIC_CELL_COLLECTION')
    seen = set()
    for row in cells:
        if set(row) != CELL_FIELDS or any(type(row[k]) is not int for k in CELL_FIELDS):
            raise ValueError('UNEXPECTED_PUBLIC_CELL_FIELD')
        if (not 0 <= row['lat_bin'] <= 35 or not 0 <= row['lon_bin'] <= 71 or
                row['records'] < MIN_CELL_RECORDS or
                not 1 <= row['taxa'] <= min(row['records'], panel['total_taxa'])):
            raise ValueError('INVALID_OR_UNSUPPRESSED_PUBLIC_CELL')
        key = (row['lat_bin'], row['lon_bin'])
        if key in seen:
            raise ValueError('DUPLICATE_PUBLIC_CELL')
        seen.add(key)
    if (panel['shown_records'] != sum(row['records'] for row in cells) or
            panel['shown_cells'] != len(cells) or
            panel['total_cells'] != panel['shown_cells'] + panel['suppressed_cells'] or
            panel['total_records'] != panel['shown_records'] + panel['suppressed_records']):
        raise ValueError('PUBLIC_COHORT_COUNT_RECONCILIATION_FAILED')
    if not panel['suppressed_cells'] <= panel['suppressed_records'] <= 4 * panel['suppressed_cells']:
        raise ValueError('INVALID_SUPPRESSION_ACCOUNTING')


def validate_coverage(data):
    expected = {'schema', 'cell_degrees', 'minimum_records_per_displayed_cell',
                'source_verification', 'panel_overlap', 'panels'}
    if (set(data) != expected or data['schema'] != 'V4_SAMPLING_COVERAGE_20261030' or
            data['cell_degrees'] != CELL_DEGREES or
            data['minimum_records_per_displayed_cell'] != MIN_CELL_RECORDS):
        raise ValueError('UNRECOGNISED_COVERAGE_SCHEMA')
    names = set()
    for panel in data['panels']:
        validate_panel(panel)
        name = panel['panel']
        if name in names:
            raise ValueError('DUPLICATE_COVERAGE_PANEL')
        names.add(name)
        if name in EXPECTED_PANEL_COUNTS and panel['total_records'] != EXPECTED_PANEL_COUNTS[name]:
            raise ValueError('FROZEN_COVERAGE_DENOMINATOR_CHANGED')
    if names != set(EXPECTED_PANEL_COUNTS) | {'h3c_primary'}:
        raise ValueError('INCOMPLETE_COVERAGE_PANELS')
    if next(p for p in data['panels'] if p['panel'] == 's1_angle')['total_taxa'] != 178:
        raise ValueError('S1_COVERAGE_TAXA_CHANGED')
    if next(p for p in data['panels'] if p['panel'] == 'h3c_primary')['total_taxa'] != 63:
        raise ValueError('H3C_COVERAGE_TAXA_CHANGED')
    def keys(node):
        if isinstance(node, dict):
            for key, value in node.items():
                yield key
                yield from keys(value)
        elif isinstance(node, list):
            for value in node:
                yield from keys(value)
    forbidden = {'obs_id', 'user_id', 'photo_id', 'photographer_id', 'accepted_key',
                 'latitude', 'longitude', 'taxon_name'}
    if forbidden & set(keys(data)):
        raise ValueError('PROTECTED_SOURCE_FIELD_IN_COVERAGE_OUTPUT')


def build_coverage(native_path, receipt_path, focal_path, extra_path, measured_path):
    from analysis.v4.s1_climate_trait_hypervolume import load_source, geo_aggregate
    from analysis.v4.s1_bio18_geographic_replication import ANGLE, reproduce_frozen_baseline
    from analysis.v4.s1_s2_extreme_lego import read_tables
    from analysis.v4.s1_s2_vif5_si10_alltraits import selected_source
    from analysis.v4.s2_original_colour_spatial_block_sensitivity import (
        EXPECTED, original_sample, _original_fit,
    )

    # This also verifies the decoded S1 table's SHA-256 against its receipt.
    h3, provenance = load_source(native_path, receipt_path)
    if len(h3) != 157822 or h3.accepted_key.nunique() != 242:
        raise ValueError('FROZEN_H3_MATCHED_SOURCE_CHANGED')
    native, focal = read_tables(native_path, focal_path, extra_path, measured_path)
    angle, status = selected_source(native, [ANGLE], min_taxa=20)
    if angle is None or len(angle) != 168907 or angle.accepted_key.nunique() != 178:
        raise ValueError('FROZEN_S1_ANGLE_COHORT_CHANGED')
    baseline = reproduce_frozen_baseline(angle)
    panels = [aggregate_panel(angle, 's1_angle', taxon_column='accepted_key')]

    primary, primary_support = geo_aggregate(h3, 1)
    secondary, secondary_support = geo_aggregate(h3, 5)
    if (primary.accepted_key.nunique(), len(primary),
            secondary.accepted_key.nunique(), len(secondary)) != (63, 7117, 31, 814):
        raise ValueError('FROZEN_H3_GEOGRAPHIC_SUPPORT_CHANGED')
    primary_records = h3.loc[h3.accepted_key.isin(primary.accepted_key)]
    panels.append(aggregate_panel(primary_records, 'h3c_primary', taxon_column='accepted_key'))

    s2_verified = []
    s2_seen = {}
    for species, state, trait, exposure, expected_n, expected_beta in EXPECTED:
        cohort, _ = original_sample(focal, species, state, trait)
        if len(cohort) != expected_n:
            raise ValueError('FROZEN_S2_COLOUR_DENOMINATOR_CHANGED')
        _, _, _, _, fit = _original_fit(cohort, trait, exposure)
        if not np.isclose(fit['beta'], expected_beta, rtol=0, atol=2e-7):
            raise ValueError('FROZEN_S2_COLOUR_COEFFICIENT_CHANGED')
        panel_id = species.split()[-1] + '_' + state
        identities = frozenset(cohort.obs_id.astype(str))
        if panel_id in s2_seen:
            if s2_seen[panel_id] != identities:
                raise ValueError('CHROMA_AND_LIGHTNESS_COVERAGE_NOT_IDENTICAL')
        else:
            s2_seen[panel_id] = identities
            panels.append(aggregate_panel(cohort, panel_id, taxon_column='taxon_name'))
        s2_verified.append({'species': species, 'range': state, 'trait': trait,
                            'exposure': exposure, 'n': expected_n,
                            'expected_beta': expected_beta, 'recovered_beta': fit['beta']})
    result = {
        'schema': 'V4_SAMPLING_COVERAGE_20261030',
        'cell_degrees': CELL_DEGREES,
        'minimum_records_per_displayed_cell': MIN_CELL_RECORDS,
        'source_verification': {
            'native_source_records': 313350,
            'focal_source_records': 297622,
            'decoded_native_sha256': provenance['decoded_sha256'],
            's1_original_baseline': baseline,
            's2_original_six_rows': s2_verified,
            'h3_joint_eligible_records': 157822,
            'h3_joint_eligible_taxa': 242,
            'h3_primary_taxa': 63,
            'h3_primary_taxon_cells': 7117,
            'h3_secondary_taxa': 31,
            'h3_secondary_taxon_cells': 814,
            'original_inference_and_correction_families_unchanged': True,
        },
        'panel_overlap': 'Analyses reuse observations; panel totals are not independent and must not be summed.',
        'panels': panels,
    }
    validate_coverage(result)
    return result


def main():
    parser = argparse.ArgumentParser()
    for name in ('native', 'receipt', 'focal', 'extra', 'measured', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    result = build_coverage(args.native, args.receipt, args.focal, args.extra, args.measured)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'EXACT_COHORT_COVERAGE_VERIFIED',
                      'panels': [{k: p[k] for k in ('panel', 'total_records', 'total_taxa',
                                                    'shown_records', 'suppressed_records')}
                                 for p in result['panels']]}), flush=True)


if __name__ == '__main__':
    main()

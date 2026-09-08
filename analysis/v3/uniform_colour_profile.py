"""Characterize a versioned colour candidate from saved numerical probes only."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

from .protected_artifacts import new_json, require
from .summarize_perturbations import strata, summarize_frame
from .workflow import ROOT, canonical_digest, digest, text_digest

RECIPE = ROOT/'analysis/v3/uniform_floral_chroma_recipe.json'
RECEIPT = ROOT/'reproducibility/v3_cached_perturbation_20260907.json'
METRICS = ('uniform_floral_chroma', 'matched_floral_chroma', 'matched_green_chroma',
           'matched_floral_minus_green_chroma', 'uniform_minus_legacy_chroma')


def chroma(value):
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value) and 0 <= value <= math.sqrt(2)*128


def extract(payload):
    """No raw value is silently promoted by missing status or context support."""
    result = payload.get('result', {})
    diagnostics = result.get('diagnostics', {})
    paired = diagnostics.get('paired_colour', {})
    flower = paired.get('floral_union', {})
    green = paired.get('green_non_head_context', {})
    rows = [r for r in result.get('endpoints', []) if r['endpoint_id'] == 'corolla_lab_chroma']
    require(len(rows) <= 1, 'Duplicate registered chroma endpoint')
    old = rows[0] if rows else {}
    floral_value, green_value = flower.get('lab_chroma'), green.get('lab_chroma')
    old_value = old.get('value')
    floral_ok = (payload.get('status') == 'measured' and old.get('status') == 'usable'
                 and chroma(floral_value) and flower.get('n_pixels', 0) > 0
                 and diagnostics.get('head_min_dimension_px', 0) >= 96)
    matched = floral_ok and green.get('support_status') == 'available' and chroma(green_value) and green.get('n_pixels', 0) > 0
    values = {
        'uniform_floral_chroma': (floral_value, floral_ok),
        'matched_floral_chroma': (floral_value, matched),
        'matched_green_chroma': (green_value, matched),
        'matched_floral_minus_green_chroma': (floral_value - green_value if chroma(floral_value) and chroma(green_value) else None, matched),
        'uniform_minus_legacy_chroma': (floral_value - old_value if chroma(floral_value) and chroma(old_value) else None, floral_ok and chroma(old_value)),
    }
    return {'condition_status': payload.get('status', 'missing'),
            'legacy_status': old.get('status', 'missing'),
            'floral_n_pixels': flower.get('n_pixels'), 'green_n_pixels': green.get('n_pixels'),
            'green_support_status': green.get('support_status', 'missing'),
            **{f'{key}_{field}': value[i] for key, value in values.items() for i, field in enumerate(('value', 'usable'))}}


def summarize(table, conditions):
    """All heads remain in each condition and stratum, even if no pair survives."""
    base = table[table.condition.eq('baseline')]
    require(not table.duplicated(['head_id', 'condition']).any(), 'Duplicate head-condition')
    require(len(table) == len(base)*len(conditions), 'Incomplete condition denominator')
    records = []
    identifiers = ['head_id', 'sha256', 'component_id', 'development_exposed']
    for condition in conditions:
        changed = table[table.condition.eq(condition)]
        require(set(changed.head_id) == set(base.head_id), 'Condition head membership differs')
        for metric in METRICS:
            b = base[identifiers + [metric+'_value', metric+'_usable']].rename(columns={metric+'_value':'baseline_value', metric+'_usable':'baseline_usable'})
            p = changed[['head_id', metric+'_value', metric+'_usable']].rename(columns={metric+'_value':'perturbed_value', metric+'_usable':'perturbed_usable'})
            frame = b.merge(p, on='head_id', validate='one_to_one')
            frame['baseline_value'] = frame.baseline_value.astype(float)
            frame['perturbed_value'] = frame.perturbed_value.astype(float)
            require(not (frame.baseline_usable & ~np.isfinite(frame.baseline_value)).any()
                    and not (frame.perturbed_usable & ~np.isfinite(frame.perturbed_value)).any(), 'Nonfinite eligible candidate')
            frame['signed_change'] = frame.perturbed_value - frame.baseline_value
            frame['absolute_change'] = frame.signed_change.abs()
            for label, subset in strata(frame):
                records.append({'condition':condition, 'metric_id':metric, 'exposure_stratum':label, **summarize_frame(subset)})
    return records


def run(perturbation, dependence, out):
    recipe = json.loads(RECIPE.read_text(encoding='utf-8'))
    source = json.loads(RECEIPT.read_text(encoding='utf-8'))
    require(tuple(recipe['metrics']) == METRICS and not recipe['ecological_fitting_authorized'] and not recipe['replaces_registered_chroma'], 'Candidate scope changed')
    require(source['status'] == 'TECHNICAL_PERTURBATION_COMPLETED', 'Pinned probe grid is incomplete')
    perturbation, dependence, out = [p.resolve() for p in (perturbation, dependence, out)]
    require(out.is_relative_to(ROOT/'local_data') and not out.exists(), 'Use a new ignored private output')
    pins = {'perturbation': source['output_database_sha256'], 'dependence':source['execution_contract']['input_database_sha256']['dependence_groups']}
    for name, path in [('perturbation', perturbation), ('dependence', dependence)]:
        require(digest(path) == pins[name], 'Pinned numerical input differs')
    for name, sha in source['execution_contract']['feature_specification']['source_sha256_text_lf'].items():
        require(text_digest(ROOT/name) == sha, 'Saved feature implementation differs')
    conditions = [r['id'] for r in source['execution_contract']['contract']['conditions']]
    out.mkdir(parents=True)
    execution = {'recipe':recipe, 'recipe_canonical_sha256':canonical_digest(recipe),
                 'input_database_sha256':pins, 'source_receipt_canonical_sha256':canonical_digest(source),
                 'implementation_sha256_text_lf':text_digest(Path(__file__)),
                 'summary_helper_sha256_text_lf':text_digest(ROOT/'analysis/v3/summarize_perturbations.py'),
                 'trait_environment_join_performed':False}
    new_json(out/'execution_contract.json', execution)
    probes = sqlite3.connect(perturbation.as_uri()+'?mode=ro', uri=True)
    groups = sqlite3.connect(dependence.as_uri()+'?mode=ro', uri=True)
    try:
        heads = pd.read_sql_query('SELECT head_id,sha256,component_id,status FROM jobs', probes)
        require(len(heads) == source['counts']['scheduled_heads'] and heads.status.eq('completed').all(), 'Head completion differs')
        exposure = pd.read_sql_query('SELECT component_id,development_pool_exposed AS development_exposed FROM components WHERE n_cached_objects>0', groups)
        heads = heads.merge(exposure, on='component_id', how='left', validate='many_to_one')
        require(heads.development_exposed.isin([0, 1]).all(), 'Missing development annotation')
        rows = []
        for head, condition, status, raw in probes.execute('SELECT head_id,condition,status,payload_json FROM results ORDER BY head_id,condition'):
            payload = json.loads(raw)
            require(payload['condition'] == condition and payload['status'] == status and condition in conditions, 'Condition identity differs')
            rows.append({'head_id':head, 'condition':condition, **extract(payload)})
        table = pd.DataFrame(rows).merge(heads.drop(columns=['status']), on='head_id', how='left', validate='many_to_one')
        require(len(table) == source['counts']['condition_rows'] and table.component_id.notna().all(), 'Saved condition membership differs')
        records = summarize(table, conditions)
    finally:
        probes.close()
        groups.close()
    require(len(records) == len(conditions)*len(METRICS)*3, 'Summary grid incomplete')
    table.to_csv(out/'candidate_colour_private.csv', index=False, lineterminator='\n')
    pd.DataFrame(records).to_csv(out/'technical_summary.csv', index=False, lineterminator='\n')
    public = {'status':'VERSIONED_UNIFORM_COLOUR_TECHNICAL_PROFILE_EXECUTED_NO_ECOLOGICAL_ADMISSION',
              'execution_contract':execution, 'scheduled_heads':len(heads), 'condition_rows':len(table),
              'summary_rows':len(records), 'candidate_table_sha256':digest(out/'candidate_colour_private.csv'),
              'technical_summary_sha256':digest(out/'technical_summary.csv'), 'records':records,
              'image_requests':0, 'environment_values_read':0, 'ecological_models_executed':0,
              'independent_accuracy_estimated':False, 'ecological_fitting_authorized':False,
              'limits':recipe['limits']}
    new_json(out/'public_report.json', public)
    return public


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('perturbation', 'dependence', 'out'):
        parser.add_argument('--'+arg, type=Path, required=True)
    args = parser.parse_args()
    report = run(args.perturbation, args.dependence, args.out)
    print(json.dumps({k:v for k,v in report.items() if k not in ('records', 'execution_contract')}))


if __name__ == '__main__':
    main()

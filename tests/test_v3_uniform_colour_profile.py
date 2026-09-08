import copy
import json

import numpy as np
import pandas as pd
import pytest

from analysis.v3 import uniform_colour_profile as profile
from analysis.v3.image_features import colour_summary


def payload(flower=40., green=20., old=60., status='usable'):
    return {'status':'measured', 'result':{
        'endpoints':[{'endpoint_id':'corolla_lab_chroma','value':old,'status':status}],
        'diagnostics':{'head_min_dimension_px':100, 'paired_colour':{
            'floral_union':{'lab_chroma':flower, 'n_pixels':100},
            'green_non_head_context':{'lab_chroma':green, 'n_pixels':100, 'support_status':'available'}}}}}


def test_candidate_is_not_legacy_dominant_colour_switch_or_mirror_average():
    a = profile.extract(payload(old=60))
    b = profile.extract(payload(old=90))
    assert a['uniform_floral_chroma_value'] == b['uniform_floral_chroma_value'] == 40
    assert a['uniform_minus_legacy_chroma_value'] == -20
    assert b['uniform_minus_legacy_chroma_value'] == -50


def test_failed_qc_retains_values_without_promoting_them():
    row = profile.extract(payload(status='low_confidence'))
    assert row['uniform_floral_chroma_value'] == 40
    assert all(not row[m+'_usable'] for m in profile.METRICS)


def test_missing_green_support_is_not_a_passed_control_or_zero():
    source = payload()
    source['result']['diagnostics']['paired_colour']['green_non_head_context'] = {}
    row = profile.extract(source)
    assert row['uniform_floral_chroma_usable']
    assert row['matched_floral_chroma_value'] == 40
    assert row['matched_green_chroma_value'] is None
    assert row['matched_floral_minus_green_chroma_value'] is None
    assert not row['matched_floral_chroma_usable']


@pytest.mark.parametrize('value',[None, float('nan'), float('inf'), -1, 182, True])
def test_invalid_chroma_is_not_eligible(value):
    row = profile.extract(payload(flower=value))
    assert not row['uniform_floral_chroma_usable']


def test_empty_error_condition_still_has_all_candidate_slots():
    row = profile.extract({'status':'worker_error'})
    assert all(m+'_value' in row and not row[m+'_usable'] for m in profile.METRICS)


def test_duplicate_registered_chroma_is_rejected():
    source = payload()
    source['result']['endpoints'] *= 2
    with pytest.raises(ValueError, match='Duplicate'):
        profile.extract(source)


def test_uniform_formula_on_synthetic_grayscale_patch_and_empty_mask():
    patch = np.full((10, 10, 3), 180, dtype=np.uint8)
    assert colour_summary(patch, np.ones((10,10), dtype=bool))['lab_chroma'] == 0
    assert colour_summary(patch, np.zeros((10,10), dtype=bool))['lab_chroma'] is None


def table():
    rows = []
    for i, (image, component) in enumerate([('a','a'),('a','a'),('b','b')]):
        for condition in ('baseline','changed'):
            result = profile.extract(payload(flower=40 if condition=='baseline' else (42 if component=='a' else 50)))
            rows.append({'head_id':str(i),'sha256':image,'component_id':component,
                         'development_exposed':int(component=='a'), 'condition':condition, **result})
    return pd.DataFrame(rows)


def test_component_weighting_does_not_overweight_many_heads_in_one_image():
    records = profile.summarize(table(), ['baseline','changed'])
    r = next(r for r in records if r['metric_id']=='uniform_floral_chroma' and r['condition']=='changed' and r['exposure_stratum']=='all_cached')
    assert r['mean_signed_change'] == 6  # mean of image means 2 and 10, not 14/3
    assert r['scheduled_heads'] == 3 and r['paired_components'] == 2
    assert r['spearman_paired_component_means'] is None


def test_total_loss_keeps_original_denominator_and_missing_estimates():
    data = table()
    for metric in profile.METRICS:
        data.loc[data.condition.eq('changed'), metric+'_usable'] = False
    records = profile.summarize(data, ['baseline','changed'])
    r = next(r for r in records if r['metric_id']=='uniform_floral_chroma' and r['condition']=='changed' and r['exposure_stratum']=='all_cached')
    assert r['lost_usable_heads'] == 3 and r['paired_usable_heads'] == 0
    assert r['component_weighted_loss_fraction'] == 1
    assert r['mean_signed_change'] is None


def test_missing_or_duplicate_condition_cannot_become_complete():
    data = table()
    with pytest.raises(ValueError, match='denominator'):
        profile.summarize(data.iloc[:-1], ['baseline','changed'])
    with pytest.raises(ValueError, match='Duplicate'):
        profile.summarize(pd.concat([data, data.iloc[:1]]), ['baseline','changed'])


def test_recipe_preserves_registered_values_and_ecological_hold():
    recipe = json.loads(profile.RECIPE.read_text())
    assert tuple(recipe['metrics']) == profile.METRICS
    assert recipe['replaces_registered_chroma'] is False
    assert recipe['ecological_fitting_authorized'] is False

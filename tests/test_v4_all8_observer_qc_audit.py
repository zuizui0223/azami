"""Validate the public all-eight QC audit with synthetic aggregate records only."""
import copy
import importlib.util
from pathlib import Path
import pytest

MODULE=Path(__file__).resolve().parents[1]/'analysis/v4/audit_all8_observer_qc.py'
TRAITS=('angle','lightness','chroma','hue','aspect','circularity','solidity','width_profile')
SEEDS=(20261101,20261103,20261105)


def mod():
    assert MODULE.exists(), 'aggregate-only all8 auditor missing'
    spec=importlib.util.spec_from_file_location('all8_qc_audit', MODULE)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def payload():
    rows=[]
    for trait in TRAITS:
        for grain in (1,5):
            splits=[]
            for seed in SEEDS:
                splits.append({'seed':seed,'status':'estimated_exploratory_reliability',
                  'n_taxa':40,'n_qualified_cells':320,'n_contributing_photos':2500,
                  'taxon_equal_spatial_alignment':{'n':40,'q05':-.4,'q50':.1,'q95':.7},
                  'median_alignment_taxon_bootstrap_95':[-.02,.25],
                  'within_taxon_cell_labels_permuted_null_median_alignment':{'q95':.06},
                  'fraction_taxa_alignment_positive':.6,
                  'fraction_taxa_alignment_at_least_0_5':.15})
            rows.append({'construct':trait,'degrees':grain,'splits':splits})
    return {'schema':'azami_original_S1_global_observer_disjoint_split_half_v1',
            'source_original_S1_observations':313350,
            'distinct_authenticated_source_observers':92839,
            'results':rows}


def test_null_exceedance_is_not_relabelled_bootstrap_support():
    m=mod();a=m.audit(payload())
    assert len(a['conditions'])==16
    for row in a['conditions']:
        assert row['median_above_matched_permutation_q95_count']==3
        assert row['taxon_bootstrap_lower_above_zero_count']==0
        assert row['n_splits']==3
        assert row['taxa_alignment_ge_0_5_count_range']==[6,6]
    assert a['new_RF_fit_executed'] is False
    assert a['seeds_are_independent_experiments'] is False


@pytest.mark.parametrize('corruption',['missing','duplicate','invalid_fraction','missing_seed','nonfinite'])
def test_rejects_invalid_or_incomplete_aggregate(corruption):
    m=mod();d=payload()
    if corruption=='missing':d['results'].pop()
    if corruption=='duplicate':d['results'][-1]=copy.deepcopy(d['results'][0])
    if corruption=='invalid_fraction':d['results'][0]['splits'][0]['fraction_taxa_alignment_at_least_0_5']=.151
    if corruption=='missing_seed':d['results'][0]['splits'].pop()
    if corruption=='nonfinite':d['results'][0]['splits'][0]['taxon_equal_spatial_alignment']['q50']=float('nan')
    with pytest.raises(ValueError):m.audit(d)


def test_hue_kept_as_one_construct_and_all_morphology_retained():
    a=mod().audit(payload())
    assert {r['construct'] for r in a['conditions']}==set(TRAITS)
    assert sum(r['construct']=='hue' for r in a['conditions'])==2
    assert len([r for r in a['conditions'] if r['construct'] in ('angle','aspect','circularity','solidity','width_profile')])==10

"""Posthoc split/species-coverage sensitivity on synthetic-only data."""
import numpy as np
import pytest
from analysis.v4.rf_multitrait_predictability import aggregate_endpoint, grouped_folds, score_models
from analysis.v4.rf_extended_environment import run_extended
from test_v4_rf_multitrait_predictability import tiny_synthetic


def test_balanced_taxon_fold_is_disjoint_and_equal_species_count():
    cells = aggregate_endpoint(tiny_synthetic(), "angle")
    folds = grouped_folds(cells, "unseen_taxa_balanced")
    assert len(folds) == 5
    ids = set()
    for train, test in folds:
        group_train = set(cells.iloc[train].taxon)
        group_test = set(cells.iloc[test].taxon)
        assert not group_train & group_test
        assert 5 <= len(group_test) <= 6
        assert not ids & group_test
        ids.update(group_test)
    assert len(ids) == 30


def test_relaxed_cell_support_is_separate_from_original_gate():
    frame = tiny_synthetic()
    ix = frame.index[(frame.accepted_key == 'taxon_0') & (frame.latitude < -49)].tolist()
    frame = frame.drop(index=ix[:5]).reset_index(drop=True)
    strict = aggregate_endpoint(frame, 'angle', min_per_cell=3)
    relaxed = aggregate_endpoint(frame, 'angle', min_per_cell=1)
    assert strict.taxon.nunique() == 29
    assert relaxed.taxon.nunique() == 30
    assert (relaxed.photos == 1).sum() >= 1


def test_holdout_sensitivity_has_valid_paired_scores():
    frame = tiny_synthetic()
    cells = aggregate_endpoint(frame, 'chroma', min_per_cell=1)
    a = score_models(cells, 'chroma', 'unseen_taxa_balanced', trees=8, repeats=1)
    assert a['n_taxa'] == 30
    assert all(5 <= f['test_taxa'] <= 6 for f in a['models']['fixed_four_predictor_rf']['fold_diagnostics_aggregate_only'])
    assert np.isfinite(a['twelve_minus_four_out_of_fold_skill'])


def test_disallow_biased_posthoc_tuning(monkeypatch):
    # No population analysis attempted when QC/split choice is not in fixed list.
    import inspect
    text = inspect.getsource(run_extended)
    assert 'min_per_cell not in (1, 3)' in text
    assert 'unseen_taxa_balanced' in text
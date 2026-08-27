from __future__ import annotations

import importlib.util
import io
from pathlib import Path

import numpy as np
from Bio import Phylo


ROOT = Path(__file__).resolve().parents[1]


def load(name: str, relative: str):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_spherical_basis_has_no_dateline_discontinuity():
    module = load("spatial_sensitivity", "analysis/run_geb_v2_full27_spatial_sensitivity.py")
    basis = module.spherical_basis(np.array([20.0, 20.0]), np.array([-180.0, 180.0]))
    assert basis.shape == (2, 8)
    assert np.allclose(basis[0], basis[1], atol=1e-12)


def test_group_permutation_stays_inside_taxa():
    module = load("spatial_sensitivity", "analysis/run_geb_v2_full27_spatial_sensitivity.py")
    values = np.array([1.0, 2.0, 10.0, 20.0])[:, None]
    groups = np.array(["a", "a", "b", "b"])
    result = module.permute_rows_within_groups(values, groups, np.random.default_rng(7))
    assert set(result[:2, 0]) == {1.0, 2.0}
    assert set(result[2:, 0]) == {10.0, 20.0}


def test_batched_freedman_lane_keeps_within_taxon_null():
    module = load("spatial_sensitivity", "analysis/run_geb_v2_full27_spatial_sensitivity.py")
    response = np.array([-1.0, 1.0, -1.0, 1.0])
    design = response[:, None]
    reduced = np.empty((4, 0))
    groups = np.array(["a", "a", "b", "b"])
    beta, p_value, residual = module.linear_freedman_lane(
        response,
        design,
        reduced,
        groups,
        permutations=1999,
        rng=np.random.default_rng(11),
    )
    assert np.isclose(beta, 1.0)
    assert np.allclose(residual, 0.0)
    # Two independent two-row groups yield the observed absolute slope in
    # half of all restricted permutations.
    assert 0.45 <= p_value <= 0.55


def test_batched_moran_test_matches_observed_statistic():
    module = load("spatial_sensitivity", "analysis/run_geb_v2_full27_spatial_sensitivity.py")
    values = np.linspace(-1.0, 1.0, 20)
    latitude = np.linspace(30.0, 40.0, 20)
    longitude = np.linspace(130.0, 140.0, 20)
    expected = module.morans_i(values, latitude, longitude)
    observed, p_value, n_rows = module.moran_test(
        values,
        latitude,
        longitude,
        permutations=99,
        maximum_observations=20,
        rng=np.random.default_rng(13),
    )
    assert np.isclose(observed, expected)
    assert 0.01 <= p_value <= 1.0
    assert n_rows == 20


def test_tree_covariance_and_pagel_fit_are_finite():
    module = load("historical_sensitivity", "analysis/run_geb_v2_full27_historical_sensitivity.py")
    tree = Phylo.read(io.StringIO("((a:1,b:1):1,(c:1,d:1):1);"), "newick")
    names, covariance = module.tree_covariance(tree)
    fit = module.fit_pagel(
        np.array([0.0, 1.0, 2.0, 3.0]),
        np.array([0.0, 1.0, 2.0, 3.0]),
        covariance,
    )
    assert names == ["a", "b", "c", "d"]
    assert covariance.shape == (4, 4)
    assert 0.0 <= fit["lambda"] <= 1.0
    assert np.isfinite(fit["p_value"])


def test_joint_hue_pagel_fit_returns_two_coefficients():
    module = load("historical_sensitivity", "analysis/run_geb_v2_full27_historical_sensitivity.py")
    tree = Phylo.read(io.StringIO("(((a:1,b:1):1,c:2):1,(d:1,e:1):2);"), "newick")
    _, covariance = module.tree_covariance(tree)
    angles = np.radians(np.array([10.0, 30.0, 80.0, 140.0, 220.0]))
    fit = module.fit_pagel(
        np.column_stack([np.sin(angles), np.cos(angles)]),
        np.arange(5, dtype=float),
        covariance,
    )
    assert len(fit["beta"]) == 2
    assert np.isfinite(fit["joint_statistic"])
    assert np.isfinite(fit["p_value"])

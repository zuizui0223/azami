"""Synthetic-only tests: never load original photographs, coordinates or keys."""
import numpy as np
import pytest

from analysis.v5.background_paired_overlap import (
    BACKGROUND_FLOOR, MIN_OCCUPANCY_BG_SUPPORT,
    block_key, correct_occupancy, overlap_matrix, regional_density,
    weighted_project_from_background,
)


def test_block_key_matches_existing_20_degree_partition():
    latitude = np.array([0.0, 1.0, 19.99, 20.01, -89.99])
    longitude = np.array([0.0, 0.1, 19.99, 20.01, -179.99])
    i = ((latitude + 90).clip(max=179.99999)).astype(int)
    j = ((longitude + 180) % 360).astype(int)
    old = (i // 20) * 18 + (j // 20)
    np.testing.assert_array_equal(block_key(latitude, longitude), old)


def test_correct_occupancy_normalises_and_accounts_for_unequal_background():
    # Equal occupied frequency on two supported bins but different availability.
    occ = np.array([[0.5, 0.5]])
    bg = np.array([[0.8, 0.2]])
    corrected, support = correct_occupancy(occ, bg)
    assert support == pytest.approx(1)
    assert corrected.sum() == pytest.approx(1)
    assert corrected[0, 1] > corrected[0, 0]
    np.testing.assert_allclose(corrected, np.array([[0.2, 0.8]]))


def test_correct_occupancy_rejects_unavailable_climate_mass():
    occ = np.array([[0.01, 0.99]])
    bg = np.array([[1, 0]])
    with pytest.raises(ValueError, match="OUTSIDE_AVAILABLE"):
        correct_occupancy(occ, bg)


def test_density_is_normalised_for_weighted_background_and_observed_cells():
    edges = [np.linspace(-3, 3, 25), np.linspace(-3, 3, 25)]
    points = np.array([[-2, -2], [1, 1], [2, 2]], float)
    obs = regional_density(points, edges, .5)
    available = regional_density(points, edges, .5, weights=[1, 2, 10])
    assert obs.shape == (24, 24)
    assert np.isfinite(obs).all() and np.isfinite(available).all()
    assert obs.sum() == pytest.approx(1)
    assert available.sum() == pytest.approx(1)
    assert not np.allclose(obs, available)


def test_identical_and_disjoint_density_overlap():
    a = np.array([[1., 0., 0.], [0., 1., 0.]])
    b = np.array([[1., 0., 0.], [0., 0., 1.]])
    d = overlap_matrix(a, b)
    assert d.shape == (2, 2)
    np.testing.assert_allclose(d, np.array([[1, 0], [0, 0]]))
    with pytest.raises(ValueError, match="UNNORMALISED"):
        overlap_matrix(a * 2, b)


def test_background_pca_keeps_2_axes_and_uses_wcvp_area_not_photo_numbers():
    import pandas as pd
    from analysis.v4.native_cached_climate import VARIABLES
    rng = np.random.default_rng(10)
    n = 80
    z = rng.normal(size=(2*n, len(VARIABLES)))
    # Both ranges have geographic support and nonzero variance.
    df = pd.DataFrame(z, columns=VARIABLES)
    df["range_state"] = ["native"]*n + ["introduced"]*n
    df["_area_weight"] = np.linspace(.1, 1., 2*n)
    project, variance = weighted_project_from_background(df)
    x = project(df[list(VARIABLES)].to_numpy())
    assert x.shape == (2*n, 2)
    assert 0 < variance <= 1
    assert np.isfinite(x).all()


def test_within_region_background_correction_not_equivalent_to_raw_occupancy():
    # Changing the available climate density can reverse which region pair
    # overlaps most despite identical sampled photo densities.
    occ_n = np.array([[.5, .5]])
    occ_i = np.array([[.5, .5]])
    bg_n = np.array([[.8, .2]])
    bg_i = np.array([[.2, .8]])
    corr_n, _ = correct_occupancy(occ_n, bg_n)
    corr_i, _ = correct_occupancy(occ_i, bg_i)
    raw = overlap_matrix(occ_n.reshape(1, -1), occ_i.reshape(1, -1))[0, 0]
    adjusted = overlap_matrix(corr_n.reshape(1, -1), corr_i.reshape(1, -1))[0, 0]
    assert raw == pytest.approx(1)
    assert adjusted == pytest.approx(.4)
    assert BACKGROUND_FLOOR > 0
    assert MIN_OCCUPANCY_BG_SUPPORT > .5

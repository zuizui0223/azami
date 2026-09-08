import numpy as np
import pytest

from analysis.v3.hierarchical_ecology import SPATIAL_TERMS, spherical_basis
from analysis.v3.nuisance_design import definition, matrix


def test_definition_matches_calendar_imaging_and_all_eight_spatial_terms():
    spec = definition()
    assert spec["status"] == "NUISANCE_DESIGN_VERIFIED_NOT_FITTED"
    assert len(spec["columns"]) == 13
    assert spec["columns"][-8:] == tuple(f"spatial_{name}" for name in SPATIAL_TERMS)
    assert not spec["ecological_fitting_authorized"]


def test_nuisance_matrix_uses_fixed_year_and_imaging_transforms():
    sin = np.array([0.0, 1.0, -1.0])
    cos = np.array([1.0, 0.0, 0.0])
    year = np.array([1975.0, 2010.0, 2025.0])
    lat = np.array([35.0, 45.0, 55.0])
    lon = np.array([-100.0, 10.0, 120.0])
    size = np.array([10.0, 100.0, 1000.0])
    sharpness = np.array([0.0, 9.0, 99.0])
    out, names = matrix(sin_doy=sin, cos_doy=cos, observed_year=year,
                        latitude=lat, longitude=lon, size=size, sharpness=sharpness)
    assert out.shape == (3, 13)
    assert names[:5] == (
        "calendar_sin_doy", "calendar_cos_doy", "observation_year_decade",
        "log_head_min_dimension_px", "log1p_head_laplacian_variance",
    )
    assert np.allclose(out[:, 0], sin)
    assert np.allclose(out[:, 1], cos)
    assert np.allclose(out[:, 2], [-3.5, 0.0, 1.5])
    assert np.allclose(out[:, 3], np.log(size))
    assert np.allclose(out[:, 4], np.log1p(sharpness))
    assert np.allclose(out[:, 5:], spherical_basis(lat, lon))


def test_nuisance_matrix_never_silently_imputes_or_drops_invalid_quality():
    common = dict(sin_doy=[0.0, 0.1], cos_doy=[1.0, 0.9], observed_year=[2010, 2011],
                  latitude=[40, 41], longitude=[10, 11])
    with pytest.raises(ValueError, match="finite"):
        matrix(**common, size=[100, np.nan], sharpness=[1, 2])
    with pytest.raises(ValueError, match="positive"):
        matrix(**common, size=[100, 0], sharpness=[1, 2])
    with pytest.raises(ValueError, match="negative"):
        matrix(**common, size=[100, 200], sharpness=[1, -1])

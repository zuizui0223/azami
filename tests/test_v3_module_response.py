import numpy as np
import pandas as pd
import pytest

from analysis.v3.module_response import definition, matrix


def test_primary_geometry_is_exactly_one_eight_four_coordinates():
    spec = definition()
    assert spec["status"] == "PRIMARY_MODULE_RESPONSE_GEOMETRY_VERIFIED_NOT_FITTED"
    assert len(spec["modules"]["orientation"]["coordinates"]) == 1
    assert len(spec["modules"]["visible_colour"]["coordinates"]) == 8
    assert len(spec["modules"]["gross_shape"]["coordinates"]) == 4
    assert not spec["ecological_fitting_authorized"]


def test_colour_uses_joint_hue_and_hellinger_closed_composition_without_renormalization():
    frame = pd.DataFrame({
        "corolla_lab_lightness": [50.0, 60.0],
        "corolla_lab_chroma": [20.0, 30.0],
        "corolla_hue_sin": [0.0, 0.6],
        "corolla_hue_cos": [1.0, 0.8],
        "corolla_white_pixel_fraction": [0.25, 0.0],
        "corolla_redmagenta_pixel_fraction": [0.25, 0.25],
        "corolla_purple_pixel_fraction": [0.25, 0.25],
        "corolla_yellow_pixel_fraction": [0.25, 0.5],
    })
    out, names = matrix(frame, "visible_colour")
    assert out.shape == (2, 8)
    assert names[2:4] == ("corolla_hue_sin", "corolla_hue_cos")
    assert np.allclose(out[:, :4], frame.iloc[:, :4].to_numpy(float))
    assert np.allclose(out[:, 4:], np.sqrt(frame.iloc[:, 4:8].to_numpy(float)))


def test_colour_rejects_unclosed_or_impossible_hue_rows_instead_of_repairing_them():
    base = {
        "corolla_lab_lightness": [50.0],
        "corolla_lab_chroma": [20.0],
        "corolla_hue_sin": [0.0],
        "corolla_hue_cos": [1.0],
        "corolla_white_pixel_fraction": [0.25],
        "corolla_redmagenta_pixel_fraction": [0.25],
        "corolla_purple_pixel_fraction": [0.25],
        "corolla_yellow_pixel_fraction": [0.25],
    }
    bad = pd.DataFrame(base)
    bad.loc[0, "corolla_yellow_pixel_fraction"] = 0.2
    with pytest.raises(ValueError, match="not closed"):
        matrix(bad, "visible_colour")
    bad = pd.DataFrame(base)
    bad.loc[0, "corolla_hue_sin"] = 1.0
    bad.loc[0, "corolla_hue_cos"] = 1.0
    with pytest.raises(ValueError, match="unit-disc"):
        matrix(bad, "visible_colour")


def test_orientation_and_shape_are_identity_transforms():
    orientation = pd.DataFrame({"orientation_image_vertical_angle": [10.0, 20.0]})
    out, names = matrix(orientation, "orientation")
    assert names == ("orientation_image_vertical_angle",)
    assert np.allclose(out[:, 0], [10.0, 20.0])

    shape = pd.DataFrame({
        "capitulum_outline_aspect_ratio": [1.1, 1.2],
        "capitulum_outline_circularity": [0.8, 0.9],
        "capitulum_outline_solidity": [0.95, 0.96],
        "capitulum_width_profile_cv": [0.2, 0.3],
    })
    out, names = matrix(shape, "gross_shape")
    assert out.shape == (2, 4)
    assert np.allclose(out, shape.to_numpy(float))

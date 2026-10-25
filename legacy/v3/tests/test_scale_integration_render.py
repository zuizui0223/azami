from pathlib import Path

import numpy as np
from PIL import Image
import pytest

from reproducibility.render_scale_integration import CORE, load_inputs, render


def test_scale_integration_sources_and_metrics_are_frozen():
    data = load_inputs()
    summary = data["summary"]
    upgrade = data["upgrade"]
    estimator = data["estimator"]
    assert len(CORE) == 9
    assert len(data["pairwise"]) == 36
    assert len(data["bootstrap"]) == 1000
    assert len(data["sources"]) == 7
    assert upgrade["common_cohort"]["observations"] == 1734
    assert upgrade["common_cohort"]["taxa"] == 42
    assert upgrade["common_cohort_matrix_alignment"]["rho"] == pytest.approx(0.4391248391248392)
    assert upgrade["common_cohort_matrix_alignment"]["qap_p_one_sided"] == pytest.approx(0.0041)
    assert summary["observed"]["relations_stronger_among"] == 33
    assert estimator["equal_n_within_resampling"]["overall_strength_gate_pass"] is True
    assert estimator["equal_n_within_resampling"]["among_minus_within_median_rv_low95"] > 0
    assert estimator["equal_n_within_resampling"]["probability_positive_median_difference"] == pytest.approx(0.983)
    assert estimator["equal_n_within_resampling"]["relations_stronger_among_median"] == pytest.approx(23.0)
    assert estimator["permutation_null_centering"]["relations_stronger_among_after_null_centering"] == 19
    assert estimator["joint_coordinate_standardization"]["relations_stronger_among"] == 34


def test_scale_integration_renderer_writes_figure_and_provenance(tmp_path: Path):
    result = render(tmp_path)
    png = tmp_path / "Figure_v3_scale_integration.png"
    pdf = tmp_path / "Figure_v3_scale_integration.pdf"
    provenance = tmp_path / "Figure_v3_scale_integration_provenance.json"
    assert png.is_file() and png.stat().st_size > 50_000
    assert pdf.is_file() and pdf.stat().st_size > 5_000
    assert provenance.is_file()
    assert result["schema_version"] == 2
    assert result["common_cohort"]["observations"] == 1734
    assert result["headline_metrics"]["raw_relations_stronger_among"] == 33
    assert result["headline_metrics"]["equal_n_difference_low95"] > 0
    assert result["headline_metrics"]["equal_n_probability_positive"] == pytest.approx(0.983)
    assert result["headline_metrics"]["equal_n_relations_stronger_among_median"] == pytest.approx(23.0)
    assert len(result["outputs"]) == 2

    image = np.asarray(Image.open(png).convert("RGB"))
    assert np.all(image[:, :6, :] >= 245), "content touches/clips the left PNG canvas edge"
    assert np.all(image[:, -6:, :] >= 245), "content touches/clips the right PNG canvas edge"

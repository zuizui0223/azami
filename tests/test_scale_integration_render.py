from pathlib import Path

import numpy as np
from PIL import Image
import pytest

from reproducibility.render_scale_integration import CORE, load_inputs, render


def test_scale_integration_sources_and_metrics_are_frozen():
    data = load_inputs()
    summary = data["summary"]
    upgrade = data["upgrade"]
    assert len(CORE) == 9
    assert len(data["pairwise"]) == 36
    assert len(data["bootstrap"]) == 1000
    assert len(data["sources"]) == 6
    assert upgrade["common_cohort"]["observations"] == 1734
    assert upgrade["common_cohort"]["taxa"] == 42
    assert upgrade["common_cohort_matrix_alignment"]["rho"] == pytest.approx(0.4391248391248392)
    assert upgrade["common_cohort_matrix_alignment"]["qap_p_one_sided"] == pytest.approx(0.0041)
    assert summary["observed"]["relations_stronger_among"] == 33
    assert summary["taxon_bootstrap"]["probability_median_among_exceeds_within"] == pytest.approx(1.0)


def test_scale_integration_renderer_writes_figure_and_provenance(tmp_path: Path):
    result = render(tmp_path)
    png = tmp_path / "Figure_v3_scale_integration.png"
    pdf = tmp_path / "Figure_v3_scale_integration.pdf"
    provenance = tmp_path / "Figure_v3_scale_integration_provenance.json"
    assert png.is_file() and png.stat().st_size > 50_000
    assert pdf.is_file() and pdf.stat().st_size > 5_000
    assert provenance.is_file()
    assert result["common_cohort"]["observations"] == 1734
    assert result["headline_metrics"]["relations_stronger_among"] == 33
    assert result["headline_metrics"]["bootstrap_probability_among_exceeds_within"] == pytest.approx(1.0)
    assert len(result["outputs"]) == 2

    # Fixed-width manuscript figures must retain a clean outer gutter. The
    # first CI-rendered revision exposed clipped y-axis labels at the left
    # canvas edge, so make this a regression test rather than relying on visual
    # inspection alone. Anti-aliased near-white pixels are allowed.
    image = np.asarray(Image.open(png).convert("RGB"))
    assert np.all(image[:, :6, :] >= 245), "content touches/clips the left PNG canvas edge"
    assert np.all(image[:, -6:, :] >= 245), "content touches/clips the right PNG canvas edge"

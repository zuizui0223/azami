from pathlib import Path

import pytest

from reproducibility.render_layout_revisions import ROOT, render, replace_once


def test_layout_anchor_is_exact_and_fail_closed():
    assert replace_once("a OLD b", "OLD", "NEW") == "a NEW b"
    for source in ["no match", "OLD OLD"]:
        with pytest.raises(ValueError, match="exactly one"):
            replace_once(source, "OLD", "NEW")


def test_frozen_figure_archive_cannot_be_overwritten():
    for output in [ROOT / "reproducibility/figures", ROOT / "reproducibility/figures/nested"]:
        with pytest.raises(ValueError, match="frozen figure archive"):
            render(output)

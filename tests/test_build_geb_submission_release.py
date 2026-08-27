from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ANALYSIS = Path(__file__).resolve().parents[1] / "analysis"
sys.path.insert(0, str(ANALYSIS))

from geb_submission_release_common import stable_pca, symmetric_matrix  # noqa: E402
from geb_submission_release_tables import abstract_body_word_count  # noqa: E402


def test_stable_pca_matches_complete_case_svd_contract() -> None:
    wide = pd.DataFrame(
        {
            "a": [1.0, 2.0, 3.0, 4.0, np.nan],
            "b": [4.0, 3.0, 2.0, 1.0, 0.0],
            "c": [1.0, 1.5, 3.5, 5.0, 2.0],
        },
        index=pd.Index(["t1", "t2", "t3", "t4", "t5"], name="taxon_name"),
    )
    result = stable_pca(wide)
    assert result.scores.shape == (4, 4)
    assert result.loadings.shape == (3, 4)
    assert np.isclose(result.explained.sum(), 1.0)
    assert result.scores["taxon_name"].tolist() == ["t1", "t2", "t3", "t4"]


def test_symmetric_matrix_reconstructs_pairwise_values_without_filling_diagonal() -> None:
    pairs = pd.DataFrame(
        {
            "left": ["a", "a", "b"],
            "right": ["b", "c", "c"],
            "value": [0.1, 0.2, 0.3],
        }
    )
    matrix = symmetric_matrix(pairs, ["a", "b", "c"])
    assert np.isnan(np.diag(matrix)).all()
    assert np.allclose(matrix[np.triu_indices(3, 1)], [0.1, 0.2, 0.3])
    assert np.allclose(matrix, matrix.T, equal_nan=True)


def test_abstract_body_count_ignores_structure_and_metadata() -> None:
    text = """# Title and abstract

## Preferred title

**A title**

**Running title:** Short

## Abstract

### Aim

One two three.

### Location

Global.

**Keywords:** a; b
"""
    assert abstract_body_word_count(text) == 4

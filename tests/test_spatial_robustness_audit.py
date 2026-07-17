import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

MODULE_PATH = Path(__file__).parents[1] / "analysis" / "audit_spatial_robustness.py"
spec = importlib.util.spec_from_file_location("spatial_audit", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_morans_i_is_high_for_clustered_values():
    xy = np.array([[0, 0], [0, 1], [1, 0], [10, 10], [10, 11], [11, 10]], dtype=float)
    values = np.array([1, 1, 1, -1, -1, -1], dtype=float)
    assert module.morans_i(values, xy, k=2) > 0.5


def test_region_summary_counts_rows_and_taxa():
    df = pd.DataFrame({
        "region": ["A", "A", "B"],
        "taxon": ["x", "y", "x"],
    })
    out = module.region_summary(df, "region", "taxon").set_index("region")
    assert out.loc["A", "n_rows"] == 2
    assert out.loc["A", "n_taxa"] == 2


def test_leave_one_region_out_reports_shared_taxa():
    df = pd.DataFrame({
        "region": ["A", "A", "A", "B", "B", "B"],
        "taxon": ["x", "y", "z", "x", "y", "z"],
        "value": [1, 2, 3, 1.1, 2.1, 3.1],
    })
    out = module.leave_one_region_out_rank(df, "region", "taxon", "value")
    assert set(out["n_shared_taxa"]) == {3}
    assert (out["spearman_rho"] > 0.99).all()

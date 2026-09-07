import importlib
from pathlib import Path
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
verification=importlib.import_module("analysis.v3.verify_perturbation_summary")


def test_comparison_distinguishes_missing_from_zero():
    row={key:"0" for key in verification.COUNTS}
    row.update({key:"" for key in verification.FLOATS})
    actual={key:None for key in row}
    verification.compare(row,actual)
    actual["mean_absolute_change"]=0
    with pytest.raises(ValueError,match="availability"):
        verification.compare(row,actual)


@pytest.mark.parametrize("field",verification.COUNTS+verification.FLOATS)
def test_comparison_detects_each_checked_error(field):
    row={key:"1" for key in verification.COUNTS+verification.FLOATS}
    actual={key:1 for key in row}
    actual[field]=2
    with pytest.raises(ValueError,match="mismatch"):
        verification.compare(row,actual)


def test_independent_sql_verifies_complete_missing_grid(tmp_path):
    import test_v3_perturbation_summary as fixtures
    roots,paths,receipt=fixtures.complete_inputs(tmp_path)
    fixtures.summary.run(*roots,tmp_path/"summary")
    result=verification.verify(*roots,tmp_path/"summary",tmp_path/"verification")
    assert result["scalar_rows_checked"]==1134
    assert all(v==0 for v in result["max_absolute_numerical_differences"].values())
    assert "rank correlations" in result["not_independently_recomputed_here"]

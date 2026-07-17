import importlib.util
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "prepare_spatial_diagnostic_input.py"
spec = importlib.util.spec_from_file_location("prepare_spatial", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_require_columns_rejects_missing():
    frame = pd.DataFrame({"a": [1]})
    with pytest.raises(ValueError, match="missing required columns"):
        module.require_columns(frame, ["a", "b"], "table")


def test_assert_unique_rejects_duplicate_keys():
    frame = pd.DataFrame({"id": ["x", "x"]})
    with pytest.raises(ValueError, match="duplicate keys"):
        module.assert_unique(frame, ["id"], "table")


def test_assert_unique_accepts_endpoint_specific_rows():
    frame = pd.DataFrame(
        {
            "observation_id": ["x", "x"],
            "endpoint": ["orientation", "colour"],
        }
    )
    module.assert_unique(frame, ["observation_id", "endpoint"], "predictions")

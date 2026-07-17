import importlib.util
import unittest
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "prepare_spatial_diagnostic_input.py"
spec = importlib.util.spec_from_file_location("prepare_spatial", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


class TestPrepareSpatialDiagnosticInput(unittest.TestCase):
    def test_require_columns_rejects_missing(self):
        frame = pd.DataFrame({"a": [1]})
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            module.require_columns(frame, ["a", "b"], "table")

    def test_assert_unique_rejects_duplicate_keys(self):
        frame = pd.DataFrame({"id": ["x", "x"]})
        with self.assertRaisesRegex(ValueError, "duplicate keys"):
            module.assert_unique(frame, ["id"], "table")

    def test_assert_unique_accepts_endpoint_specific_rows(self):
        frame = pd.DataFrame(
            {
                "observation_id": ["x", "x"],
                "endpoint": ["orientation", "colour"],
            }
        )
        module.assert_unique(frame, ["observation_id", "endpoint"], "predictions")


if __name__ == "__main__":
    unittest.main()

import unittest

import pandas as pd

from azami_ch1.tabular import assert_unique, require_columns, require_complete_text


class TestTabularHelpers(unittest.TestCase):
    def test_require_columns_reports_missing(self):
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            require_columns(pd.DataFrame({"a": [1]}), ["a", "b"], "table")

    def test_assert_unique_reports_duplicate_examples(self):
        with self.assertRaisesRegex(ValueError, "duplicate keys"):
            assert_unique(pd.DataFrame({"id": ["x", "x"]}), ["id"], "table")

    def test_require_complete_text_rejects_blank_rows(self):
        frame = pd.DataFrame({"region": ["A", " ", None]})
        with self.assertRaisesRegex(ValueError, "2 rows lack reviewed region"):
            require_complete_text(frame, "region", "reviewed region")


if __name__ == "__main__":
    unittest.main()

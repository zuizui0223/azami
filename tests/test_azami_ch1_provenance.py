import json
import tempfile
import unittest
from pathlib import Path

from azami_ch1.provenance import utc_now_iso, write_json


class TestProvenance(unittest.TestCase):
    def test_utc_timestamp_is_explicit(self):
        self.assertTrue(utc_now_iso().endswith("+00:00"))

    def test_write_json_creates_parent_and_final_newline(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nested" / "value.json"
            write_json(path, {"label": "アザミ"})
            text = path.read_text(encoding="utf-8")
            self.assertTrue(text.endswith("\n"))
            self.assertEqual(json.loads(text)["label"], "アザミ")

    def test_generated_timestamp_is_optional(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "value.json"
            write_json(path, {"status": "ok"}, include_generated_utc=True)
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.assertIn("generated_utc", payload)


if __name__ == "__main__":
    unittest.main()

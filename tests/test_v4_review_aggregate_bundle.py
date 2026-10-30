"""Behavioral checks for the narrowly scoped aggregate release builder."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile


BUILDER = (
    Path(__file__).resolve().parents[1]
    / "manuscripts/v4_s1_s2_frozen_20261025/build_review_aggregate_bundle.py"
)


class ReviewAggregateBundleTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(BUILDER.is_file(), "aggregate builder has not been implemented")
        spec = importlib.util.spec_from_file_location("review_aggregate_builder", BUILDER)
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)

    def test_nested_observation_coordinates_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "protected"):
            self.module.check_data("aggregate.json", json.dumps({
                "results": [{"beta": 0.04, "latitude": 38.1}]
            }).encode())

    def test_observation_links_are_rejected_without_an_id_column(self):
        with self.assertRaisesRegex(ValueError, "observation"):
            self.module.check_data("aggregate.json", json.dumps({
                "note": "https://www.inaturalist.org/observations/123456"
            }).encode())

    def test_uploader_names_are_rejected_without_numeric_ids(self):
        with self.assertRaisesRegex(ValueError, "protected"):
            self.module.check_data("aggregate.json", json.dumps({
                "summary": {"uploader_name": "Example Person"}
            }).encode())

    def test_environment_provider_urls_and_trait_coordinates_are_permitted(self):
        raw = json.dumps({"coordinate": "corolla_lab_chroma", "beta": 0.1345,
                          "source_url": "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip"}).encode()
        self.module.check_data("aggregate.json", raw)

    def test_identity_redaction_preserves_all_scientific_numbers(self):
        original = {"beta": 0.04544684307686103, "n": 168907,
                    "source": {"run_id": 40841007115,
                               "run_url": "https://github.com/zuizui0223/azami/actions/runs/40841007115"},
                    "p": [0.009206350580895119, 1.0]}
        ledger = []
        result = self.module.review_json(original, {40841007115: "audit-001"}, {}, ledger)
        self.assertEqual(result["beta"], 0.04544684307686103)
        self.assertEqual(result["n"], 168907)
        self.assertEqual(result["p"], [0.009206350580895119, 1.0])
        self.assertEqual(result["source"]["run_id"], "audit-001")
        self.assertNotIn("github.com", json.dumps(result))
        self.assertNotIn("zuizui0223", json.dumps(result))
        self.assertTrue(ledger)
        self.assertEqual(original["source"]["run_id"], 40841007115)

    def test_changed_numeric_result_is_detected(self):
        with self.assertRaisesRegex(ValueError, "scientific"):
            self.module.assert_scientific_numbers({"beta": 0.04}, {"beta": 0.05})

    def test_unexpected_zip_member_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("aggregate.json", '{"beta":0.04}')
                archive.writestr("observations.csv", "value\n1\n")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, "members"):
                self.module.read_archive(path, digest, {"aggregate.json"})

    def test_changed_archive_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("aggregate.json", '{"beta":0.04}')
            with self.assertRaisesRegex(ValueError, "SHA"):
                self.module.read_archive(path, "0" * 64, {"aggregate.json"})


if __name__ == "__main__":
    unittest.main()

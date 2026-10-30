"""Tests for the no-publish Zenodo V4 aggregate handoff checks."""

import hashlib
import importlib.util
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/validate_zenodo_v4_release.py"
META = ROOT / "manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/ZENODO_V4_METADATA_WORKING_DRAFT.json"
spec = importlib.util.spec_from_file_location("v4_zenodo_guard", SCRIPT)
guard = importlib.util.module_from_spec(spec)
spec.loader.exec_module(guard)


class ZenodoGuardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.meta = json.loads(META.read_text(encoding="utf-8"))

    def test_working_draft_is_explicitly_blocked(self):
        report = guard.check_metadata(self.meta)
        self.assertEqual(report["status"], "PREPARED_NOT_PUBLISHED")
        self.assertEqual(len(report["blockers"]), 5)
        self.assertFalse(report["doi_assigned"])
        self.assertTrue(all(not p["recommended_publication"] for p in self.meta["candidate_files"]))

    def test_no_existing_v2_doi_as_v4(self):
        self.assertIsNone(self.meta["published_doi"])
        self.assertNotEqual(self.meta.get("version_label"), "v2")
        self.assertIn("10.5281/zenodo.22295791", self.meta["prior_archive_note"])

    def _make_zip(self, root: Path, name: str, entry: str, payload: bytes = b"aggregate"):
        path = root / name
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr(entry, payload)
        altered = json.loads(json.dumps(self.meta))
        altered["candidate_files"][1]["name"] = name
        altered["candidate_files"][1]["bytes"] = path.stat().st_size
        altered["candidate_files"][1]["members"] = 1
        altered["candidate_files"][1]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path, altered

    def test_valid_exact_bundle_is_identity_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._make_zip(Path(directory), "synthetic-review.zip", "safe_data.csv")
            result = guard.inspect_bundle(path, meta)
            self.assertEqual(result["member_count"], 1)

    def test_wrong_sha_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._make_zip(Path(directory), "synthetic-review.zip", "safe_data.csv")
            meta["candidate_files"][1]["sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                guard.inspect_bundle(path, meta)

    def test_parent_path_rejected_even_with_correct_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._make_zip(Path(directory), "synthetic-review.zip", "../private.csv")
            with self.assertRaisesRegex(ValueError, "Unsafe"):
                guard.inspect_bundle(path, meta)

    def test_private_provenance_rejected_even_with_correct_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            path, meta = self._make_zip(Path(directory), "synthetic-review.zip", "internal_provenance.json")
            with self.assertRaisesRegex(ValueError, "identifying"):
                guard.inspect_bundle(path, meta)

    def test_publication_gate_requires_file_and_approvals(self):
        self.assertEqual(guard.main(["--metadata", str(META)]), 0)
        self.assertEqual(guard.main(["--metadata", str(META), "--require-ready"]), 2)

    def test_script_has_no_mutating_network_or_publish_calls(self):
        content = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("requests.post(", content)
        self.assertNotIn("urlopen(", content)
        self.assertNotIn("zenodo.org/api/", content)
        self.assertNotIn("actions/publish", content)


if __name__ == "__main__":
    unittest.main()

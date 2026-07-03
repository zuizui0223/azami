from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "ch1_global" / "v2" / "00_build_annotation_manifest.py"
FIXTURE = ROOT / "tests" / "fixtures" / "ch1_crop_metadata_small.csv"


class TestAnnotationManifest(unittest.TestCase):
    def test_manifest_is_observation_deduplicated_and_cv_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "out"
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input", str(FIXTURE),
                    "--out-dir", str(out_dir),
                    "--seed", "123",
                    "--min-heads-per-species", "2",
                    "--max-heads-per-species", "3",
                    "--n-main", "10",
                    "--n-calibration", "4",
                    "--n-segmentation", "3",
                    "--n-cv-folds", "3",
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            self.assertIn("[OK] Wrote Ch.1 annotation manifest", result.stdout)

            manifest = pd.read_csv(out_dir / "annotation_manifest.csv")
            assignments = pd.read_csv(out_dir / "annotation_assignments.csv")
            run_info = pd.read_csv(out_dir / "annotation_manifest_run_info.csv")

            # Low-confidence and obscured records are excluded; only two species meet the quota.
            self.assertEqual(len(manifest), 6)
            self.assertSetEqual(set(manifest["species"]), {"Cirsium a", "Cirsium b"})
            self.assertTrue(manifest["coordinate_usable_for_environment"].all())

            # Annotation units are one per observation, and group keys survive for later CV/LOSO.
            self.assertTrue(manifest["observation_group"].is_unique)
            self.assertTrue(manifest["annotation_unit_id"].is_unique)
            self.assertTrue(manifest["species_cv_fold"].between(1, 3).all())
            self.assertTrue(manifest["spatial_cv_fold"].between(1, 3).all())
            self.assertTrue((manifest.groupby("species")["species_cv_fold"].nunique() == 1).all())

            self.assertEqual((manifest["annotation_batch"] == "calibration_double_label").sum(), 4)
            self.assertEqual(manifest["needs_segmentation"].sum(), 3)
            self.assertEqual(len(assignments), 10)  # 4 double-labelled + 2 single-labelled
            self.assertEqual(int(run_info.loc[0, "manifest_units"]), 6)


if __name__ == "__main__":
    unittest.main()

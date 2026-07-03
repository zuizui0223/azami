from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_SCRIPT = ROOT / "ch1_global" / "v2" / "00_build_annotation_manifest.py"
PACKET_SCRIPT = ROOT / "ch1_global" / "v2" / "01_build_blinded_annotation_packets.py"
AGREEMENT_SCRIPT = ROOT / "ch1_global" / "v2" / "02_compile_double_annotations.py"
FIXTURE = ROOT / "tests" / "fixtures" / "ch1_crop_metadata_small.csv"

MANIFEST_ARGS = [
    "--seed", "123",
    "--min-heads-per-species", "2",
    "--max-heads-per-species", "3",
    "--n-main", "10",
    "--n-calibration", "4",
    "--n-segmentation", "3",
    "--n-cv-folds", "3",
]


class TestAnnotationFoundation(unittest.TestCase):
    def build_manifest(self, out_dir: Path) -> pd.DataFrame:
        result = subprocess.run(
            [
                sys.executable,
                str(MANIFEST_SCRIPT),
                "--input", str(FIXTURE),
                "--out-dir", str(out_dir),
                *MANIFEST_ARGS,
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        self.assertIn("[OK] Wrote Ch.1 annotation manifest", result.stdout)
        return pd.read_csv(out_dir / "annotation_manifest.csv")

    def test_manifest_is_observation_deduplicated_and_cv_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "manifest"
            manifest = self.build_manifest(out_dir)
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

    def test_blinded_packets_and_adjudication_queue(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            manifest_dir = tmp_path / "manifest"
            manifest = self.build_manifest(manifest_dir)
            packet_dir = tmp_path / "packets"
            subprocess.run(
                [
                    sys.executable,
                    str(PACKET_SCRIPT),
                    "--manifest", str(manifest_dir / "annotation_manifest.csv"),
                    "--out-dir", str(packet_dir),
                    "--seed", "123",
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            primary = pd.read_csv(packet_dir / "annotator_1_packet.csv")
            secondary = pd.read_csv(packet_dir / "annotator_2_packet.csv")
            self.assertEqual(len(primary), len(manifest))
            self.assertEqual(len(secondary), 4)
            self.assertFalse({"species", "latitude", "longitude", "spatial_block"}.intersection(primary.columns))

            shared_ids = sorted(set(primary["annotation_unit_id"]).intersection(secondary["annotation_unit_id"]))
            defaults = {
                "flowering_qc": "usable",
                "detector_box_qc": "accept",
                "camera_tilt_flag": "vertical_reference_clear",
                "view_angle_class": "lateral",
                "head_orientation_class": "upright",
                "flower_colour_class": "white",
                "colour_qc": "usable",
                "involucral_cover_score_0_to_4": "2",
                "head_shape_class": "globose",
                "corolla_display_score_0_to_4": "3",
                "annotation_complete": "yes",
            }
            for frame in (primary, secondary):
                for column, value in defaults.items():
                    frame[column] = value
            # One deliberate conflict confirms that the compiler produces an adjudication row.
            secondary.loc[secondary["annotation_unit_id"].eq(shared_ids[0]), "head_orientation_class"] = "nodding"
            primary_path = tmp_path / "annotator_1_completed.csv"
            secondary_path = tmp_path / "annotator_2_completed.csv"
            primary.to_csv(primary_path, index=False)
            secondary.to_csv(secondary_path, index=False)

            agreement_dir = tmp_path / "agreement"
            subprocess.run(
                [
                    sys.executable,
                    str(AGREEMENT_SCRIPT),
                    "--primary", str(primary_path),
                    "--secondary", str(secondary_path),
                    "--out-dir", str(agreement_dir),
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            summary = pd.read_csv(agreement_dir / "annotation_agreement_summary.csv")
            conflicts = pd.read_csv(agreement_dir / "annotation_adjudication_queue.csv")
            orientation = summary.loc[summary["trait"].eq("head_orientation_class")].iloc[0]
            self.assertEqual(int(orientation["n_joint_completed"]), 4)
            self.assertAlmostEqual(float(orientation["percent_agreement"]), 0.75)
            self.assertEqual(len(conflicts), 1)
            self.assertEqual(conflicts.loc[0, "trait"], "head_orientation_class")


if __name__ == "__main__":
    unittest.main()

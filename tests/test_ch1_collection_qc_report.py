from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "ch1_global" / "v2" / "08_summarize_metadata_collection.py"
SPEC = importlib.util.spec_from_file_location("ch1_collection_qc_report", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


class TestCollectionQCReport(unittest.TestCase):
    def test_writes_metrics_tables_and_pilot_warning(self) -> None:
        table = pd.DataFrame(
            {
                "obs_id": ["1", "2", "3"],
                "photo_id": ["11", "12", "13"],
                "taxon_name": ["Cirsium alpha", "Cirsium beta", "Cirsium beta"],
                "taxon_rank": ["species", "species", "genus"],
                "coordinate_usable_for_environment": [True, False, True],
                "inat_flowers_annotation": [True, False, False],
                "quality_grade": ["research", "needs_id", "casual"],
                "photo_license_code": ["cc-by", "", "cc0"],
                "geoprivacy": ["", "obscured", ""],
            }
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "photo_metadata.csv"
            table.to_csv(metadata, index=False)
            queue_dir = root / "queue"
            queue_dir.mkdir()
            pd.DataFrame(
                {
                    "queue_id": ["q1", "q2"],
                    "taxon_name": ["Cirsium alpha", "Cirsium beta"],
                    "spatial_block": ["a", "b"],
                }
            ).to_csv(queue_dir / "image_screening_queue.csv", index=False)
            out_dir = root / "qc"
            with patch.object(sys, "argv", [
                "report",
                "--metadata", str(metadata),
                "--queue-dir", str(queue_dir),
                "--out-dir", str(out_dir),
                "--technical-pilot",
                "--minimum-species-for-final", "3",
            ]):
                REPORT.main()
            metrics = json.loads((out_dir / "collection_qc_metrics.json").read_text(encoding="utf-8"))
            self.assertEqual(metrics["n_photo_rows"], 3)
            self.assertEqual(metrics["n_observations"], 3)
            self.assertEqual(metrics["n_species_total"], 2)
            self.assertEqual(metrics["queue_rows"], 2)
            report = (out_dir / "collection_qc_report.md").read_text(encoding="utf-8")
            self.assertIn("technical pilot", report)
            self.assertIn("below the final-stage target", report)
            self.assertTrue((out_dir / "top_species_photo_rows.csv").exists())


if __name__ == "__main__":
    unittest.main()

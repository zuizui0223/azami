from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
APP_PATH = ROOT / "ch1_global" / "v2" / "03_trait_annotation_app.py"
SPEC = importlib.util.spec_from_file_location("ch1_trait_annotation_app", APP_PATH)
assert SPEC is not None and SPEC.loader is not None
APP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(APP)


class TestTraitAnnotationAppCore(unittest.TestCase):
    def packet(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "annotation_unit_id": ["ch1_000002", "ch1_000001"],
                "source_image": ["source2.jpg", "source1.jpg"],
                "crop_path": ["crop2.jpg", "crop1.jpg"],
                "needs_trait_labels": [True, True],
            }
        )

    def test_response_ready_requires_every_trait_to_be_considered(self) -> None:
        blank = {field: "" for field in APP.RESPONSE_FIELDS}
        self.assertFalse(APP.response_ready(blank))
        assessable = {field: options[-1] for field, options in APP.RESPONSE_FIELDS.items()}
        self.assertTrue(APP.response_ready(assessable))

    def test_existing_responses_join_by_id_and_keep_packet_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "responses.csv"
            existing = pd.DataFrame(
                {
                    "annotation_unit_id": ["ch1_000001"],
                    "flowering_qc": ["usable"],
                    "annotation_complete": ["yes"],
                    "notes": ["resumed"],
                }
            )
            existing.to_csv(path, index=False)
            merged = APP.merge_existing(self.packet(), path)
            self.assertListEqual(merged["annotation_unit_id"].tolist(), ["ch1_000002", "ch1_000001"])
            second = merged.loc[merged["annotation_unit_id"].eq("ch1_000001")].iloc[0]
            self.assertEqual(second["flowering_qc"], "usable")
            self.assertEqual(second["annotation_complete"], "yes")
            self.assertEqual(second["notes"], "resumed")
            self.assertTrue(APP.is_complete(second))

    def test_atomic_writes_create_reusable_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            csv_path = root / "nested" / "responses.csv"
            json_path = root / "nested" / "status.json"
            frame = self.packet()
            APP.atomic_write_csv(frame, csv_path)
            APP.atomic_write_json({"n_total": 2, "last_index_1_based": 1}, json_path)
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())
            self.assertEqual(len(pd.read_csv(csv_path)), 2)
            self.assertIn('"n_total": 2', json_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import importlib.util
import io
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "ch1_global" / "v2" / "31_trait_audit_app.py"
spec = importlib.util.spec_from_file_location("trait_audit_app", SCRIPT)
assert spec is not None and spec.loader is not None
APP = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = APP
spec.loader.exec_module(APP)


class TestTraitAuditCloudPacket(unittest.TestCase):
    def _task_frame(self) -> pd.DataFrame:
        return pd.DataFrame([
            {
                "task_id": "task_1",
                "annotation_unit_id": "unit_1",
                "trait_id": "capitulum_orientation",
                "source_image": "source_images/1.jpg",
                "crop_path": "head_crops/1.jpg",
                "context_crop_path": "context_crops/1.jpg",
            }
        ])

    def test_extract_and_discover_public_packet(self) -> None:
        payload = io.BytesIO()
        with zipfile.ZipFile(payload, "w") as archive:
            archive.writestr("bundle/blinded_trait_audit_tasks.csv", self._task_frame().to_csv(index=False))
            archive.writestr("bundle/source_images/1.jpg", b"image")
            archive.writestr("bundle/head_crops/1.jpg", b"image")
            archive.writestr("bundle/context_crops/1.jpg", b"image")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            APP.extract_public_packet_zip(payload.getvalue(), root)
            packet = APP.discover_public_packet(root)
            self.assertEqual(packet.task_path.name, "blinded_trait_audit_tasks.csv")
            self.assertTrue(packet.task_path.is_file())

    def test_rejects_zip_path_traversal(self) -> None:
        payload = io.BytesIO()
        with zipfile.ZipFile(payload, "w") as archive:
            archive.writestr("../private.txt", "no")
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "Unsafe path"):
                APP.extract_public_packet_zip(payload.getvalue(), Path(tmp))

    def test_rejects_private_fields_even_under_benign_filename(self) -> None:
        leaked = self._task_frame()
        leaked["ai_candidate_state"] = "upright"
        payload = io.BytesIO()
        with zipfile.ZipFile(payload, "w") as archive:
            archive.writestr("bundle/tasks.csv", leaked.to_csv(index=False))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            APP.extract_public_packet_zip(payload.getvalue(), root)
            with self.assertRaisesRegex(ValueError, "not blinded"):
                APP.discover_public_packet(root)

    def test_resume_csv_must_match_packet_tasks(self) -> None:
        response = pd.DataFrame([
            {"task_id": "other", "human_state": "upright", "task_complete": "yes", "notes": ""}
        ]).to_csv(index=False).encode()
        with self.assertRaisesRegex(ValueError, "absent from this packet"):
            APP.validate_resume_csv(response, {"task_1"})


if __name__ == "__main__":
    unittest.main()

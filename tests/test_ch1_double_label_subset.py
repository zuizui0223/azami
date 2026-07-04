from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILD_SUBSET = load_module("34_build_blinded_double_label_subset.py", "ch1_double_label_subset")


class TestDoubleLabelSubset(unittest.TestCase):
    def test_subset_is_exactly_designated_and_public_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            public_root = root / "public"
            for directory in ("source_images", "head_crops", "context_crops"):
                (public_root / directory).mkdir(parents=True, exist_ok=True)
            rows = []
            for index in range(1, 4):
                for relative in (
                    f"source_images/source_{index}.jpg",
                    f"head_crops/head_{index}.jpg",
                    f"context_crops/context_{index}.jpg",
                ):
                    Image.new("RGB", (16, 12), color=(index, index, index)).save(public_root / relative)
                rows.append({
                    "task_id": f"task_{index}", "annotation_unit_id": f"unit_{index}",
                    "trait_id": "capitulum_orientation", "source_image": f"source_images/source_{index}.jpg",
                    "crop_path": f"head_crops/head_{index}.jpg", "context_crop_path": f"context_crops/context_{index}.jpg",
                })
            pd.DataFrame(rows).to_csv(public_root / "blinded_trait_audit_tasks.csv", index=False)
            key = root / "private_key.csv"
            pd.DataFrame([
                {"task_id": "task_1", "trait_id": "capitulum_orientation", "double_label": "true", "ai_candidate_state": "upright"},
                {"task_id": "task_2", "trait_id": "capitulum_orientation", "double_label": "false", "ai_candidate_state": "nodding"},
                {"task_id": "task_3", "trait_id": "capitulum_orientation", "double_label": "true", "ai_candidate_state": "upright"},
            ]).to_csv(key, index=False)
            output = root / "out"
            with patch.object(sys, "argv", [
                "subset", "--public-packet-root", str(public_root), "--private-audit-key", str(key),
                "--out-dir", str(output), "--batch-name", "fixture",
            ]):
                BUILD_SUBSET.main()

            subset_root = output / "blinded_second_annotator_packet"
            tasks = pd.read_csv(subset_root / "blinded_trait_audit_tasks.csv", dtype=str, keep_default_na=False)
            self.assertEqual(tasks["task_id"].tolist(), ["task_1", "task_3"])
            self.assertEqual(
                set(tasks.columns),
                {"task_id", "annotation_unit_id", "trait_id", "source_image", "crop_path", "context_crop_path"},
            )
            self.assertEqual(
                {"ai_candidate_state", "double_label", "selection_stratum"}.intersection(tasks.columns),
                set(),
            )
            for field in ("source_image", "crop_path", "context_crop_path"):
                self.assertTrue(all((subset_root / value).is_file() for value in tasks[field]))
            report = json.loads((output / "blinded_second_annotator_subset_report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["n_tasks"], 2)
            self.assertEqual(report["n_unique_annotation_units"], 2)


if __name__ == "__main__":
    unittest.main()

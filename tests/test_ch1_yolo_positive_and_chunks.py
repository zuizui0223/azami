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
V2 = ROOT / "ch1_global" / "v2"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


POOL = load_module("09_build_yolo_positive_pool.py", "ch1_yolo_positive_pool")
MERGE = load_module("10_merge_metadata_chunks.py", "ch1_merge_metadata_chunks")


class TestYOLOPositivePoolAndChunkMerge(unittest.TestCase):
    def test_positive_pool_keeps_photo_level_detected_head_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            queue = root / "queue.csv"
            screening = root / "screening.csv"
            crops = root / "crops.csv"
            pd.DataFrame([
                {"queue_id": "q1", "obs_id": "o1", "photo_id": "p1", "taxon_name": "Cirsium alpha"},
                {"queue_id": "q2", "obs_id": "o2", "photo_id": "p2", "taxon_name": "Cirsium beta"},
            ]).to_csv(queue, index=False)
            pd.DataFrame([
                {"queue_id": "q1", "screen_status": "detected", "n_detections": 2, "max_yolo_conf": 0.93},
                {"queue_id": "q2", "screen_status": "no_detection", "n_detections": 0, "max_yolo_conf": ""},
            ]).to_csv(screening, index=False)
            pd.DataFrame([
                {"queue_id": "q1", "yolo_conf": 0.93, "crop_path": "crop1.jpg"},
                {"queue_id": "q1", "yolo_conf": 0.71, "crop_path": "crop2.jpg"},
            ]).to_csv(crops, index=False)
            out = root / "out"
            with patch.object(sys, "argv", ["pool", "--queue", str(queue), "--screening", str(screening), "--crops", str(crops), "--out-dir", str(out), "--min-yolo-conf", "0.8"]):
                POOL.main()
            image_summary = pd.read_csv(out / "yolo_image_summary.csv")
            heads = pd.read_csv(out / "yolo_positive_heads.csv")
            payload = json.loads((out / "yolo_positive_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(int(image_summary.loc[image_summary.queue_id.eq("q1"), "visible_head_count_yolo"].iloc[0]), 2)
            self.assertEqual(len(heads), 1)
            self.assertEqual(int(heads.loc[0, "source_image_head_count"]), 2)
            self.assertEqual(payload["n_yolo_positive_images"], 1)
            self.assertEqual(payload["n_yolo_positive_head_crops"], 1)

    def test_chunk_merger_deduplicates_photo_ids_and_retains_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for folder, rows in {
                "chunk_a": [
                    {"obs_id": "1", "photo_id": "11", "taxon_name": "Cirsium alpha", "taxon_rank": "species", "source_file_stem": "inat_photo_11"},
                    {"obs_id": "2", "photo_id": "12", "taxon_name": "Cirsium beta", "taxon_rank": "species", "source_file_stem": "inat_photo_12"},
                ],
                "chunk_b": [
                    {"obs_id": "2", "photo_id": "12", "taxon_name": "Cirsium beta", "taxon_rank": "species", "source_file_stem": "inat_photo_12"},
                    {"obs_id": "3", "photo_id": "13", "taxon_name": "Cirsium gamma", "taxon_rank": "species", "source_file_stem": "inat_photo_13"},
                ],
            }.items():
                directory = root / folder
                directory.mkdir()
                pd.DataFrame(rows).to_csv(directory / "photo_metadata.csv", index=False)
            out = root / "merged"
            with patch.object(sys, "argv", ["merge", "--input-dir", str(root), "--out-dir", str(out)]):
                MERGE.main()
            merged = pd.read_csv(out / "photo_metadata_merged.csv", dtype=str)
            provenance = json.loads((out / "merge_provenance.json").read_text(encoding="utf-8"))
            self.assertEqual(len(merged), 3)
            self.assertEqual(provenance["n_duplicate_photo_rows_removed"], 1)
            self.assertIn("metadata_chunk_source", merged.columns)


if __name__ == "__main__":
    unittest.main()

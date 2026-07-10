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


ADAPT = load_module("50_build_representative_holdout_inputs.py", "ch1_holdout_inputs")
BUILD = load_module("47_build_representative_trait_holdout.py", "ch1_build_holdout_from_inputs")


class TestHoldoutInputs(unittest.TestCase):
    def test_join_drops_unassessable_and_missing_block_then_feeds_47(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # head-level table: 3 units, one unassessable, one whose obs has no block
            head = pd.DataFrame([
                {"annotation_unit_id": "u1", "trait_id": "corolla_colour_class", "taxon_name": "Cirsium_a",
                 "obs_id": "o1", "analysis_state_ai_conservative": "purple_blue", "analysis_state_ai_all": "purple_blue"},
                {"annotation_unit_id": "u2", "trait_id": "corolla_colour_class", "taxon_name": "Cirsium_b",
                 "obs_id": "o2", "analysis_state_ai_conservative": "unassessable", "analysis_state_ai_all": "pink"},
                {"annotation_unit_id": "u3", "trait_id": "corolla_colour_class", "taxon_name": "Cirsium_c",
                 "obs_id": "o3", "analysis_state_ai_conservative": "white_or_cream", "analysis_state_ai_all": "white_or_cream"},
            ])
            head_csv = root / "head.csv"
            head.to_csv(head_csv, index=False)
            # environment: o1/o3 have blocks, o2 missing entirely
            env = pd.DataFrame([
                {"obs_id": "o1", "spatial_block_10deg": "blkA"},
                {"obs_id": "o3", "spatial_block_10deg": "blkB"},
            ])
            env_csv = root / "env.csv"
            env.to_csv(env_csv, index=False)

            out = root / "inputs"
            with patch.object(sys, "argv", [
                "adapt", "--head-traits", str(head_csv), "--environment", str(env_csv),
                "--out-dir", str(out), "--measurement-mode", "conservative",
            ]):
                ADAPT.main()
            candidate = pd.read_csv(out / "candidate_units.csv")
            report = json.loads((out / "candidate_units_report.json").read_text())

            # u2 dropped (unassessable); all survivors have the 47-required columns and a block
            self.assertEqual(set(candidate["annotation_unit_id"]), {"u1", "u3"})
            self.assertEqual(
                list(candidate.columns),
                ["annotation_unit_id", "trait_id", "taxon_name", "spatial_block_10deg", "ai_candidate_state"],
            )
            self.assertEqual(report["n_after_dropping_unassessable"], 2)

            # the output feeds 47 without schema errors
            build_out = root / "holdout"
            with patch.object(sys, "argv", [
                "build", "--candidate-units", str(out / "candidate_units.csv"), "--out-dir", str(build_out),
                "--batch-name", "t", "--n-tasks-per-trait", "2", "--min-block-taxa", "1",
            ]):
                BUILD.main()
            key = pd.read_csv(build_out / "representative_trait_holdout_key" / "representative_trait_holdout_key.csv")
            self.assertGreaterEqual(len(key), 1)

    def test_packet_source_normalised_to_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            head = pd.DataFrame([
                {"annotation_unit_id": "u1", "trait_id": "corolla_colour_class", "taxon_name": "Cirsium_a",
                 "obs_id": "o1", "analysis_state_ai_conservative": "purple_blue", "analysis_state_ai_all": "purple_blue"},
            ])
            env = pd.DataFrame([{"obs_id": "o1", "spatial_block_10deg": "blkA"}])
            packet = pd.DataFrame([{"annotation_unit_id": "u1", "source_image": "s/u1.jpg",
                                    "crop_path": "c/u1.jpg", "context_crop_path": "x/u1.jpg", "extra": "ignored"}])
            for name, frame in (("head.csv", head), ("env.csv", env), ("packet.csv", packet)):
                frame.to_csv(root / name, index=False)
            out = root / "inputs"
            with patch.object(sys, "argv", [
                "adapt", "--head-traits", str(root / "head.csv"), "--environment", str(root / "env.csv"),
                "--out-dir", str(out), "--packet-source", str(root / "packet.csv"),
            ]):
                ADAPT.main()
            manifest = pd.read_csv(out / "packet_manifest.csv")
            self.assertEqual(list(manifest.columns),
                             ["annotation_unit_id", "source_image", "crop_path", "context_crop_path"])


if __name__ == "__main__":
    unittest.main()

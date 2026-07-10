from __future__ import annotations

import importlib.util
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


BUILD = load_module("51_build_representative_holdout_packet.py", "ch1_simple_holdout_builder")
EVAL = load_module("48_evaluate_representative_trait_holdout.py", "ch1_simple_holdout_eval")


class TestSimpleRepresentativeHoldout(unittest.TestCase):
    def test_sample_is_proportional_and_deterministic(self) -> None:
        rows = []
        for index in range(80):
            rows.append({
                "annotation_unit_id": f"common_{index:03d}",
                "taxon_name": "Cirsium_common",
                "spatial_block_10deg": "block_common",
            })
        for index in range(20):
            rows.append({
                "annotation_unit_id": f"rare_{index:03d}",
                "taxon_name": "Cirsium_rare",
                "spatial_block_10deg": "block_rare",
            })
        source = pd.DataFrame(rows)
        sample = BUILD.representative_sample(source, 20, "fixture")
        self.assertEqual(sample["taxon_name"].value_counts().to_dict(), {
            "Cirsium_common": 16,
            "Cirsium_rare": 4,
        })
        repeated = BUILD.representative_sample(source, 20, "fixture")
        self.assertEqual(list(sample["annotation_unit_id"]), list(repeated["annotation_unit_id"]))

    def test_default_gate_uses_one_sample_and_human_agreement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame([{
                "trait_id": "capitulum_orientation",
                "allow_multiple": "false",
            }]).to_csv(root / "ontology.csv", index=False)

            key_rows = []
            unit_rows = []
            annotation_rows = []
            for index in range(12):
                task_id = f"task_{index:03d}"
                unit_id = f"unit_{index:03d}"
                double = index < 4
                key_rows.append({
                    "task_id": task_id,
                    "annotation_unit_id": unit_id,
                    "trait_id": "capitulum_orientation",
                    "taxon_name": f"Cirsium_{index % 3}",
                    "double_label": str(double).lower(),
                    "ai_candidate_state": "upright",
                })
                unit_rows.append({
                    "annotation_unit_id": unit_id,
                    "trait_id": "capitulum_orientation",
                    "taxon_name": f"Cirsium_{index % 3}",
                    "analysis_state_ai_conservative": "upright",
                })
                human_state = "nodding" if index == 0 else "upright"
                annotation_rows.append({
                    "task_id": task_id,
                    "annotation_unit_id": unit_id,
                    "trait_id": "capitulum_orientation",
                    "human_state": human_state,
                    "annotator_id": "primary",
                })
                if double:
                    annotation_rows.append({
                        "task_id": task_id,
                        "annotation_unit_id": unit_id,
                        "trait_id": "capitulum_orientation",
                        "human_state": human_state,
                        "annotator_id": "secondary",
                    })

            pd.DataFrame(key_rows).to_csv(root / "key.csv", index=False)
            pd.DataFrame(unit_rows).to_csv(root / "units.csv", index=False)
            pd.DataFrame(annotation_rows).to_csv(root / "annotations.csv", index=False)
            output = root / "output"
            argv = [
                "evaluate",
                "--canonical-annotations", str(root / "annotations.csv"),
                "--holdout-key", str(root / "key.csv"),
                "--analysis-units", str(root / "units.csv"),
                "--ontology", str(root / "ontology.csv"),
                "--out-dir", str(output),
                "--primary-annotator", "primary",
                "--min-assessable-per-trait", "10",
                "--accuracy-wilson-floor", "0",
                "--human-human-wilson-floor", "0",
            ]
            with patch.object(sys, "argv", argv):
                EVAL.main()

            gate = pd.read_csv(output / "representative_holdout_trait_gate.csv")
            row = gate.iloc[0]
            self.assertEqual(int(row["n_scoreable"]), 12)
            self.assertAlmostEqual(float(row["ai_human_accuracy"]), 11 / 12)
            self.assertEqual(int(row["human_human_pairs"]), 4)
            self.assertEqual(row["gate_decision"], "accept_for_analysis")
            self.assertNotIn("selection_stratum", pd.read_csv(root / "key.csv").columns)


if __name__ == "__main__":
    unittest.main()

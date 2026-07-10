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


BUILD = load_module("51_build_representative_holdout_packet.py", "ch1_representative_primary_packet")
EVAL = load_module("48_evaluate_representative_trait_holdout.py", "ch1_representative_primary_eval")


class TestRepresentativePrimaryPacket(unittest.TestCase):
    def test_population_sample_is_proportional_across_cells(self) -> None:
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
        sample = BUILD.representative_sample(pd.DataFrame(rows), 20, "fixture")
        counts = sample["taxon_name"].value_counts().to_dict()
        self.assertEqual(len(sample), 20)
        self.assertEqual(counts, {"Cirsium_common": 16, "Cirsium_rare": 4})
        repeated = BUILD.representative_sample(pd.DataFrame(rows), 20, "fixture")
        self.assertEqual(list(sample["annotation_unit_id"]), list(repeated["annotation_unit_id"]))

    def test_certification_separates_population_and_species_diagnostic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame([{
                "trait_id": "capitulum_orientation",
                "allow_multiple": "false",
            }]).to_csv(root / "ontology.csv", index=False)

            key_rows = []
            unit_rows = []
            annotation_rows = []
            task_index = 0
            for species_index in range(3):
                for within_species in range(4):
                    task_index += 1
                    task_id = f"task_{task_index:03d}"
                    unit_id = f"unit_{task_index:03d}"
                    stratum = "representative_population" if within_species < 2 else "species_diagnostic_topup"
                    key_rows.append({
                        "task_id": task_id,
                        "annotation_unit_id": unit_id,
                        "trait_id": "capitulum_orientation",
                        "taxon_name": f"Cirsium_{species_index}",
                        "double_label": "true" if within_species == 0 else "false",
                        "selection_stratum": stratum,
                        "ai_candidate_state": "upright",
                    })
                    unit_rows.append({
                        "annotation_unit_id": unit_id,
                        "trait_id": "capitulum_orientation",
                        "taxon_name": f"Cirsium_{species_index}",
                        "analysis_state_ai_conservative": "upright",
                    })
                    human_state = (
                        "nodding"
                        if species_index == 0 and within_species == 0
                        else "upright"
                    )
                    annotation_rows.append({
                        "task_id": task_id,
                        "annotation_unit_id": unit_id,
                        "trait_id": "capitulum_orientation",
                        "human_state": human_state,
                        "annotator_id": "primary",
                    })
                    if within_species == 0:
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
                "--min-assessable-per-trait", "5",
                "--accuracy-wilson-floor", "0",
                "--min-species-per-trait", "3",
                "--min-records-per-species", "4",
                "--worst-species-accuracy-floor", "0",
                "--human-human-wilson-floor", "0",
            ]
            with patch.object(sys, "argv", argv):
                EVAL.main()

            gate = pd.read_csv(output / "representative_holdout_trait_gate.csv")
            row = gate.iloc[0]
            self.assertEqual(int(row["n_population_scoreable"]), 6)
            self.assertEqual(int(row["n_species_diagnostic_scoreable"]), 12)
            self.assertAlmostEqual(float(row["population_accuracy"]), 5 / 6)
            self.assertEqual(int(row["n_species_eligible"]), 3)
            self.assertEqual(row["gate_decision"], "accept_for_analysis")


if __name__ == "__main__":
    unittest.main()

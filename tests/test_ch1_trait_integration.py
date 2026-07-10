from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"
ONTOLOGY = V2 / "ontology" / "ch1_trait_ontology.csv"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


INT = load_module("49_run_trait_integration_and_constraint.py", "ch1_trait_integration")


def write_long(path: Path, integrated: bool, n_species: int = 40, obs_per: int = 4) -> None:
    """Two ontology traits; when integrated, colour state is a deterministic function
    of the armature state (perfect coordination -> forbidden combinations)."""
    rng = np.random.default_rng(0)
    colours = ["white_or_cream", "purple_blue"]
    armature = ["unarmed", "long_spined"]
    rows = []
    for s in range(n_species):
        arm = armature[s % 2]
        if integrated:
            col = colours[armature.index(arm)]  # colour locked to armature
        else:
            col = colours[rng.integers(0, 2)]
        for o in range(obs_per):
            base = {"obs_id": f"s{s}_o{o}", "taxon_name": f"Cirsium_{s}"}
            rows.append({**base, "trait_id": "corolla_colour_class",
                         "observation_ai_conservative_state": col, "observation_ai_all_state": col})
            rows.append({**base, "trait_id": "bract_tip_armature",
                         "observation_ai_conservative_state": arm, "observation_ai_all_state": arm})
    pd.DataFrame(rows).to_csv(path, index=False)


class TestTraitIntegration(unittest.TestCase):
    def _run(self, integrated: bool):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            long = root / "long.csv"
            write_long(long, integrated=integrated)
            out = root / "out"
            with patch.object(sys, "argv", [
                "int", "--observation-long", str(long), "--ontology", str(ONTOLOGY),
                "--out-dir", str(out), "--min-species", "10", "--min-obs-per-species-trait", "2",
                "--n-permutations", "199",
            ]):
                INT.main()
            return json.loads((out / "trait_integration_constraint_report.json").read_text()), \
                pd.read_csv(out / "trait_combination_occupancy.csv")

    def test_integrated_traits_are_detected(self) -> None:
        report, occupancy = self._run(integrated=True)
        self.assertGreater(report["integration"]["mean_pairwise_cramers_v"], 0.8)
        self.assertLessEqual(report["integration"]["permutation_p_one_sided_high"], 0.05)
        # perfect coordination -> some combinations never realised
        self.assertGreaterEqual(int(occupancy["forbidden"].sum()), 1)

    def test_independent_traits_show_weak_integration(self) -> None:
        report, _ = self._run(integrated=False)
        self.assertLess(report["integration"]["mean_pairwise_cramers_v"], 0.5)


if __name__ == "__main__":
    unittest.main()

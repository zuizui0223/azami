from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT / "ch1_global" / "v2"
ONTOLOGY = V2 / "ontology" / "ch1_trait_ontology.csv"
PROMPTS = V2 / "ontology" / "clip_trait_prompt_spec.json"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ENSEMBLE = load_module("35_score_traits_ai_ensemble.py", "ch1_ai_ensemble")


class TestAutomatedEnsembleMeasurement(unittest.TestCase):
    def test_prompt_spec_matches_image_traits_and_retains_mucilage_exclusion(self) -> None:
        ontology = pd.read_csv(ONTOLOGY, dtype=str, keep_default_na=False)
        prompts = json.loads(PROMPTS.read_text(encoding="utf-8"))
        traits = ENSEMBLE.load_trait_specs(ontology, prompts)
        by_trait = {trait["trait_id"]: trait for trait in traits}
        self.assertEqual(by_trait["capitulum_orientation"]["mode"], "score")
        self.assertEqual(by_trait["capitulum_orientation"]["view"], "context")
        self.assertEqual(by_trait["corolla_colour_class"]["mode"], "score")
        self.assertEqual(by_trait["external_mucilage_visible"]["mode"], "not_scored")

    def test_ensemble_statuses_distinguish_strict_consensus_and_disagreement(self) -> None:
        strict_rows = pd.DataFrame([
            {"annotation_unit_id": "u1", "trait_id": "capitulum_orientation", "model_id": "m1", "variant_candidate_state": "upright", "variant_top_similarity": 0.8, "variant_similarity_margin": 0.2},
            {"annotation_unit_id": "u1", "trait_id": "capitulum_orientation", "model_id": "m1", "variant_candidate_state": "upright", "variant_top_similarity": 0.7, "variant_similarity_margin": 0.1},
            {"annotation_unit_id": "u1", "trait_id": "capitulum_orientation", "model_id": "m2", "variant_candidate_state": "upright", "variant_top_similarity": 0.6, "variant_similarity_margin": 0.3},
            {"annotation_unit_id": "u1", "trait_id": "capitulum_orientation", "model_id": "m2", "variant_candidate_state": "upright", "variant_top_similarity": 0.65, "variant_similarity_margin": 0.25},
        ])
        strict = ENSEMBLE.aggregate_group(strict_rows)
        self.assertEqual(strict["ai_candidate_state"], "upright")
        self.assertEqual(strict["n_variant_votes"], 4)
        self.assertEqual(strict["vote_fraction"], 1.0)
        self.assertEqual(
            ENSEMBLE.status_from_consensus(strict["vote_fraction"], 0.7, 0.5, 0.75),
            "strict_ensemble_consensus_higher_margin",
        )

        disagreement_rows = strict_rows.copy()
        disagreement_rows.loc[2:, "variant_candidate_state"] = "nodding"
        disagreement = ENSEMBLE.aggregate_group(disagreement_rows)
        self.assertEqual(disagreement["vote_fraction"], 0.5)
        self.assertEqual(
            ENSEMBLE.status_from_consensus(disagreement["vote_fraction"], 0.9, 0.5, 0.75),
            "model_variant_disagreement",
        )


if __name__ == "__main__":
    unittest.main()

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
ONTOLOGY = V2 / "ontology" / "ch1_trait_ontology.csv"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILD = load_module("47_build_representative_trait_holdout.py", "ch1_build_representative_holdout")
EVAL = load_module("48_evaluate_representative_trait_holdout.py", "ch1_evaluate_representative_holdout")


def make_candidates() -> pd.DataFrame:
    rows = []
    # 4 taxa x 3 spatial blocks x 4 units = 48 candidate units, single trait
    for taxon in ("Cirsium_a", "Cirsium_b", "Cirsium_c", "Cirsium_d"):
        for block in ("b1", "b2", "b3"):
            for unit in range(4):
                rows.append({
                    "annotation_unit_id": f"{taxon}_{block}_{unit}",
                    "trait_id": "corolla_colour_class",
                    "taxon_name": taxon,
                    "spatial_block_10deg": block,
                    "ai_candidate_state": "purple_blue",
                })
    return pd.DataFrame(rows)


class TestRepresentativeHoldout(unittest.TestCase):
    def test_build_is_stratified_and_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cand = root / "candidates.csv"
            make_candidates().to_csv(cand, index=False)
            out = root / "build"
            argv = ["build", "--candidate-units", str(cand), "--out-dir", str(out),
                    "--batch-name", "test", "--n-tasks-per-trait", "24", "--seed", "s1"]
            with patch.object(sys, "argv", argv):
                BUILD.main()
            key = pd.read_csv(out / "representative_trait_holdout_key" / "representative_trait_holdout_key.csv")
            public = pd.read_csv(out / "representative_trait_holdout_tasks" / "representative_trait_holdout_tasks.csv")
            self.assertEqual(len(key), 24)
            # public packet is blinded: no taxon / AI candidate leaked
            self.assertNotIn("taxon_name", public.columns)
            self.assertNotIn("ai_candidate_state", public.columns)
            # stratification spread all four taxa and all three blocks
            self.assertEqual(key["taxon_name"].nunique(), 4)
            self.assertEqual(key["spatial_block_10deg"].nunique(), 3)
            # determinism
            out2 = root / "build2"
            with patch.object(sys, "argv", ["build", "--candidate-units", str(cand), "--out-dir", str(out2),
                                            "--batch-name", "test", "--n-tasks-per-trait", "24", "--seed", "s1"]):
                BUILD.main()
            key2 = pd.read_csv(out2 / "representative_trait_holdout_key" / "representative_trait_holdout_key.csv")
            self.assertEqual(list(key["annotation_unit_id"]), list(key2["annotation_unit_id"]))

    def test_app_ready_packet_copies_images_and_matches_app_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cand = root / "candidates.csv"
            make_candidates().to_csv(cand, index=False)
            # fake packet manifest + image files for every candidate unit
            packet_root = root / "packet_src"
            manifest_rows = []
            for unit in make_candidates()["annotation_unit_id"].unique():
                rels = {}
                for col in ("source_image", "crop_path", "context_crop_path"):
                    rel = f"{col}/{unit}.jpg"
                    (packet_root / rel).parent.mkdir(parents=True, exist_ok=True)
                    (packet_root / rel).write_bytes(b"fake-image-bytes")
                    rels[col] = rel
                manifest_rows.append({"annotation_unit_id": unit, **rels})
            manifest = root / "packet_manifest.csv"
            pd.DataFrame(manifest_rows).to_csv(manifest, index=False)
            out = root / "build"
            with patch.object(sys, "argv", [
                "build", "--candidate-units", str(cand), "--out-dir", str(out), "--batch-name", "t",
                "--n-tasks-per-trait", "12", "--seed", "s1",
                "--packet-manifest", str(manifest), "--packet-root", str(packet_root),
            ]):
                BUILD.main()
            public = pd.read_csv(out / "representative_trait_holdout_tasks" / "representative_trait_holdout_tasks.csv")
            # app (31) REQUIRED columns present, and still blinded (no taxon/AI leak)
            for col in ("task_id", "annotation_unit_id", "trait_id", "source_image", "crop_path", "context_crop_path"):
                self.assertIn(col, public.columns)
            self.assertNotIn("taxon_name", public.columns)
            self.assertNotIn("ai_candidate_state", public.columns)
            # images actually copied into the packet
            first = public.iloc[0]
            self.assertTrue((out / "representative_trait_holdout_tasks" / first["crop_path"]).is_file())

    def _run_eval(self, root: Path, ai_correct: bool, worst_species_wrong: bool):
        cand = root / "candidates.csv"
        make_candidates().to_csv(cand, index=False)
        build = root / "build"
        with patch.object(sys, "argv", ["build", "--candidate-units", str(cand), "--out-dir", str(build),
                                         "--batch-name", "t", "--n-tasks-per-trait", "48",
                                         "--double-label-fraction", "0.5", "--seed", "s1"]):
            BUILD.main()
        key = pd.read_csv(build / "representative_trait_holdout_key" / "representative_trait_holdout_key.csv", dtype=str)

        # Build AI analysis units and matching human labels
        unit_rows, human_rows = [], []
        for _, row in key.iterrows():
            unit = row["annotation_unit_id"]
            taxon = row["taxon_name"]
            human_state = "purple_blue"
            ai = "purple_blue"
            if not ai_correct:
                ai = "pink"
            if worst_species_wrong and taxon == "Cirsium_d":
                ai = "pink"  # AI wrong only for one species -> species generalization gap
            unit_rows.append({
                "annotation_unit_id": unit, "trait_id": "corolla_colour_class", "taxon_name": taxon,
                "observation_ai_conservative_state": ai, "observation_ai_all_state": ai,
            })
            human_rows.append({
                "task_id": row["task_id"], "annotation_unit_id": unit, "trait_id": "corolla_colour_class",
                "human_state": human_state, "annotator_id": "annotator_1",
            })
            if row["double_label"] == "true":
                human_rows.append({
                    "task_id": row["task_id"], "annotation_unit_id": unit, "trait_id": "corolla_colour_class",
                    "human_state": human_state, "annotator_id": "annotator_2",
                })
        units_csv = root / "units.csv"
        ann_csv = root / "annotations.csv"
        pd.DataFrame(unit_rows).to_csv(units_csv, index=False)
        pd.DataFrame(human_rows).to_csv(ann_csv, index=False)
        eval_out = root / "eval"
        with patch.object(sys, "argv", [
            "eval", "--canonical-annotations", str(ann_csv), "--holdout-key", str(build / "representative_trait_holdout_key" / "representative_trait_holdout_key.csv"),
            "--analysis-units", str(units_csv), "--ontology", str(ONTOLOGY), "--out-dir", str(eval_out),
            "--min-assessable-per-trait", "20", "--min-species-per-trait", "3", "--min-records-per-species", "2",
        ]):
            EVAL.main()
        return json.loads((eval_out / "representative_holdout_certification.json").read_text())

    def test_gate_accepts_accurate_trait(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = self._run_eval(Path(tmp), ai_correct=True, worst_species_wrong=False)
            self.assertIn("corolla_colour_class", report["traits_accepted_for_analysis"])

    def test_gate_flags_species_generalization_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = self._run_eval(Path(tmp), ai_correct=True, worst_species_wrong=True)
            self.assertNotIn("corolla_colour_class", report["traits_accepted_for_analysis"])


if __name__ == "__main__":
    unittest.main()

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
ONTOLOGY = V2 / "ontology" / "ch1_trait_ontology.csv"


def load_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, V2 / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILD_AUDIT = load_module("30_build_blinded_trait_audit_packet.py", "ch1_build_trait_audit")
COMPILE_AUDIT = load_module("32_compile_trait_audit_responses.py", "ch1_compile_trait_audit")
EVALUATE_AUDIT = load_module("33_evaluate_trait_audit_against_clip.py", "ch1_evaluate_trait_audit")


class TestTraitAuditLoop(unittest.TestCase):
    def _make_packet(self, root: Path) -> tuple[Path, Path, Path]:
        packet_root = root / "full_packet"
        for dirname in ("source_images", "head_crops", "context_crops"):
            (packet_root / dirname).mkdir(parents=True, exist_ok=True)
        rows = []
        for index in range(1, 5):
            image = Image.new("RGB", (40, 30), color=(index * 20, index * 20, index * 20))
            source = f"source_images/source_{index}.jpg"
            head = f"head_crops/head_{index}.jpg"
            context = f"context_crops/context_{index}.jpg"
            image.save(packet_root / source)
            image.save(packet_root / head)
            image.save(packet_root / context)
            rows.append({
                "annotation_unit_id": f"unit_{index:02d}",
                "annotation_batch": "fixture",
                "source_image": source,
                "crop_path": head,
                "context_crop_path": context,
            })
        manifest = packet_root / "blinded_trait_annotation_manifest.csv"
        pd.DataFrame(rows).to_csv(manifest, index=False)

        proposals = root / "proposals_long.csv"
        proposal_rows = []
        for trait_id, states in {
            "capitulum_orientation": ("upright", "nodding"),
            "corolla_colour_class": ("pink", "white_or_cream"),
        }.items():
            for index in range(1, 5):
                proposal_rows.append({
                    "annotation_unit_id": f"unit_{index:02d}",
                    "trait_id": trait_id,
                    "proposal_status": "zero_shot_candidate_not_validated",
                    "ai_candidate_state": states[0],
                    "runner_up_state": states[1],
                    "similarity_margin": str(0.01 * index),
                    "margin_percentile_within_trait": str(index / 4),
                    "model_id": "fixture-clip",
                    "prompt_spec_version": "fixture-v1",
                })
        pd.DataFrame(proposal_rows).to_csv(proposals, index=False)

        prompt_spec = root / "prompt_spec.json"
        prompt_spec.write_text(json.dumps({
            "traits": {
                "capitulum_orientation": {"review_weight": 1.0},
                "corolla_colour_class": {"review_weight": 1.0},
            }
        }), encoding="utf-8")
        return packet_root, manifest, proposals, prompt_spec

    def test_blinded_audit_selection_then_compile_and_evaluate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            packet_root, manifest, proposals, prompt_spec = self._make_packet(root)
            audit_out = root / "audit_build"
            with patch.object(sys, "argv", [
                "build", "--packet-manifest", str(manifest), "--packet-root", str(packet_root),
                "--proposals-long", str(proposals), "--prompt-spec", str(prompt_spec),
                "--out-dir", str(audit_out), "--batch-name", "fixture", "--n-uncertainty-tasks", "2",
                "--n-calibration-tasks", "2", "--double-label-fraction", "1.0", "--seed", "7",
            ]):
                BUILD_AUDIT.main()

            public_tasks = pd.read_csv(
                audit_out / "blinded_trait_audit_packet" / "blinded_trait_audit_tasks.csv",
                dtype=str,
                keep_default_na=False,
            )
            private_key = pd.read_csv(
                audit_out / "private_audit_key" / "trait_audit_selection_key.csv",
                dtype=str,
                keep_default_na=False,
            )
            self.assertEqual(len(public_tasks), 4)
            self.assertEqual(len(private_key), 4)
            self.assertTrue(public_tasks["task_id"].is_unique)
            self.assertEqual(
                {"ai_candidate_state", "runner_up_state", "similarity_margin", "selection_stratum", "double_label"}.intersection(public_tasks.columns),
                set(),
            )
            self.assertTrue(private_key["double_label"].eq("true").all())
            for field in ("source_image", "crop_path", "context_crop_path"):
                self.assertTrue(all((audit_out / "blinded_trait_audit_packet" / value).is_file() for value in public_tasks[field]))

            responses_a = root / "responses_a.csv"
            responses_b = root / "responses_b.csv"
            response_rows = []
            for row in public_tasks.to_dict("records"):
                state = "upright" if row["trait_id"] == "capitulum_orientation" else "pink"
                response_rows.append({"task_id": row["task_id"], "human_state": state, "task_complete": "yes", "notes": "fixture"})
            pd.DataFrame(response_rows).to_csv(responses_a, index=False)
            pd.DataFrame(response_rows).to_csv(responses_b, index=False)

            compiled = root / "compiled"
            with patch.object(sys, "argv", [
                "compile", "--tasks", str(audit_out / "blinded_trait_audit_packet" / "blinded_trait_audit_tasks.csv"),
                "--ontology", str(ONTOLOGY), "--response", f"annotator_a={responses_a}",
                "--response", f"annotator_b={responses_b}", "--out-dir", str(compiled),
            ]):
                COMPILE_AUDIT.main()
            canonical = pd.read_csv(compiled / "trait_audit_annotations_canonical.csv", dtype=str, keep_default_na=False)
            self.assertEqual(len(canonical), 8)
            self.assertEqual(canonical["task_id"].nunique(), 4)
            self.assertEqual(canonical["annotator_id"].nunique(), 2)

            evaluation = root / "evaluation"
            with patch.object(sys, "argv", [
                "evaluate", "--canonical-annotations", str(compiled / "trait_audit_annotations_canonical.csv"),
                "--private-audit-key", str(audit_out / "private_audit_key" / "trait_audit_selection_key.csv"),
                "--ontology", str(ONTOLOGY), "--out-dir", str(evaluation),
            ]):
                EVALUATE_AUDIT.main()
            policy = pd.read_csv(evaluation / "trait_candidate_use_policy.csv", dtype=str, keep_default_na=False)
            agreement = pd.read_csv(evaluation / "human_human_agreement_by_trait.csv", dtype=str, keep_default_na=False)
            self.assertTrue(policy["hard_limit"].str.contains("No automatic trait acceptance", regex=False).all())
            self.assertTrue((pd.to_numeric(agreement["n_pairwise_records_designated_double_label"]) > 0).all())
            report = json.loads((evaluation / "trait_audit_evaluation_report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["n_pairwise_human_records_designated_double_label"], 4)


if __name__ == "__main__":
    unittest.main()

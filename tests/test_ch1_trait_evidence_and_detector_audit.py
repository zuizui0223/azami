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


VALIDATE_ONTOLOGY = load_module("12_validate_trait_ontology.py", "ch1_validate_ontology")
PACKETS = load_module("13_build_literature_extraction_packets.py", "ch1_literature_packets")
VALIDATE_EVIDENCE = load_module("14_validate_literature_evidence.py", "ch1_validate_evidence")
COMPILE_EVIDENCE = load_module("15_compile_literature_trait_summary.py", "ch1_compile_evidence")
AUDIT_MANIFEST = load_module("16_build_detector_audit_manifest.py", "ch1_detector_audit_manifest")
EVALUATE_AUDIT = load_module("17_evaluate_detector_audit.py", "ch1_detector_audit_evaluate")


class TestTraitEvidenceAndDetectorAudit(unittest.TestCase):
    def test_ontology_packets_evidence_and_conflict_preservation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ontology_out = root / "ontology_check"
            with patch.object(sys, "argv", ["validate", "--ontology", str(ONTOLOGY), "--out-dir", str(ontology_out)]):
                VALIDATE_ONTOLOGY.main()
            report = json.loads((ontology_out / "trait_ontology_validation.json").read_text(encoding="utf-8"))
            self.assertGreaterEqual(report["n_traits"], 5)

            passages = root / "passages.csv"
            pd.DataFrame([{
                "passage_id": "p1", "source_id": "s1", "accepted_taxon": "Cirsium alpha",
                "source_taxon_name": "C. alpha", "page_or_figure": "p. 1", "language": "en",
                "passage_text": "Heads nodding; corollas purple.", "rights_checked": "true",
            }]).to_csv(passages, index=False)
            packet_out = root / "packets"
            with patch.object(sys, "argv", ["packets", "--ontology", str(ONTOLOGY), "--passages", str(passages), "--out-dir", str(packet_out)]):
                PACKETS.main()
            task = json.loads((packet_out / "literature_extraction_tasks.jsonl").read_text(encoding="utf-8").splitlines()[0])
            self.assertEqual(task["passage_id"], "p1")
            self.assertIn("exact contiguous quote", task["prompt"])

            evidence = root / "evidence.csv"
            columns = pd.read_csv(packet_out / "literature_evidence_template.csv").columns
            rows = [
                {
                    "evidence_id": "e1", "passage_id": "p1", "source_id": "s1", "accepted_taxon": "Cirsium alpha",
                    "source_taxon_name": "C. alpha", "trait_id": "capitulum_orientation", "trait_state": "nodding",
                    "scope": "typical heads", "life_stage": "anthesis", "geographic_scope": "", "evidence_quote": "Heads nodding",
                    "evidence_translation": "", "llm_confidence": "high", "extraction_method": "llm_assisted",
                    "review_status": "accepted", "reviewer_id": "reviewer_1", "notes": "",
                },
                {
                    "evidence_id": "e2", "passage_id": "p1", "source_id": "s1", "accepted_taxon": "Cirsium alpha",
                    "source_taxon_name": "C. alpha", "trait_id": "capitulum_orientation", "trait_state": "upright",
                    "scope": "", "life_stage": "", "geographic_scope": "", "evidence_quote": "Heads nodding",
                    "evidence_translation": "", "llm_confidence": "low", "extraction_method": "llm_assisted",
                    "review_status": "accepted", "reviewer_id": "reviewer_1", "notes": "deliberate conflict fixture",
                },
            ]
            pd.DataFrame(rows, columns=columns).to_csv(evidence, index=False)
            evidence_out = root / "validated"
            with patch.object(sys, "argv", ["evidence", "--ontology", str(ONTOLOGY), "--passages", str(passages), "--evidence", str(evidence), "--out-dir", str(evidence_out)]):
                VALIDATE_EVIDENCE.main()
            valid = pd.read_csv(evidence_out / "literature_evidence_valid.csv", dtype=str, keep_default_na=False)
            self.assertEqual(len(valid), 2)
            summary_out = root / "summary"
            with patch.object(sys, "argv", ["summary", "--evidence", str(evidence_out / "literature_evidence_valid.csv"), "--out-dir", str(summary_out)]):
                COMPILE_EVIDENCE.main()
            summary = pd.read_csv(summary_out / "literature_trait_summary.csv", dtype=str, keep_default_na=False)
            self.assertEqual(summary.loc[0, "interpretation"], "conflicting_reviewed_evidence")
            self.assertEqual(summary.loc[0, "trait_state_consensus"], "")

    def test_detector_manifest_and_evaluation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            metadata = root / "photo_metadata.csv"
            rows = []
            for i, species in enumerate(["Cirsium alpha", "Cirsium beta", "Cirsium gamma", "Cirsium delta"], start=1):
                rows.append({
                    "obs_id": str(i), "photo_id": str(10 + i), "taxon_name": species, "taxon_rank": "species",
                    "source_file_stem": f"inat_photo_{10+i}", "captive": "false", "latitude": str(30 + i), "longitude": str(130 + i),
                    "small_image_url": f"https://example.org/{i}/small.jpg", "medium_image_url": f"https://example.org/{i}/medium.jpg",
                })
            pd.DataFrame(rows).to_csv(metadata, index=False)
            audit_out = root / "audit"
            with patch.object(sys, "argv", [
                "audit", "--metadata", str(metadata), "--out-dir", str(audit_out), "--n-images", "4",
                "--min-images-per-species", "1", "--max-images-per-species", "2", "--double-label-fraction", "0.5", "--seed", "7",
            ]):
                AUDIT_MANIFEST.main()
            queue = pd.read_csv(audit_out / "detector_audit_queue.csv", dtype=str, keep_default_na=False)
            self.assertEqual(len(queue), 4)
            self.assertTrue(queue["screen_download_filename"].str.contains("Cirsium", regex=False).eq(False).all())
            first_id, second_id = queue.loc[0, "queue_id"], queue.loc[1, "queue_id"]
            annotations = root / "annotations.csv"
            pd.DataFrame([
                {"audit_id": first_id, "annotation_status": "adjudicated", "assessability": "assessable", "gt_index": "1", "gt_label": "visible_capitulum", "gt_x1": "10", "gt_y1": "10", "gt_x2": "30", "gt_y2": "30"},
                {"audit_id": second_id, "annotation_status": "adjudicated", "assessability": "assessable", "gt_index": "0", "gt_label": "none", "gt_x1": "", "gt_y1": "", "gt_x2": "", "gt_y2": ""},
            ]).to_csv(annotations, index=False)
            predictions = root / "predictions.csv"
            pd.DataFrame([
                {"queue_id": first_id, "yolo_conf": "0.9", "bbox_x1": "10", "bbox_y1": "10", "bbox_x2": "30", "bbox_y2": "30"},
                {"queue_id": second_id, "yolo_conf": "0.9", "bbox_x1": "40", "bbox_y1": "40", "bbox_x2": "50", "bbox_y2": "50"},
            ]).to_csv(predictions, index=False)
            eval_out = root / "evaluation"
            with patch.object(sys, "argv", ["evaluate", "--annotations", str(annotations), "--predictions", str(predictions), "--out-dir", str(eval_out)]):
                EVALUATE_AUDIT.main()
            metrics = json.loads((eval_out / "detector_audit_metrics.json").read_text(encoding="utf-8"))
            self.assertEqual(metrics["box_true_positive"], 1)
            self.assertEqual(metrics["image_presence_false_positive"], 1)


if __name__ == "__main__":
    unittest.main()

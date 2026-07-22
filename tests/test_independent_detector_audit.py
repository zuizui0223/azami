from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class DetectionAuditUtilitiesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = load_module(
            "detection_audit", ROOT / "src" / "azami_ch1" / "detection_audit.py"
        )

    def test_greedy_matching_and_metrics(self):
        truth = [(0, 0, 10, 10), (20, 20, 30, 30)]
        predictions = [(1, 1, 9, 9), (40, 40, 50, 50)]
        matches, missing_truth, extra_predictions = self.audit.greedy_match(
            truth, predictions, 0.5
        )
        self.assertEqual(len(matches), 1)
        self.assertEqual(missing_truth, [1])
        self.assertEqual(extra_predictions, [1])
        precision, recall, f1 = self.audit.precision_recall_f1(1, 1, 1)
        self.assertAlmostEqual(precision, 0.5)
        self.assertAlmostEqual(recall, 0.5)
        self.assertAlmostEqual(f1, 0.5)

    def test_wilson_and_area_classification(self):
        low, high = self.audit.wilson_interval(8, 10)
        self.assertLess(low, 0.8)
        self.assertGreater(high, 0.8)
        self.assertEqual(
            self.audit.classify_relative_area((0, 0, 5, 5), 100, 100),
            "small_lt_1pct",
        )
        self.assertEqual(
            self.audit.classify_relative_area((0, 0, 30, 30), 100, 100),
            "large_ge_5pct",
        )


class IndependentAuditPacketTest(unittest.TestCase):
    def test_builder_excludes_prior_photo_and_observation_and_blinds_packet(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = pd.DataFrame(
                [
                    {"obs_id": "o1", "photo_id": "p1", "taxon_name": "Cirsium a", "taxon_rank": "species", "source_file_stem": "p1", "captive": "false", "latitude": "40", "longitude": "10", "medium_image_url": "https://example.org/p1.jpg"},
                    {"obs_id": "o2", "photo_id": "p2", "taxon_name": "Cirsium a", "taxon_rank": "species", "source_file_stem": "p2", "captive": "false", "latitude": "41", "longitude": "11", "medium_image_url": "https://example.org/p2.jpg"},
                    {"obs_id": "o3", "photo_id": "p3", "taxon_name": "Cirsium b", "taxon_rank": "species", "source_file_stem": "p3", "captive": "false", "latitude": "20", "longitude": "20", "medium_image_url": "https://example.org/p3.jpg"},
                    {"obs_id": "o4", "photo_id": "p4", "taxon_name": "Cirsium b", "taxon_rank": "species", "source_file_stem": "p4", "captive": "false", "latitude": "21", "longitude": "21", "medium_image_url": "https://example.org/p4.jpg"},
                    {"obs_id": "o5", "photo_id": "p5", "taxon_name": "Cirsium c", "taxon_rank": "species", "source_file_stem": "p5", "captive": "false", "latitude": "-10", "longitude": "30", "medium_image_url": "https://example.org/p5.jpg"},
                    {"obs_id": "o6", "photo_id": "p6", "taxon_name": "Cirsium c", "taxon_rank": "species", "source_file_stem": "p6", "captive": "false", "latitude": "-11", "longitude": "31", "medium_image_url": "https://example.org/p6.jpg"},
                ]
            )
            metadata_path = root / "metadata.csv"
            metadata.to_csv(metadata_path, index=False)
            exclusion_path = root / "exclude.csv"
            pd.DataFrame([{"obs_id": "o1", "photo_id": "p1"}]).to_csv(
                exclusion_path, index=False
            )
            out = root / "out"
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "analysis" / "build_independent_detector_audit_packet.py"),
                    "--metadata", str(metadata_path),
                    "--exclude-queue", str(exclusion_path),
                    "--out-dir", str(out),
                    "--n-images", "3",
                    "--target-images-per-species", "1",
                    "--max-images-per-species", "2",
                    "--double-label-fraction", "0.34",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            queue = pd.read_csv(out / "detector_independent_audit_queue.csv", dtype=str)
            self.assertNotIn("p1", set(queue["photo_id"]))
            self.assertNotIn("o1", set(queue["obs_id"]))
            self.assertEqual(len(queue), 3)
            blinded = pd.read_csv(
                out / "detector_independent_audit_blinded_manifest.csv", dtype=str
            )
            self.assertNotIn("taxon_name", blinded.columns)
            self.assertNotIn("latitude", blinded.columns)
            self.assertNotIn("selected_image_url", blinded.columns)
            report = (out / "detector_independent_audit_report.json").read_text()
            self.assertIn('"photo_overlap_with_exclusions": 0', report)
            self.assertIn('"observation_overlap_with_exclusions": 0', report)


class HumanAnnotationAgreementTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.finalizer = load_module(
            "finalize_detector_audit_annotations",
            ROOT / "analysis" / "finalize_detector_audit_annotations.py",
        )

    def annotation_frame(self, annotator: str, shift: float = 0) -> pd.DataFrame:
        return pd.DataFrame(
            [{
                "audit_id": "det_ind_00001", "annotator_id": annotator,
                "annotation_status": "complete", "assessability": "assessable",
                "image_quality": "good", "image_width": "100", "image_height": "100",
                "gt_index": "1", "gt_label": "visible_capitulum",
                "life_stage": "anthesis", "occlusion": "none",
                "edge_truncated": "false", "gt_x1": str(10 + shift),
                "gt_y1": str(10 + shift), "gt_x2": str(50 + shift),
                "gt_y2": str(50 + shift), "notes": "",
            }]
        )

    def test_double_annotations_agree_by_iou(self):
        first = self.finalizer.validate_group(self.annotation_frame("annotator_1"))
        second = self.finalizer.validate_group(self.annotation_frame("annotator_2", 2))
        comparison = self.finalizer.compare_double(first, second, 0.5)
        self.assertTrue(comparison["agreement"])
        self.assertEqual(comparison["matched"], 1)

    def test_double_annotations_require_adjudication_when_disjoint(self):
        first = self.finalizer.validate_group(self.annotation_frame("annotator_1"))
        second = self.finalizer.validate_group(self.annotation_frame("annotator_2", 45))
        comparison = self.finalizer.compare_double(first, second, 0.5)
        self.assertFalse(comparison["agreement"])


class IndependentDetectorEvaluationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evaluator = load_module(
            "evaluate_independent_detector_audit",
            ROOT / "analysis" / "evaluate_independent_detector_audit.py",
        )

    def test_threshold_evaluation_counts_tp_fp_fn(self):
        annotations = pd.DataFrame(
            [
                {"image_id": "a", "object_id": "a1", "x1": 0, "y1": 0, "x2": 10, "y2": 10},
                {"image_id": "b", "object_id": "b1", "x1": 20, "y1": 20, "x2": 30, "y2": 30},
            ]
        )
        annotations["assessable_bool"] = True
        annotations["image_width_num"] = 100.0
        annotations["image_height_num"] = 100.0
        annotations["annotation_state"] = "assessable"
        annotations["image_quality"] = "good"
        annotations["life_stage"] = "anthesis"
        annotations["occlusion"] = "none"
        annotations["edge_truncated"] = "false"
        predictions = pd.DataFrame(
            [
                {"queue_id": "a", "bbox_x1": 0, "bbox_y1": 0, "bbox_x2": 10, "bbox_y2": 10, "yolo_conf": 0.9},
                {"queue_id": "a", "bbox_x1": 40, "bbox_y1": 40, "bbox_x2": 50, "bbox_y2": 50, "yolo_conf": 0.8},
                {"queue_id": "b", "bbox_x1": 60, "bbox_y1": 60, "bbox_x2": 70, "bbox_y2": 70, "yolo_conf": 0.7},
            ]
        )
        predictions["confidence"] = predictions["yolo_conf"].astype(float)
        predictions["prediction_id"] = ["p1", "p2", "p3"]
        metrics, _, _ = self.evaluator.evaluate_threshold(
            annotations, predictions, threshold=0.25, iou_threshold=0.5
        )
        self.assertEqual(metrics["true_positive"], 1)
        self.assertEqual(metrics["false_positive"], 2)
        self.assertEqual(metrics["false_negative"], 1)
        self.assertAlmostEqual(metrics["precision"], 1 / 3)
        self.assertAlmostEqual(metrics["recall"], 1 / 2)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Validate and fingerprint a detector model package before Chapter 1 inference.

A valid package contains exactly one `.pt` weight file and one `model_card.json`.
The model card records the biological target separately from the weight's raw
class label, so legacy labels such as `flower` can be retained without silently
renaming the trained model's output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REQUIRED_CARD_FIELDS = {
    "model_id",
    "semantic_target",
    "detector_output_class",
    "training_data_id",
    "training_commit_or_run",
    "training_split_definition",
    "class_definition",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate a single-class YOLO detector package.")
    parser.add_argument("--package-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--expected-semantic-target", default="visible_capitulum")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_names(names: Any) -> dict[str, str]:
    if isinstance(names, dict):
        return {str(key): str(value) for key, value in names.items()}
    if isinstance(names, list):
        return {str(index): str(value) for index, value in enumerate(names)}
    raise ValueError(f"Unsupported Ultralytics model.names type: {type(names).__name__}")


def main() -> None:
    args = parse_args()
    package_dir = Path(args.package_dir)
    out_dir = Path(args.out_dir)
    if not package_dir.is_dir():
        raise FileNotFoundError(f"--package-dir is not a directory: {package_dir}")
    weights = sorted(package_dir.rglob("*.pt"))
    if len(weights) != 1:
        raise ValueError(f"Expected exactly one .pt weights file below {package_dir}; found {len(weights)}: {weights}")
    card_path = package_dir / "model_card.json"
    if not card_path.is_file():
        raise FileNotFoundError("Missing model_card.json at the root of the detector package")
    card = json.loads(card_path.read_text(encoding="utf-8"))
    missing = REQUIRED_CARD_FIELDS.difference(card)
    if missing:
        raise ValueError(f"model_card.json missing required fields: {sorted(missing)}")
    if card["semantic_target"] != args.expected_semantic_target:
        raise ValueError(
            f"model_card semantic_target={card['semantic_target']!r}; expected {args.expected_semantic_target!r}"
        )
    if not isinstance(card["class_definition"], str) or not card["class_definition"].strip():
        raise ValueError("model_card.class_definition must be a non-empty description")

    try:
        import ultralytics
        from ultralytics import YOLO
    except ImportError as error:
        raise SystemExit("Ultralytics is required to inspect detector weights") from error

    weight_path = weights[0]
    model = YOLO(str(weight_path))
    names = normalize_names(getattr(model, "names", {}))
    if len(names) != 1:
        raise ValueError(f"Chapter 1 detector must be single-class; weights expose {len(names)} classes: {names}")
    model_output_class = next(iter(names.values()))
    if str(card["detector_output_class"]).strip() != model_output_class:
        raise ValueError(
            "model_card.detector_output_class does not match the class name embedded in weights: "
            f"card={card['detector_output_class']!r}, weights={model_output_class!r}"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    report = {
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "package_dir": str(package_dir.resolve()),
        "weight_file": weight_path.name,
        "weight_sha256": sha256(weight_path),
        "weight_size_bytes": weight_path.stat().st_size,
        "ultralytics_version": getattr(ultralytics, "__version__", "unknown"),
        "model_names": names,
        "model_card": card,
        "semantic_target": args.expected_semantic_target,
        "inference_policy": "Detector outputs remain predictions until evaluated against independent adjudicated source-image annotations.",
    }
    (out_dir / "detector_model_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

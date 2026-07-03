#!/usr/bin/env python3
"""Train a provisional visible-capitulum YOLO detector and export stable artifacts.

The trainer's own save directory may be controlled by the installed Ultralytics
runtime. This wrapper copies the resulting weights and selected diagnostics into
a predictable project directory, records SHA256 fingerprints, and makes clear
that pseudo-label validation is not independent detector accuracy.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train and export a bootstrap visible-capitulum detector.")
    parser.add_argument("--dataset-yaml", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument("--imgsz", type=int, default=640)
    parser.add_argument("--batch", type=int, default=8)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--model", default="yolo11n.pt")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def candidate_run_dirs(result: Any, model: Any) -> list[Path]:
    values = [
        getattr(result, "save_dir", None),
        getattr(getattr(model, "trainer", None), "save_dir", None),
        getattr(model, "save_dir", None),
    ]
    seen: set[str] = set()
    output: list[Path] = []
    for value in values:
        if value is None:
            continue
        path = Path(value)
        key = str(path)
        if key not in seen:
            seen.add(key)
            output.append(path)
    return output


def select_run_dir(result: Any, model: Any) -> Path:
    candidates = candidate_run_dirs(result, model)
    for directory in candidates:
        if (directory / "weights" / "best.pt").is_file():
            return directory
    raise FileNotFoundError(
        "Ultralytics training completed but no best.pt could be found in candidate run directories: "
        + ", ".join(str(path) for path in candidates)
    )


def read_last_metrics(results_csv: Path) -> dict[str, str]:
    if not results_csv.is_file():
        return {}
    with results_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows[-1] if rows else {}


def main() -> None:
    args = parse_args()
    if args.epochs < 1 or args.imgsz < 32 or args.batch < 1 or args.workers < 0:
        raise ValueError("Invalid training parameters")
    dataset_yaml = Path(args.dataset_yaml)
    if not dataset_yaml.is_file():
        raise FileNotFoundError(dataset_yaml)
    try:
        import ultralytics
        from ultralytics import YOLO
    except ImportError as error:
        raise SystemExit("Ultralytics is required") from error

    model = YOLO(args.model)
    result = model.train(
        data=str(dataset_yaml),
        epochs=args.epochs,
        imgsz=args.imgsz,
        batch=args.batch,
        device=args.device,
        workers=args.workers,
        project="bootstrap/runtime_runs",
        name="visible_capitulum_bootstrap",
        exist_ok=True,
        pretrained=True,
        verbose=False,
    )
    run_dir = select_run_dir(result, model)
    best = run_dir / "weights" / "best.pt"
    last = run_dir / "weights" / "last.pt"

    out_dir = Path(args.out_dir)
    weights_dir = out_dir / "weights"
    diagnostics_dir = out_dir / "diagnostics"
    weights_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    stable_best = weights_dir / "best.pt"
    stable_last = weights_dir / "last.pt"
    shutil.copy2(best, stable_best)
    if last.is_file():
        shutil.copy2(last, stable_last)
    for filename in ("results.csv", "args.yaml", "labels.jpg", "results.png", "PR_curve.png", "F1_curve.png", "P_curve.png", "R_curve.png"):
        source = run_dir / filename
        if source.is_file():
            shutil.copy2(source, diagnostics_dir / filename)

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "ultralytics_version": getattr(ultralytics, "__version__", "unknown"),
        "initialization": args.model,
        "dataset_yaml": str(dataset_yaml.resolve()),
        "actual_ultralytics_run_dir": str(run_dir.resolve()),
        "stable_best_weights": str(stable_best.resolve()),
        "best_weight_sha256": sha256(stable_best),
        "best_weight_size_bytes": stable_best.stat().st_size,
        "stable_last_weights": str(stable_last.resolve()) if stable_last.is_file() else "",
        "last_weight_sha256": sha256(stable_last) if stable_last.is_file() else "",
        "epochs": args.epochs,
        "imgsz": args.imgsz,
        "batch": args.batch,
        "device": args.device,
        "target": "visible_capitulum",
        "training_labels": "open-vocabulary pseudo-label proposals only",
        "metric_interpretation": "Training/validation metrics quantify agreement with pseudo-labels derived from the same open-vocabulary proposal process. They are not independent accuracy, precision, recall, or mAP against human annotations.",
        "status": "bootstrap proposal model; not validated and not for biological inference",
        "last_results_row": read_last_metrics(run_dir / "results.csv"),
    }
    (out_dir / "bootstrap_training_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

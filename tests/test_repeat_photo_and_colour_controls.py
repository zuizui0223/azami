from __future__ import annotations

import importlib.util
import tempfile
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_module(relative: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REPEAT = load_module(
    "analysis/build_repeat_photo_bias_control_cohort.py", "repeat_photo_control"
)
COLOUR = load_module(
    "analysis/run_colour_negative_control.py", "colour_negative_control"
)


def test_repeat_photo_selection_uses_only_photo_availability() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        strict = pd.DataFrame({
            "obs_id": ["1", "2", "3"],
            "taxon_name": ["Cirsium a", "Cirsium b", "Cirsium c"],
            "observed_on": ["2020-01-01", "2020-01-02", "2020-01-03"],
        })
        metadata = pd.DataFrame([
            {"obs_id": obs, "photo_id": photo, "photo_index": index,
             "photo_license_code": "cc-by", "medium_image_url": url,
             "large_image_url": "", "small_image_url": "", "raw_image_url": ""}
            for obs, photo, index, url in [
                ("1", "p1", "0", "https://example.org/p1.jpg"),
                ("1", "p2", "1", "https://example.org/p2.jpg"),
                ("2", "p3", "0", "https://example.org/p3.jpg"),
                ("3", "p4", "0", "https://example.org/p4.jpg"),
                ("3", "p5", "1", ""),
            ]
        ])
        metadata_path = root / "metadata.csv"
        metadata.to_csv(metadata_path, index=False)
        manifest, summary = REPEAT.build_repeat_cohort(strict, metadata_path, 2)
        assert set(manifest["obs_id"]) == {"1"}
        assert len(manifest) == 2
        assert int(summary.iloc[0]["n_additional_photos"]) == 1


def test_colour_controls_recover_outer_green_pixels() -> None:
    image = np.full((160, 160, 3), (40, 130, 40), dtype=np.uint8)
    cv2.circle(image, (80, 80), 42, (180, 60, 220), -1)
    result = COLOUR.measure_controls(image, minimum=200)
    assert result["status"] == "usable"
    assert result["green_status"] == "usable"
    assert int(result["green_n_pixels"]) >= 200
    assert float(result["green_lab_chroma"]) > 0


def test_paired_interaction_detects_background_difference() -> None:
    rng = np.random.default_rng(4)
    rows = []
    for taxon_index in range(12):
        taxon = f"Cirsium {taxon_index}"
        for observation in range(12):
            climate = observation - 5.5
            rows.append({
                "obs_id": f"{taxon_index}_{observation}",
                "taxon_name": taxon,
                "chelsa_bio12": climate,
                "corolla_lab_chroma_median": 30 + 1.8 * climate + rng.normal(0, 0.3),
                "border_lab_chroma_median": 20 + 0.1 * climate + rng.normal(0, 0.3),
            })
    result = COLOUR.paired_region_fit(
        pd.DataFrame(rows), "chroma", "border", minimum_taxon_n=5
    )
    assert result["status"] == "ok"
    assert float(result["interaction_p_value"]) < 0.001
    assert float(result["flower_slope_lab_units"]) > float(result["background_slope_lab_units"])

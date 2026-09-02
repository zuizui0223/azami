#!/usr/bin/env python3
"""Export blinded CSV packets for Chapter 1 image-level annotation.

The input manifest retains taxonomy, coordinates, and validation-group fields.
The exported packets intentionally do not: annotators should judge what is in
an image, not what species it belongs to or where it was photographed.

Example
-------
python ch1_global/v2/01_build_blinded_annotation_packets.py \
  --manifest <annotation_manifest.csv> \
  --out-dir <annotation_packets>
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_MANIFEST_COLUMNS = {
    "annotation_unit_id",
    "annotation_batch",
    "source_image",
    "crop_path",
    "needs_detector_box",
    "needs_keypoints",
    "needs_segmentation",
    "needs_trait_labels",
}

FORBIDDEN_PACKET_COLUMNS = {
    "species",
    "taxon_name",
    "scientific_name",
    "latitude",
    "longitude",
    "spatial_block",
    "species_cv_fold",
    "spatial_cv_fold",
    "observation_group",
    "photo_group",
    "detector_confidence",
    "raw_image_url",
    "geoprivacy",
    "obscured",
    "positional_accuracy",
}

TASK_COLUMNS = [
    "annotation_unit_id",
    "source_image",
    "crop_path",
    "needs_detector_box",
    "needs_keypoints",
    "needs_segmentation",
    "needs_trait_labels",
]

RESPONSE_COLUMNS = [
    "flowering_qc",
    "detector_box_qc",
    "camera_tilt_flag",
    "view_angle_class",
    "head_orientation_class",
    "flower_colour_class",
    "colour_qc",
    "involucral_cover_score_0_to_4",
    "head_shape_class",
    "corolla_display_score_0_to_4",
    "pedicel_x_px",
    "pedicel_y_px",
    "capitulum_apex_x_px",
    "capitulum_apex_y_px",
    "annotation_complete",
    "notes",
]

INSTRUCTIONS = """# Ch.1 image annotation packet

Do not infer any trait from the file name, remembered locality, or species identity. Score what is visible.

Allowed values
--------------
- flowering_qc: usable | not_flowering | too_small_or_blurred | occluded | wrong_taxon_or_non-Cirsium | uncertain
- detector_box_qc: accept | revise | reject | not_assessable
- camera_tilt_flag: vertical_reference_clear | vertical_reference_uncertain | tilted_or_rotated | not_assessable
- view_angle_class: lateral | oblique | frontal | not_assessable
- head_orientation_class: upright | inclined | nodding | not_assessable
- flower_colour_class: white | pale | dark | not_assessable
- colour_qc: usable | illumination_or_white_balance_unusable | not_assessable
- involucral_cover_score_0_to_4: 0, 1, 2, 3, 4, or not_assessable
- head_shape_class: globose | depressed | cylindrical | not_assessable
- corolla_display_score_0_to_4: 0, 1, 2, 3, 4, or not_assessable
- annotation_complete: yes | no

Orientation requires a reliable biological vertical reference (upright stem or horizon), not merely the image frame. Cover and display are primary only for lateral images. Leave pixel-coordinate fields blank unless the image can support the keypoint task.

Do not alter task columns. Add responses only in the blank response columns.
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build blinded Ch.1 annotation packets.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--seed", type=int, default=20260703)
    parser.add_argument("--primary-annotator", default="annotator_1")
    parser.add_argument("--secondary-annotator", default="annotator_2")
    parser.add_argument(
        "--allow-potential-filename-leak",
        action="store_true",
        help="Allow a source/crop path containing a taxon name; discouraged.",
    )
    return parser.parse_args()


def compact(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def potential_filename_leaks(df: pd.DataFrame) -> pd.DataFrame:
    """Conservative audit: paths should not visibly embed a species/taxon string."""
    if "species" not in df.columns:
        return df.iloc[0:0].copy()
    flagged = []
    for _, row in df.iterrows():
        species = compact(row.get("species", ""))
        if len(species) < 6:
            continue
        path_text = compact(f"{row.get('source_image', '')} {row.get('crop_path', '')}")
        if species and species in path_text:
            flagged.append(row)
    return pd.DataFrame(flagged, columns=df.columns)


def make_packet(df: pd.DataFrame, seed: int) -> pd.DataFrame:
    packet = df[TASK_COLUMNS].copy()
    for column in RESPONSE_COLUMNS:
        packet[column] = ""
    packet = packet.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    overlap = FORBIDDEN_PACKET_COLUMNS.intersection(packet.columns)
    if overlap:
        raise RuntimeError(f"Blinding failure: forbidden columns in packet: {sorted(overlap)}")
    return packet


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, low_memory=False)
    missing = REQUIRED_MANIFEST_COLUMNS.difference(manifest.columns)
    if missing:
        raise ValueError(f"Manifest is missing required columns: {sorted(missing)}")
    if manifest["annotation_unit_id"].duplicated().any():
        raise ValueError("Manifest has duplicate annotation_unit_id values")

    leaks = potential_filename_leaks(manifest)
    if not leaks.empty and not args.allow_potential_filename_leak:
        examples = leaks[["annotation_unit_id", "species", "source_image", "crop_path"]].head(10)
        raise ValueError(
            "Potential taxonomy leak in image path(s). Rename/relink image files before annotation, "
            "or explicitly use --allow-potential-filename-leak. Examples:\n"
            + examples.to_string(index=False)
        )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    calibration = manifest.loc[manifest["annotation_batch"].eq("calibration_double_label")].copy()
    if calibration.empty:
        raise ValueError("No calibration_double_label rows found in manifest")

    primary = make_packet(manifest, args.seed)
    secondary = make_packet(calibration, args.seed + 1)
    primary_path = out_dir / f"{args.primary_annotator}_packet.csv"
    secondary_path = out_dir / f"{args.secondary_annotator}_packet.csv"
    primary.to_csv(primary_path, index=False, encoding="utf-8-sig")
    secondary.to_csv(secondary_path, index=False, encoding="utf-8-sig")
    (out_dir / "README_annotation_instructions.md").write_text(INSTRUCTIONS, encoding="utf-8")

    qc = pd.DataFrame([
        {
            "packet": primary_path.name,
            "annotator": args.primary_annotator,
            "n_rows": len(primary),
            "n_calibration_rows": int(manifest["annotation_batch"].eq("calibration_double_label").sum()),
            "blinded_taxonomy": True,
            "blinded_coordinates": True,
            "potential_filename_leaks": len(leaks),
        },
        {
            "packet": secondary_path.name,
            "annotator": args.secondary_annotator,
            "n_rows": len(secondary),
            "n_calibration_rows": len(secondary),
            "blinded_taxonomy": True,
            "blinded_coordinates": True,
            "potential_filename_leaks": len(leaks),
        },
    ])
    qc.to_csv(out_dir / "packet_qc.csv", index=False, encoding="utf-8-sig")
    print(f"[OK] Wrote blinded packets to {out_dir}")
    print(qc.to_string(index=False))


if __name__ == "__main__":
    main()

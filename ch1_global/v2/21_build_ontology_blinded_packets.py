#!/usr/bin/env python3
"""Export taxonomy-blinded packets for the ontology-driven annotation app.

Unlike the legacy packet format, this packet preserves an optional
`context_crop_path` so capitulum orientation can be judged from the head plus
peduncle rather than from a tight head crop alone.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

REQUIRED = {"annotation_unit_id", "annotation_batch", "source_image", "crop_path"}
FORBIDDEN = {
    "species", "taxon_name", "scientific_name", "latitude", "longitude", "spatial_block",
    "species_cv_fold", "spatial_cv_fold", "observation_group", "photo_group", "detector_confidence",
    "raw_image_url", "geoprivacy", "obscured", "positional_accuracy",
}
TASK_COLUMNS = ["annotation_unit_id", "source_image", "crop_path", "context_crop_path", "needs_trait_labels"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build context-aware blinded ontology annotation packets.")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--seed", type=int, default=20260703)
    parser.add_argument("--primary-annotator", default="annotator_1")
    parser.add_argument("--secondary-annotator", default="annotator_2")
    parser.add_argument("--allow-potential-filename-leak", action="store_true")
    return parser.parse_args()


def compact(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def potential_leaks(frame: pd.DataFrame) -> pd.DataFrame:
    if "species" not in frame.columns:
        return frame.iloc[0:0].copy()
    flagged = []
    for _, row in frame.iterrows():
        species = compact(row.get("species", ""))
        if len(species) < 6:
            continue
        paths = compact(" ".join(str(row.get(column, "")) for column in ("source_image", "crop_path", "context_crop_path")))
        if species in paths:
            flagged.append(row)
    return pd.DataFrame(flagged, columns=frame.columns)


def make_packet(frame: pd.DataFrame, seed: int) -> pd.DataFrame:
    packet = frame.copy()
    for column in TASK_COLUMNS:
        if column not in packet.columns:
            packet[column] = "" if column != "needs_trait_labels" else True
    packet = packet[TASK_COLUMNS].sample(frac=1.0, random_state=seed).reset_index(drop=True)
    overlap = FORBIDDEN.intersection(packet.columns)
    if overlap:
        raise RuntimeError(f"Blinding failure: forbidden columns in packet: {sorted(overlap)}")
    return packet


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, dtype=str, keep_default_na=False)
    missing = REQUIRED.difference(manifest.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {sorted(missing)}")
    if manifest["annotation_unit_id"].duplicated().any():
        raise ValueError("Manifest has duplicate annotation_unit_id")
    leaks = potential_leaks(manifest)
    if not leaks.empty and not args.allow_potential_filename_leak:
        raise ValueError("Potential taxonomy leak in paths. Rename/relink files or use --allow-potential-filename-leak deliberately.")
    calibration = manifest.loc[manifest["annotation_batch"].eq("calibration_double_label")].copy()
    if calibration.empty:
        raise ValueError("No calibration_double_label rows found")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    primary = make_packet(manifest, args.seed)
    secondary = make_packet(calibration, args.seed + 1)
    primary.to_csv(out_dir / f"{args.primary_annotator}_ontology_packet.csv", index=False, encoding="utf-8-sig")
    secondary.to_csv(out_dir / f"{args.secondary_annotator}_ontology_packet.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([{
        "primary_rows": len(primary), "secondary_rows": len(secondary),
        "potential_filename_leaks": len(leaks), "taxonomy_blinded": True, "coordinates_blinded": True,
    }]).to_csv(out_dir / "ontology_packet_qc.csv", index=False, encoding="utf-8-sig")
    print(f"[OK] wrote ontology blinded packets to {out_dir}")


if __name__ == "__main__":
    main()

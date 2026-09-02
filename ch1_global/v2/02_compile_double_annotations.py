#!/usr/bin/env python3
"""Compile blinded double annotations into agreement and adjudication outputs.

This script never joins taxonomy or coordinates back into annotator responses.
It compares only shared calibration units and writes a conflict queue that can
be adjudicated without revealing species identity or locality.

Example
-------
python ch1_global/v2/02_compile_double_annotations.py \
  --primary <annotator_1_packet_completed.csv> \
  --secondary <annotator_2_packet_completed.csv> \
  --out-dir <annotation_agreement>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


CATEGORICAL_TRAITS = [
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
]

ORDINAL_ORDERS = {
    "head_orientation_class": ["upright", "inclined", "nodding"],
    "flower_colour_class": ["white", "pale", "dark"],
    "involucral_cover_score_0_to_4": ["0", "1", "2", "3", "4"],
    "corolla_display_score_0_to_4": ["0", "1", "2", "3", "4"],
}

ALLOWED = {
    "flowering_qc": {"usable", "not_flowering", "too_small_or_blurred", "occluded", "wrong_taxon_or_non-Cirsium", "uncertain"},
    "detector_box_qc": {"accept", "revise", "reject", "not_assessable"},
    "camera_tilt_flag": {"vertical_reference_clear", "vertical_reference_uncertain", "tilted_or_rotated", "not_assessable"},
    "view_angle_class": {"lateral", "oblique", "frontal", "not_assessable"},
    "head_orientation_class": {"upright", "inclined", "nodding", "not_assessable"},
    "flower_colour_class": {"white", "pale", "dark", "not_assessable"},
    "colour_qc": {"usable", "illumination_or_white_balance_unusable", "not_assessable"},
    "involucral_cover_score_0_to_4": {"0", "1", "2", "3", "4", "not_assessable"},
    "head_shape_class": {"globose", "depressed", "cylindrical", "not_assessable"},
    "corolla_display_score_0_to_4": {"0", "1", "2", "3", "4", "not_assessable"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compile Ch.1 double annotations.")
    parser.add_argument("--primary", required=True, help="Completed primary packet CSV")
    parser.add_argument("--secondary", required=True, help="Completed secondary packet CSV")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--primary-name", default="annotator_1")
    parser.add_argument("--secondary-name", default="annotator_2")
    parser.add_argument(
        "--include-incomplete",
        action="store_true",
        help="Compare rows even when annotation_complete is not yes; default is completed rows only.",
    )
    return parser.parse_args()


def clean(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def validate_packet(df: pd.DataFrame, packet_name: str) -> pd.DataFrame:
    errors: list[dict[str, object]] = []
    if "annotation_unit_id" not in df.columns:
        raise ValueError(f"{packet_name}: missing annotation_unit_id")
    if "annotation_complete" not in df.columns:
        errors.append({"packet": packet_name, "annotation_unit_id": "", "trait": "annotation_complete", "value": "", "error": "missing_column"})
    duplicate_ids = df.loc[df["annotation_unit_id"].duplicated(keep=False), "annotation_unit_id"]
    for value in duplicate_ids.unique():
        errors.append({"packet": packet_name, "annotation_unit_id": value, "trait": "annotation_unit_id", "value": value, "error": "duplicate_id"})

    for trait, allowed in ALLOWED.items():
        if trait not in df.columns:
            errors.append({"packet": packet_name, "annotation_unit_id": "", "trait": trait, "value": "", "error": "missing_column"})
            continue
        values = df[trait].map(clean)
        invalid = ~values.isin(allowed | {""})
        for row_index in df.index[invalid]:
            errors.append({
                "packet": packet_name,
                "annotation_unit_id": df.at[row_index, "annotation_unit_id"],
                "trait": trait,
                "value": values.at[row_index],
                "error": "invalid_value",
            })
    return pd.DataFrame(errors, columns=["packet", "annotation_unit_id", "trait", "value", "error"])


def unweighted_kappa(left: pd.Series, right: pd.Series) -> float:
    if len(left) == 0:
        return float("nan")
    observed = float((left == right).mean())
    left_freq = left.value_counts(normalize=True)
    right_freq = right.value_counts(normalize=True)
    expected = sum(float(left_freq.get(label, 0.0) * right_freq.get(label, 0.0)) for label in set(left_freq.index) | set(right_freq.index))
    if np.isclose(1.0 - expected, 0.0):
        return float("nan")
    return (observed - expected) / (1.0 - expected)


def weighted_kappa(left: pd.Series, right: pd.Series, order: list[str]) -> float:
    if len(left) == 0:
        return float("nan")
    index = {value: i for i, value in enumerate(order)}
    if not set(left).issubset(index) or not set(right).issubset(index):
        return float("nan")
    k = len(order)
    matrix = np.zeros((k, k), dtype=float)
    for a, b in zip(left, right):
        matrix[index[a], index[b]] += 1.0
    observed = matrix / matrix.sum()
    expected = np.outer(observed.sum(axis=1), observed.sum(axis=0))
    weights = 1.0 - (np.abs(np.arange(k)[:, None] - np.arange(k)[None, :]) / max(k - 1, 1))
    p_observed = float((weights * observed).sum())
    p_expected = float((weights * expected).sum())
    if np.isclose(1.0 - p_expected, 0.0):
        return float("nan")
    return (p_observed - p_expected) / (1.0 - p_expected)


def main() -> None:
    args = parse_args()
    primary = pd.read_csv(args.primary, low_memory=False)
    secondary = pd.read_csv(args.secondary, low_memory=False)
    primary_errors = validate_packet(primary, args.primary_name)
    secondary_errors = validate_packet(secondary, args.secondary_name)
    errors = pd.concat([primary_errors, secondary_errors], ignore_index=True)

    if primary["annotation_unit_id"].duplicated().any() or secondary["annotation_unit_id"].duplicated().any():
        raise ValueError("Duplicate annotation_unit_id values cannot be compared")
    if "annotation_complete" not in primary.columns or "annotation_complete" not in secondary.columns:
        raise ValueError("Both packets need annotation_complete")

    image_cols = [column for column in ["source_image", "crop_path"] if column in primary.columns]
    left_columns = ["annotation_unit_id", "annotation_complete", *image_cols, *[t for t in CATEGORICAL_TRAITS if t in primary.columns]]
    right_columns = ["annotation_unit_id", "annotation_complete", *[t for t in CATEGORICAL_TRAITS if t in secondary.columns]]
    left = primary[left_columns].copy().add_suffix("_primary")
    left = left.rename(columns={"annotation_unit_id_primary": "annotation_unit_id"})
    right = secondary[right_columns].copy().add_suffix("_secondary")
    right = right.rename(columns={"annotation_unit_id_secondary": "annotation_unit_id"})
    merged = left.merge(right, on="annotation_unit_id", how="inner", validate="one_to_one")

    if not args.include_incomplete:
        keep = merged["annotation_complete_primary"].map(clean).eq("yes") & merged["annotation_complete_secondary"].map(clean).eq("yes")
        merged = merged.loc[keep].copy()

    summaries: list[dict[str, object]] = []
    conflicts: list[dict[str, object]] = []
    for trait in CATEGORICAL_TRAITS:
        left_name, right_name = f"{trait}_primary", f"{trait}_secondary"
        if left_name not in merged.columns or right_name not in merged.columns:
            continue
        pair = merged[["annotation_unit_id", left_name, right_name, *[f"{c}_primary" for c in image_cols]]].copy()
        pair[left_name] = pair[left_name].map(clean)
        pair[right_name] = pair[right_name].map(clean)
        pair = pair.loc[pair[left_name].ne("") & pair[right_name].ne("")].copy()
        values_left, values_right = pair[left_name], pair[right_name]
        assessable = pair.loc[values_left.ne("not_assessable") & values_right.ne("not_assessable")].copy()
        score_left, score_right = assessable[left_name], assessable[right_name]
        if trait in ORDINAL_ORDERS:
            kappa = weighted_kappa(score_left, score_right, ORDINAL_ORDERS[trait])
            kappa_kind = "linear_weighted_assessable_only"
        else:
            kappa = unweighted_kappa(values_left, values_right)
            kappa_kind = "unweighted"
        summaries.append({
            "trait": trait,
            "n_joint_completed": len(pair),
            "n_joint_assessable": len(assessable),
            "percent_agreement": float((values_left == values_right).mean()) if len(pair) else float("nan"),
            "cohen_kappa": kappa,
            "kappa_type": kappa_kind,
        })
        for _, row in pair.loc[values_left.ne(values_right)].iterrows():
            conflict = {
                "annotation_unit_id": row["annotation_unit_id"],
                "trait": trait,
                f"{args.primary_name}_value": row[left_name],
                f"{args.secondary_name}_value": row[right_name],
                "adjudicated_value": "",
                "adjudication_status": "pending",
                "adjudicator": "",
                "adjudication_notes": "",
            }
            for column in image_cols:
                conflict[column] = row.get(f"{column}_primary", "")
            conflicts.append(conflict)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(out_dir / "annotation_agreement_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(conflicts).to_csv(out_dir / "annotation_adjudication_queue.csv", index=False, encoding="utf-8-sig")
    errors.to_csv(out_dir / "annotation_validation_errors.csv", index=False, encoding="utf-8-sig")
    shared = pd.DataFrame([{
        "n_primary_rows": len(primary),
        "n_secondary_rows": len(secondary),
        "n_shared_units": len(merged),
        "n_conflicts": len(conflicts),
        "n_validation_errors": len(errors),
        "include_incomplete": args.include_incomplete,
    }])
    shared.to_csv(out_dir / "annotation_comparison_run_info.csv", index=False, encoding="utf-8-sig")
    print(f"[OK] Wrote agreement and adjudication outputs to {out_dir}")
    print(shared.to_string(index=False))


if __name__ == "__main__":
    main()

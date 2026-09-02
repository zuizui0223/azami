#!/usr/bin/env python3
"""Summarize a real Ch.1 iNaturalist metadata collection before image download.

This report is a QC gate, not a biological analysis. It records how many records
are eligible for screening and explicitly warns when a limited ID-ordered pilot
is being treated as a representative global sample.

Example
-------
python ch1_global/v2/08_summarize_metadata_collection.py \
  --metadata "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_raw_pilot\\photo_metadata.csv" \
  --collection-summary "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_raw_pilot\\collection_summary.csv" \
  --queue-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_screen_queue_pilot" \
  --out-dir "C:\\Users\\zuizui\\cirsium_inat\\ch1_v2_raw_pilot\\qc"
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="QC report for one real iNaturalist metadata collection.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--collection-summary", default="")
    parser.add_argument("--queue-dir", default="")
    parser.add_argument("--minimum-species-for-final", type=int, default=36)
    parser.add_argument("--technical-pilot", action="store_true", help="Marks output as nonrepresentative technical validation.")
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return value is True or text(value).lower() in {"true", "1", "yes"}


def summarize_series(df: pd.DataFrame, column: str, n: int = 30) -> pd.DataFrame:
    if column not in df.columns:
        return pd.DataFrame(columns=[column, "n", "fraction"])
    counts = df[column].map(text).replace("", "missing").value_counts(dropna=False).rename_axis(column).reset_index(name="n")
    counts["fraction"] = counts["n"] / max(1, len(df))
    return counts.head(n)


def main() -> None:
    args = parse_args()
    metadata_path = Path(args.metadata)
    data = pd.read_csv(metadata_path, low_memory=False, dtype={"obs_id": str, "photo_id": str})
    required = {"obs_id", "photo_id", "taxon_name", "taxon_rank", "coordinate_usable_for_environment", "inat_flowers_annotation"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required metadata columns: {sorted(missing)}")
    if data.empty:
        raise ValueError("Metadata input is empty")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data["coordinate_usable_bool"] = data["coordinate_usable_for_environment"].map(as_bool)
    data["flower_annotation_bool"] = data["inat_flowers_annotation"].map(as_bool)
    species_rows = data.loc[data["taxon_rank"].map(text).str.lower().eq("species")].copy()
    species_rows["taxon_name"] = species_rows["taxon_name"].map(text)
    open_geo_species = species_rows.loc[species_rows["coordinate_usable_bool"] & species_rows["taxon_name"].ne("")].copy()

    collection_summary = None
    summary_path = Path(args.collection_summary) if args.collection_summary else None
    if summary_path and summary_path.exists():
        collection_summary = pd.read_csv(summary_path).to_dict(orient="records")

    queue_metrics: dict[str, Any] = {}
    queue_dir = Path(args.queue_dir) if args.queue_dir else None
    if queue_dir and (queue_dir / "image_screening_queue.csv").exists():
        queue = pd.read_csv(queue_dir / "image_screening_queue.csv", low_memory=False)
        queue_metrics = {
            "queue_rows": int(len(queue)),
            "queue_species": int(queue["taxon_name"].nunique()) if "taxon_name" in queue.columns else 0,
            "queue_spatial_blocks": int(queue["spatial_block"].nunique()) if "spatial_block" in queue.columns else 0,
        }

    metrics = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "metadata_path": str(metadata_path.resolve()),
        "technical_pilot": bool(args.technical_pilot),
        "n_photo_rows": int(len(data)),
        "n_observations": int(data["obs_id"].nunique()),
        "n_species_rank_photo_rows": int(len(species_rows)),
        "n_species_total": int(species_rows["taxon_name"].nunique()),
        "n_open_georeferenced_species": int(open_geo_species["taxon_name"].nunique()),
        "n_open_georeferenced_observations": int(open_geo_species["obs_id"].nunique()),
        "coordinate_usable_fraction": float(data["coordinate_usable_bool"].mean()),
        "flower_annotation_fraction": float(data["flower_annotation_bool"].mean()),
        "collection_summary_file": collection_summary,
        **queue_metrics,
    }
    (out_dir / "collection_qc_metrics.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")

    summaries = {
        "taxon_rank_counts.csv": summarize_series(data, "taxon_rank"),
        "quality_grade_counts.csv": summarize_series(data, "quality_grade"),
        "photo_license_counts.csv": summarize_series(data, "photo_license_code"),
        "geoprivacy_counts.csv": summarize_series(data, "geoprivacy"),
        "top_species_photo_rows.csv": (
            species_rows.loc[species_rows["taxon_name"].ne("")]
            .groupby("taxon_name")
            .agg(n_photo_rows=("photo_id", "size"), n_observations=("obs_id", "nunique"), n_open_georeferenced=("coordinate_usable_bool", "sum"))
            .reset_index()
            .sort_values(["n_photo_rows", "taxon_name"], ascending=[False, True])
            .head(100)
        ),
    }
    for filename, table in summaries.items():
        table.to_csv(out_dir / filename, index=False, encoding="utf-8-sig")

    warnings: list[str] = []
    if args.technical_pilot:
        warnings.append("This is a technical pilot, not a geographically or taxonomically representative biological sample. Do not estimate trait frequencies, climate associations, or species prevalence from it.")
    if metrics["n_open_georeferenced_species"] < args.minimum_species_for_final:
        warnings.append(f"Only {metrics['n_open_georeferenced_species']} open, georeferenced species are present; below the final-stage target of {args.minimum_species_for_final}.")
    if metrics["coordinate_usable_fraction"] < 0.5:
        warnings.append("Fewer than half of photo rows have open usable coordinates. Reassess the environmental-analysis sampling frame before image download.")
    if metrics["flower_annotation_fraction"] < 0.05:
        warnings.append("Flower annotations are rare. This is expected for opportunistic records; retain them only as QC metadata and use human flowering assessment later.")
    if not queue_metrics:
        warnings.append("No screening queue was supplied. Build and inspect a balanced queue before downloading images.")

    lines = [
        "# Ch.1 iNaturalist metadata collection QC",
        "",
        "## Scope",
        "",
        "This report audits data availability and sampling-frame construction. It is not a trait result and does not support biological inference.",
        "",
        "## Core counts",
        "",
        f"- Photo metadata rows: **{metrics['n_photo_rows']:,}**",
        f"- Unique observations: **{metrics['n_observations']:,}**",
        f"- Species-rank taxa: **{metrics['n_species_total']:,}**",
        f"- Open, georeferenced species: **{metrics['n_open_georeferenced_species']:,}**",
        f"- Open-coordinate fraction: **{metrics['coordinate_usable_fraction']:.1%}**",
        f"- iNaturalist flower-annotation fraction: **{metrics['flower_annotation_fraction']:.1%}**",
    ]
    if queue_metrics:
        lines.extend([
            f"- Candidate queue rows: **{metrics['queue_rows']:,}**",
            f"- Candidate queue species: **{metrics['queue_species']:,}**",
            f"- Candidate queue spatial blocks: **{metrics['queue_spatial_blocks']:,}**",
        ])
    lines.extend(["", "## QC decision", ""])
    if warnings:
        lines.extend([f"- ⚠️ {warning}" for warning in warnings])
    else:
        lines.append("- ✅ No automatic blocking warning. Manually inspect the species and spatial QC tables before image download.")
    lines.extend([
        "",
        "## Output tables",
        "",
        "- `top_species_photo_rows.csv` — concentration of records by species.",
        "- `taxon_rank_counts.csv` — taxonomic-resolution composition.",
        "- `quality_grade_counts.csv`, `photo_license_counts.csv`, `geoprivacy_counts.csv` — provenance and usability context.",
    ])
    (out_dir / "collection_qc_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("[OK] wrote collection QC report", out_dir.resolve())
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

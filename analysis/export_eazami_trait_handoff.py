#!/usr/bin/env python3
"""Export the frozen Chapter 1 taxon-summary table for downstream EAzami work.

This exporter is intentionally a *summary bridge*: it preserves taxon medians,
MADs, coverage and image-proxy semantics for topology/crosswalk work, but it does
not preserve the full within-taxon observation distribution. Downstream analyses
that depend on within-taxon structure must use an observation-level export from
the frozen exhaustive Chapter 1 source (for example the Japan-38 exhaustive
handoff added in PR57).

The export does not assign phylogenetic placements, ancestral states or botanical
meanings to image proxies.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


PRIMARY_COLUMNS = [
    "taxon_name",
    "n_heads",
    "n_observations",
    "n_usable_colour",
    "n_usable_shape",
    "n_usable_orientation",
    "n_spatial_blocks",
    "corolla_lab_lightness_species_median",
    "corolla_lab_lightness_species_mad",
    "corolla_lab_chroma_species_median",
    "corolla_lab_chroma_species_mad",
    "corolla_hue_sin_species_median",
    "corolla_hue_sin_species_mad",
    "corolla_hue_cos_species_median",
    "corolla_hue_cos_species_mad",
    "orientation_angle_degrees_species_median",
    "orientation_angle_degrees_species_mad",
    "shape_aspect_ratio_species_median",
    "shape_aspect_ratio_species_mad",
    "shape_circularity_species_median",
    "shape_circularity_species_mad",
    "shape_solidity_species_median",
    "shape_solidity_species_mad",
    "shape_width_cv_species_median",
    "shape_width_cv_species_mad",
]

AUXILIARY_COLUMNS = [
    "n_heads_species",
    "n_usable_heads_species",
    "involucre_projection_roughness_species_median",
    "involucre_projection_roughness_species_mad",
    "involucre_spread_fraction_species_median",
    "involucre_spread_fraction_species_mad",
    "spine_peak_count_proxy_species_median",
    "spine_peak_count_proxy_species_mad",
    "spine_relative_length_max_proxy_species_median",
    "spine_relative_length_max_proxy_species_mad",
]

FINAL_AUXILIARY_INFERENTIAL_ENDPOINTS = [
    "involucre_projection_roughness",
    "involucre_spread_fraction",
    "spine_relative_length_max_proxy",
]

AUXILIARY_PROVENANCE_ONLY_ENDPOINTS = [
    "spine_peak_count_proxy",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_handoff(table: pd.DataFrame) -> pd.DataFrame:
    required = set(PRIMARY_COLUMNS + AUXILIARY_COLUMNS)
    missing = sorted(required.difference(table.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if table["taxon_name"].duplicated().any():
        raise ValueError("Species table must be unique by taxon_name")

    out = table[PRIMARY_COLUMNS + AUXILIARY_COLUMNS].copy()
    out.insert(1, "evidence_scope", "azami_ch1_public_image_taxon_summary")
    out.insert(2, "trait_state_scope", "continuous_taxon_summary_median_MAD")
    out.insert(3, "distribution_scope", "summary_bridge_not_full_within_taxon_distribution")
    out.insert(
        4,
        "distribution_requirement",
        "use_frozen_observation_level_export_when_within_taxon_structure_matters",
    )
    out.insert(5, "orientation_reference", "EXIF_oriented_image_vertical_not_gravity")
    out.insert(6, "colour_calibration", "uncalibrated_public_image_visible_colour")
    out.insert(7, "involucre_spine_proxy_status", "image_geometry_proxy_not_direct_botanical_state")
    out.insert(
        8,
        "auxiliary_inference_scope",
        "roughness_spread_max_projection_final;peak_count_provenance_only",
    )
    out.insert(9, "phylogeny_status", "unassigned_downstream_EAzami_crosswalk_required")
    return out.sort_values("taxon_name").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-table", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--source-artifact-id", required=True)
    parser.add_argument("--source-artifact-digest", required=True)
    args = parser.parse_args()

    source = Path(args.species_table)
    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(source)
    handoff = build_handoff(table)
    out_csv = outdir / "azami_ch1_eazami_trait_handoff_v1.csv"
    handoff.to_csv(out_csv, index=False, encoding="utf-8")

    summary = {
        "contract_version": "azami_ch1_eazami_trait_handoff_v1",
        "bridge_type": "taxon_summary_bridge",
        "full_within_taxon_distribution_preserved": False,
        "observation_level_distribution_required_when_within_taxon_structure_matters": True,
        "final_auxiliary_inferential_endpoints": FINAL_AUXILIARY_INFERENTIAL_ENDPOINTS,
        "auxiliary_provenance_only_endpoints": AUXILIARY_PROVENANCE_ONLY_ENDPOINTS,
        "n_taxa": int(len(handoff)),
        "n_taxa_with_colour": int(handoff["n_usable_colour"].fillna(0).gt(0).sum()),
        "n_taxa_with_shape": int(handoff["n_usable_shape"].fillna(0).gt(0).sum()),
        "n_taxa_with_orientation": int(handoff["n_usable_orientation"].fillna(0).gt(0).sum()),
        "n_taxa_with_auxiliary_involucre_spine": int(handoff["n_usable_heads_species"].fillna(0).gt(0).sum()),
        "source_species_table_sha256": sha256(source),
        "output_csv_sha256": sha256(out_csv),
        "source_artifact_id": str(args.source_artifact_id),
        "source_artifact_digest": args.source_artifact_digest,
        "claim_boundary": (
            "Taxon-level public-image median/MAD summaries for downstream crosswalk and topology-ensemble analyses. "
            "This file is not a replacement for the Chapter 1 observation-level within-taxon distribution. "
            "No ancestral state, adaptation, evolutionary rate, genetic variance or direct botanical "
            "phyllary/spine state is inferred by this export."
        ),
    }
    (outdir / "azami_ch1_eazami_trait_handoff_v1.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

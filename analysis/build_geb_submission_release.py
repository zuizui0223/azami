#!/usr/bin/env python3
"""Build the frozen Azami Chapter 1 GEB submission release."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from geb_submission_release_common import as_bool, endpoint_module_map, read_csv, require_file, setup_matplotlib, stable_pca
from geb_submission_release_main_figures import figure_4, figure_5, figure_6
from geb_submission_release_supplement import figure_s110, figure_s111, figure_s112, figure_s113, figure_s115, figure_s116, figure_s117
from geb_submission_release_tables import assemble_manuscript, build_readme, sha256, validate_inputs, write_tables


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--continuous-root", required=True, type=Path)
    p.add_argument("--cross-root", required=True, type=Path)
    p.add_argument("--manuscript-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--release-date", default="2026-08-27")
    p.add_argument("--continuous-run", type=int, default=32975451732)
    p.add_argument("--continuous-artifact", type=int, default=9612943217)
    p.add_argument("--continuous-digest", default="sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e")
    p.add_argument("--cross-run", type=int, default=33035785120)
    p.add_argument("--cross-artifact", type=int, default=9632715852)
    p.add_argument("--cross-digest", default="sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6")
    return p.parse_args()

def main() -> None:
    args = parse_args()
    setup_matplotlib()
    c = args.continuous_root
    x = args.cross_root
    out = args.out_dir
    if out.exists():
        shutil.rmtree(out)
    main_fig_dir = out / "figures" / "main"
    supp_fig_dir = out / "figures" / "supplement"
    table_dir = out / "tables"
    metadata_dir = out / "metadata"
    text_dir = out / "submission_text"
    for path in [main_fig_dir, supp_fig_dir, table_dir, metadata_dir, text_dir]:
        path.mkdir(parents=True, exist_ok=True)

    contract = read_csv(c / "universe" / "continuous_trait_contract_frozen.csv")
    geography = read_csv(c / "trait_geography" / "geb_v2_endpoint_geography_summary.csv")
    scores_all = read_csv(c / "trait_geography" / "geb_v2_taxon_morphospace_scores.csv")
    loadings_all = read_csv(c / "trait_geography" / "geb_v2_taxon_morphospace_loadings.csv")
    pca_report = json.loads(require_file(c / "trait_geography" / "geb_v2_trait_geography_report.json").read_text(encoding="utf-8"))
    atlas = read_csv(c / "geb_v2" / "geb_v2_environmental_evidence_atlas.csv")
    quality_adjusted = read_csv(c / "v2_quality" / "candidate_quality_adjusted_coefficients_v2.csv")
    quality_sensitivity = read_csv(c / "v2_quality" / "candidate_quality_sensitivity_v2.csv")
    species_long = read_csv(c / "universe" / "continuous_trait_universe_species_long.csv", low_memory=False)

    alignment = read_csv(x / "alignment" / "geb_v2_within_among_alignment.csv")
    unit_matrices = read_csv(x / "capitulum_space" / "capitulum_space_unit_strength_matrices.csv")
    module_integration = read_csv(x / "capitulum_space" / "capitulum_space_module_integration.csv")
    scope_summary = read_csv(x / "capitulum_space" / "capitulum_space_scope_summary.csv")
    eigenspectra = read_csv(x / "capitulum_space" / "capitulum_space_eigenspectra.csv")
    block_tests = read_csv(x / "capitulum_environment" / "capitulum_environment_block_tests.csv")
    geometry = read_csv(x / "capitulum_environment" / "capitulum_environment_cross_scale_geometry.csv")
    incremental = read_csv(x / "capitulum_environment_incremental" / "capitulum_environment_incremental_tests.csv")

    validation = validate_inputs(
        contract,
        geography,
        scores_all,
        loadings_all,
        atlas,
        alignment,
        scope_summary,
        unit_matrices,
        block_tests,
        geometry,
        incremental,
    )
    module_map = endpoint_module_map(contract)
    expanded_scores = scores_all[scores_all["scope"].eq("all_inferential_endpoints")].copy()
    expanded_loadings = loadings_all[loadings_all["scope"].eq("all_inferential_endpoints")].copy()
    all_report = next(item for item in pca_report["pca_scopes"] if item["scope"] == "all_inferential_endpoints")
    expanded_explained = np.asarray([all_report["explained_variance"][f"PC{i}"] for i in [1,2,3]], dtype=float)

    species_long["value"] = pd.to_numeric(species_long["value"], errors="coerce")
    species_long["available"] = as_bool(species_long["measurement_available"]) & species_long["value"].notna()
    primary_ids = contract.loc[contract["analysis_tier"].eq("primary"), "endpoint_id"].tolist()
    primary_wide = (
        species_long[species_long["endpoint_id"].isin(primary_ids) & species_long["available"]]
        .pivot(index="taxon_name", columns="endpoint_id", values="value")
    )
    primary = stable_pca(primary_wide)
    if len(primary.scores) != 249:
        raise ValueError(f"Expected 249 complete primary-PCA taxa, observed {len(primary.scores)}")

    figure_paths: list[Path] = []
    figure_paths += figure_4(expanded_scores, expanded_loadings, expanded_explained, module_map, main_fig_dir)
    figure_paths += figure_5(unit_matrices, scope_summary, contract, main_fig_dir)
    figure_paths += figure_6(atlas, alignment, incremental, contract, main_fig_dir)
    figure_paths += figure_s110(primary, expanded_scores, expanded_loadings, expanded_explained, module_map, supp_fig_dir)
    figure_paths += figure_s111(scores_all, pca_report, supp_fig_dir)
    figure_paths += figure_s112(quality_sensitivity, contract, supp_fig_dir)
    figure_paths += figure_s113(alignment, contract, supp_fig_dir)
    figure_paths += figure_5(
        unit_matrices,
        scope_summary,
        contract,
        supp_fig_dir,
        scope="complete18_min2",
        number="Figure_S1_14_complete18_min2_organization",
        title_prefix="Figure S1.14",
    )
    figure_paths += figure_s115(block_tests, supp_fig_dir)
    figure_paths += figure_s116(incremental, supp_fig_dir)
    figure_paths += figure_s117(geometry, supp_fig_dir)

    table_index = write_tables(
        c,
        x,
        contract,
        geography,
        scores_all,
        loadings_all,
        quality_adjusted,
        quality_sensitivity,
        alignment,
        scope_summary,
        module_integration,
        unit_matrices,
        eigenspectra,
        block_tests,
        geometry,
        incremental,
        atlas,
        pca_report,
        table_dir,
    )
    metrics = assemble_manuscript(args.manuscript_dir, text_dir)
    if int(metrics["abstract_body_words"]) > 300:
        raise ValueError(f"Structured abstract exceeds 300 words: {metrics['abstract_body_words']}")

    figure_index_rows = []
    for path in sorted(figure_paths):
        figure_index_rows.append(
            {
                "filename": str(path.relative_to(out)),
                "format": path.suffix.lstrip("."),
                "bytes": path.stat().st_size,
            }
        )
    pd.DataFrame(figure_index_rows).to_csv(metadata_dir / "FIGURE_INDEX.csv", index=False)

    validation.update(
        {
            "primary_pca_complete_taxa": int(len(primary.scores)),
            "primary_pca_explained_variance": {f"PC{i+1}": float(v) for i, v in enumerate(primary.explained)},
            "expanded_pca_explained_variance": {f"PC{i+1}": float(v) for i, v in enumerate(expanded_explained)},
            "generated_main_figures": 3,
            "generated_supplement_figures": 8,
            "generated_figure_files": len(figure_paths),
            "supplement_table_parts": int(len(table_index)),
            "abstract_body_words": int(metrics["abstract_body_words"]),
        }
    )
    (metadata_dir / "release_validation_report.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")

    release_manifest = {
        "schema_version": 1,
        "release_date": args.release_date,
        "purpose": "GEB submission presentation release from frozen Azami Chapter 1 artifacts",
        "frozen_sources": {
            "continuous": {"run": args.continuous_run, "artifact": args.continuous_artifact, "digest": args.continuous_digest},
            "cross_scale": {"run": args.cross_run, "artifact": args.cross_artifact, "digest": args.cross_digest},
        },
        "release_counts": {
            "main_figures": 3,
            "supplement_figures": 8,
            "formats_per_figure": 3,
            "supplement_table_parts": int(len(table_index)),
        },
        "computational_boundary": (
            "Presentation-only release plus exact-method descriptive reconstruction of the nine-primary PCA. "
            "No inferential family, decision, endpoint, cohort or evidence grade is changed."
        ),
        "submission_status": "HOLD for six independent scientific validation gates",
    }
    (metadata_dir / "release_manifest.json").write_text(json.dumps(release_manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "README.md").write_text(build_readme(args, metrics), encoding="utf-8")

    all_files = sorted(path for path in out.rglob("*") if path.is_file() and path.name != "SHA256SUMS.txt")
    hash_lines = [f"{sha256(path)}  {path.relative_to(out).as_posix()}" for path in all_files]
    (metadata_dir / "SHA256SUMS.txt").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Tables, manuscript assembly and validation for the GEB release."""
from __future__ import annotations

from geb_submission_release_common import *

def copy_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)


def write_tables(
    continuous_root: Path,
    cross_root: Path,
    contract: pd.DataFrame,
    geography: pd.DataFrame,
    scores: pd.DataFrame,
    loadings: pd.DataFrame,
    quality_adjusted: pd.DataFrame,
    quality_sensitivity: pd.DataFrame,
    alignment: pd.DataFrame,
    scope_summary: pd.DataFrame,
    module_integration: pd.DataFrame,
    unit_matrices: pd.DataFrame,
    eigenspectra: pd.DataFrame,
    block_tests: pd.DataFrame,
    geometry: pd.DataFrame,
    incremental: pd.DataFrame,
    atlas: pd.DataFrame,
    pca_report: dict[str, object],
    out: Path,
) -> pd.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    obs_path = continuous_root / "universe" / "continuous_trait_universe_observation_long.csv"
    obs = read_csv(obs_path, usecols=["taxon_name", "endpoint_id", "measurement_available", "value"], low_memory=False)
    obs["available"] = as_bool(obs["measurement_available"]) & pd.to_numeric(obs["value"], errors="coerce").notna()
    measured = (
        obs[obs["available"]]
        .groupby("endpoint_id")
        .agg(n_available_observations=("endpoint_id", "size"), n_taxa_with_measurement=("taxon_name", "nunique"))
        .reset_index()
    )
    coverage = contract.merge(measured, on="endpoint_id", how="left")
    coverage["n_available_observations"] = coverage["n_available_observations"].fillna(0).astype(int)
    coverage["n_taxa_with_measurement"] = coverage["n_taxa_with_measurement"].fillna(0).astype(int)
    coverage["measurement_executed"] = coverage["n_available_observations"].gt(0)
    geo_keep = geography.drop(columns=["module", "analysis_tier", "unit"], errors="ignore")
    coverage = coverage.merge(geo_keep, on="endpoint_id", how="left", suffixes=("", "_inferential"))

    parts: list[tuple[str, str, str, pd.DataFrame, str]] = [
        ("S1.15", "a", "S1_15_endpoint_coverage.csv", coverage, "27-endpoint contract plus frozen execution coverage"),
        ("S1.16", "a", "S1_16a_endpoint_geography.csv", geography, "expanded endpoint geography summary"),
        ("S1.16", "b", "S1_16b_taxon_morphospace_scores.csv", scores, "taxon morphospace scores"),
        ("S1.16", "c", "S1_16c_taxon_morphospace_loadings.csv", loadings, "taxon morphospace loadings"),
        ("S1.17", "a", "S1_17a_candidate_quality_adjusted.csv", quality_adjusted, "candidate quality-adjusted coefficients"),
        ("S1.17", "b", "S1_17b_candidate_resolution_sensitivity.csv", quality_sensitivity, "candidate resolution-stratum sensitivity"),
        ("S1.18", "a", "S1_18_matched_within_among.csv", alignment, "68 matched within/among endpoint-climate rows"),
        ("S1.19", "a", "S1_19a_capitulum_scope_summary.csv", scope_summary, "whole-capitulum scope summary"),
        ("S1.19", "b", "S1_19b_module_integration.csv", module_integration, "registered-module contrasts"),
        ("S1.19", "c", "S1_19c_unit_strength_matrices.csv", unit_matrices, "17-unit association matrices"),
        ("S1.19", "d", "S1_19d_eigenspectra.csv", eigenspectra, "whole-capitulum eigenspectra"),
        ("S1.20", "a", "S1_20a_environment_block_tests.csv", block_tests, "24 environmental-block tests"),
        ("S1.20", "b", "S1_20b_cross_scale_geometry.csv", geometry, "12 coefficient-matrix cosine rows"),
        ("S1.21", "a", "S1_21_incremental_tests.csv", incremental, "20 nested core-four sufficiency tests"),
        ("S1.22", "a", "S1_22_environmental_evidence_atlas.csv", atlas, "31-row evidence atlas"),
    ]
    index_rows: list[dict[str, object]] = []
    for table_id, part, filename, frame, source in parts:
        path = out / filename
        copy_table(frame, path)
        index_rows.append(
            {
                "table_id": table_id,
                "part": part,
                "filename": filename,
                "n_rows": int(len(frame)),
                "n_columns": int(len(frame.columns)),
                "content": source,
            }
        )
    pca_rows = []
    for item in pca_report["pca_scopes"]:
        row = {k: v for k, v in item.items() if k != "explained_variance"}
        for key, value in item.get("explained_variance", {}).items():
            row[f"explained_{key}"] = value
        pca_rows.append(row)
    pca_frame = pd.DataFrame(pca_rows)
    copy_table(pca_frame, out / "S1_16d_PCA_scope_summary.csv")
    index_rows.append(
        {
            "table_id": "S1.16",
            "part": "d",
            "filename": "S1_16d_PCA_scope_summary.csv",
            "n_rows": int(len(pca_frame)),
            "n_columns": int(len(pca_frame.columns)),
            "content": "PCA scope and explained-variance summary",
        }
    )
    index = pd.DataFrame(index_rows)
    index.to_csv(out / "SUPPLEMENT_TABLE_INDEX.csv", index=False)
    return index


def markdown_word_count(text: str) -> int:
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = re.sub(r"\[[^\]]+\]\([^\)]+\)", " ", text)
    text = re.sub(r"[*_#>|]", " ", text)
    return len(re.findall(r"\b[\wÀ-ÖØ-öø-ÿ–-]+\b", text))


def abstract_body_word_count(text: str) -> int:
    lines: list[str] = []
    in_abstract = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped == "## Abstract":
            in_abstract = True
            continue
        if not in_abstract:
            continue
        if stripped.startswith("## "):
            break
        if stripped.startswith("### "):
            continue
        if stripped.startswith("**Keywords:") or stripped.startswith("**Running title:"):
            continue
        lines.append(line)
    return markdown_word_count("\n".join(lines))


def assemble_manuscript(manuscript_dir: Path, out: Path) -> dict[str, object]:
    order = [
        "00_title_abstract.md",
        "01_introduction.md",
        "02_methods.md",
        "03_results.md",
        "04_discussion.md",
        "05_conclusion_and_declarations.md",
        "07_data_code_availability.md",
        "06_references.md",
        "MAIN_FIGURE_CAPTIONS.md",
    ]
    sections = []
    metrics: dict[str, object] = {}
    for name in order:
        path = require_file(manuscript_dir / name)
        text = path.read_text(encoding="utf-8")
        sections.append(f"<!-- SOURCE: {name} -->\n\n{text.strip()}\n")
        metrics[name] = markdown_word_count(text)
    assembled = "\n\n---\n\n".join(sections) + "\n"
    out.mkdir(parents=True, exist_ok=True)
    (out / "GEB_BLINDED_MANUSCRIPT_SOURCE.md").write_text(assembled, encoding="utf-8")
    abstract = (manuscript_dir / "00_title_abstract.md").read_text(encoding="utf-8")
    metrics["abstract_body_words"] = abstract_body_word_count(abstract)
    main_names = ["01_introduction.md", "02_methods.md", "03_results.md", "04_discussion.md", "05_conclusion_and_declarations.md"]
    metrics["main_text_words"] = int(sum(int(metrics[name]) for name in main_names))
    refs = (manuscript_dir / "06_references.md").read_text(encoding="utf-8")
    metrics["reference_entries_approx"] = len(re.findall(r"^\s*[-*]\s+", refs, flags=re.M))
    lower = assembled.lower()
    scan_terms = ["zuizui0223", "rachelzhang", "@gmail.com", "github.com/zuizui0223"]
    metrics["anonymization_scan_hits"] = {term: lower.count(term) for term in scan_terms}
    (out / "submission_text_metrics.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")
    return metrics


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_inputs(
    contract: pd.DataFrame,
    geography: pd.DataFrame,
    scores: pd.DataFrame,
    loadings: pd.DataFrame,
    atlas: pd.DataFrame,
    alignment: pd.DataFrame,
    scope_summary: pd.DataFrame,
    unit_matrices: pd.DataFrame,
    block_tests: pd.DataFrame,
    geometry: pd.DataFrame,
    incremental: pd.DataFrame,
) -> dict[str, object]:
    checks: dict[str, tuple[object, object]] = {
        "contract_endpoints": (len(contract), 27),
        "measured_inferential_endpoint_geography_rows": (len(geography), 18),
        "expanded_complete_taxa": (int(scores[scores["scope"].eq("all_inferential_endpoints")]["taxon_name"].nunique()), 127),
        "expanded_loading_rows": (len(loadings[loadings["scope"].eq("all_inferential_endpoints")]), 18),
        "evidence_atlas_rows": (len(atlas), 31),
        "matched_rows": (len(alignment), 68),
        "whole_capitulum_scopes": (len(scope_summary), 2),
        "unit_strength_rows": (len(unit_matrices), 544),
        "environment_block_rows": (len(block_tests), 24),
        "cross_scale_geometry_rows": (len(geometry), 12),
        "incremental_test_rows": (len(incremental), 20),
    }
    class_counts = alignment["cross_scale_class"].value_counts().to_dict()
    expected_classes = {"both_scales": 3, "within_only": 8, "among_only": 1, "neither": 56}
    failures = {key: {"observed": obs, "expected": exp} for key, (obs, exp) in checks.items() if obs != exp}
    if class_counts != expected_classes:
        failures["cross_scale_class_counts"] = {"observed": class_counts, "expected": expected_classes}
    main = scope_summary[scope_summary["scope"].eq("complete18_min5")].iloc[0]
    numeric_expected = {
        "within_module_contrast": 0.1645023304673242,
        "among_module_contrast": 0.08847536372583811,
        "cross_scale_strength_matrix_spearman": 0.36629931778064023,
    }
    for key, expected in numeric_expected.items():
        observed = float(main[key])
        if not math.isclose(observed, expected, rel_tol=0, abs_tol=1e-12):
            failures[key] = {"observed": observed, "expected": expected}
    if failures:
        raise ValueError("Frozen input validation failed:\n" + json.dumps(failures, indent=2, default=str))
    return {
        "status": "passed",
        "checks": {key: {"observed": obs, "expected": exp} for key, (obs, exp) in checks.items()},
        "cross_scale_class_counts": class_counts,
        "main_scope_values": {key: float(main[key]) for key in numeric_expected},
    }


def build_readme(args: argparse.Namespace, metrics: dict[str, object]) -> str:
    return f"""# Azami Chapter 1 GEB submission release

**Release date:** {args.release_date}

This packet renders the final presentation layer from two immutable workflow artifacts. It does not change any numerical result, cohort, endpoint definition, multiplicity family or evidence grade.

## Frozen provenance

- continuous-trait run `{args.continuous_run}`, artifact `{args.continuous_artifact}`, digest `{args.continuous_digest}`;
- multilevel/process run `{args.cross_run}`, artifact `{args.cross_artifact}`, digest `{args.cross_digest}`.

## Contents

- `figures/main/`: final main Figures 4–6 in PNG, PDF and SVG;
- `figures/supplement/`: Figures S1.10–S1.17 in PNG, PDF and SVG;
- `tables/`: Tables S1.15–S1.22, with numbered multi-part exports where required;
- `submission_text/`: assembled blinded manuscript source and word-count/anonymization audit;
- `metadata/`: validation report, figure/table indices, release manifest and SHA-256 hashes.

## Text audit

- structured abstract body: {metrics['abstract_body_words']} words;
- main text from Introduction through declarations: approximately {metrics['main_text_words']} words.

## Claim boundary

The release reconstructs the present phenotypic field. It does not establish independent detector or trait accuracy, functional or genetic modularity, plasticity, adaptation, environmental causation or resolved evolutionary history. The six independent scientific validation gates recorded in `manuscript/EXTERNAL_COMPLETION_GATES.md` remain submission blocking.
"""

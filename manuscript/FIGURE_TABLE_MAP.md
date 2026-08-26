# Submission figure, table and supplement map

This file records the submission-facing relationship among the main manuscript, Appendix S1 and the repository assets from which the figures and tables are built. The main narrative contains five figures and two tables. GEB v2 changes Figures 4–5 from a nine-endpoint-only presentation to the expanded continuous-trait morphospace and evidence atlas; historical-placement detail remains in Appendix S1.

## Main figures

1. **Figure 1 — image to multiscale analysis:** open-licensed photographs → capitulum localization → established and expanded continuous measurements → repeated numerical phenotypes → nested variance, morphospace and within-/among-taxon environmental analyses.
2. **Figure 2 — geographic sampling and analytical domain:** detector-positive observations, the 46,276-observation spatially thinned primary cohort, balanced/exhaustive streams, the trait-blind expanded high-resolution intersection and CHELSA environmental domain.
3. **Figure 3 — nested visible variance:** source-assigned taxon → photograph → head decomposition across the nine established primary endpoints. This remains the frozen quantitative demonstration that taxon means discard substantial repeated visible variation.
4. **Figure 4 — expanded taxon-level morphospace:** 127 taxa complete for all 18 measured inferential endpoints, with standardized trait loadings and PC1–PC3 variance of 18.5%, 12.0% and 11.8%. The purpose is to show multidimensional phenotype structure rather than one syndrome.
5. **Figure 5 — environmental evidence atlas:** established primary effect sizes and spatial context together with the three registry-wide expanded candidate signals, annotated by A/B/C evidence grade and range/image-quality context.

## Main tables

1. **Table 1:** canonical cohorts and analysis roles.
2. **Table 2:** frozen 27-endpoint continuous-trait contract, analysis tier, module, measurement scope and current execution status.

## Appendix S1 tables

Submission labels are listed first; repository filenames retain their stable source identifiers. Existing primary-layer tables keep their stable labels, while the completed v2 artifact adds an expanded-results packet that should be exported into the next supplement build.

- **Table S1.1:** `S01_canonical_cohorts.csv` — canonical cohorts and inferential roles.
- **Table S1.2:** `S02_trait_scope_and_assessability.csv` — established trait scope and assessability boundaries; the 27-endpoint contract is the v2 extension used for Table 2.
- **Table S1.3:** `S03_primary_within_taxon_BH_supported_rows.csv` — eight BH-supported primary within-taxon component rows.
- **Table S1.4:** `S04_nested_visible_variance.csv` — nine-endpoint nested visible-variance decomposition and sampling sensitivities.
- **Table S1.5:** `S05_grouped_SPDE_stable_supported_patterns.csv` — stable grouped SPDE-INLA supported directions.
- **Table S1.6:** `S06_auxiliary_involucre_resolution_audit.csv` — legacy three-endpoint contour screen and failed confirmatory image-quality audit; retained as provenance rather than erased.
- **Table S1.7:** `S07_niche_permutation_summary.csv` — niche-permutation support.
- **Table S1.8:** `S08_randomized_PGLS_retained_rows.csv` — randomized-tree PGLS retained rows.
- **Table S1.9:** `S09_spatial_robustness.csv` — residual spatial and broad-region-omission robustness.
- **Table S1.10:** `S10_wcvp_synonym_candidates.csv` — WCVP synonym-collapse sensitivity.
- **Table S1.11:** `S13_primary_bias_control_audit.csv` — cyclic-season and dominant-taxon audit of the four non-circular primary rows.
- **Table S1.12:** `S14_repeat_photo_cohort_preflight.csv` — outcome-blind repeat-photo cohort counts and interpretation boundary.
- **Table S1.13:** `S15_hemisphere_season_sensitivity.csv` — phase-aligned and hemisphere-specific cyclic-season sensitivity for the four non-circular primary rows.
- **Table S1.14:** `S16_native_range_sensitivity.csv` — native-only coefficients and frozen retention decisions for the four non-circular primary rows.
- **GEB v2 expanded trait-universe table:** export endpoint coverage and taxon/observation counts from `continuous_trait_universe_report.json`.
- **GEB v2 trait-geography table:** export endpoint geography and all-inferential/module PCA summaries from `geb_v2_trait_geography_report.json` and associated CSVs.
- **GEB v2 candidate-quality table:** export the 36 quality-adjusted candidate climate rows and fixed resolution-stratum sensitivities, retaining the distinction between registry-wide candidate FDR support and quality-adjusted sensitivity rows.
- **GEB v2 environmental evidence atlas:** export the 31-row evidence atlas with A/B/C/D context and claim boundary.

The immutable final-run provenance and submission-facing values are recorded in [`results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](results/GEB_V2_FINAL_RESULTS_2026-08-26.md).

## Appendix S1 figures

The existing primary-layer filenames keep their stable asset numbers. The expanded figure build should preserve those diagnostics and add v2-specific supplementary views rather than overwrite the historical audit trail.

- **Figure S1.1:** `Figure_S1_nested_variance_sensitivities.*` — full, one-head-per-photo and balanced-10-photo below-taxon variance sensitivities.
- **Figure S1.2:** `Figure_S2_primary_coefficient_map.*` — complete primary within-taxon CHELSA coefficient map.
- **Figure S1.3:** `Figure_S3_spde_delta_waic.*` — endpoint-specific ΔWAIC across four grouped SPDE-INLA predictor sets.
- **Figure S1.4:** `Figure_S4_residual_morans_I.*` — residual Moran's I after spatial modelling.
- **Figure S1.5:** `Figure_S5_leave_one_region_out_stability.*` — minimum taxon-rank stability under broad-region omission.
- **Figure S1.6:** `Figure_S6_sampling_geography_across_filters.*` — Mollweide equal-area sampling density across detector-positive, coordinate-usable, ≤10 km and primary-thinned stages.
- **Figure S1.7:** `Figure_S7_geographic_trait_assessability.*` — Mollweide equal-area cell-level image-trait assessability.
- **Figure S1.8:** `Figure_S8_environmental_sorting.*` — environmental separation and overlap between low- and high-trait taxa.
- **Figure S1.9:** `Figure_S1_9_historical_placement_sensitivity.*` — randomized-tree PGLS placement sensitivity.
- **Expanded morphospace supplement:** module-specific PCAs and the frozen nine-primary-endpoint PCA for direct comparison with Figure 4.
- **Expanded candidate-quality supplement:** full quality-adjusted coefficients and 150–199, 200–299 and ≥300-pixel strata supporting the evidence grade shown in Figure 5.

Figure provenance and claim boundaries remain recorded in `manuscript/supplement/FIGURE_MANIFEST.md`, `manuscript/FIGURE2_GEOGRAPHY_MANIFEST.md` and `manuscript/INTERPRETIVE_FIGURES_MANIFEST.md`; those manifests should be regenerated when the v2 Figure 4–5 assets are rendered.

## Main manuscript ↔ Appendix S1 crosswalk

| Main item / result layer | Appendix S1 support | Purpose |
|---|---|---|
| Figure 1 — image to multiscale analysis | Table S1.2 + frozen 27-endpoint contract | defines measurement scope, tier and assessability boundaries |
| Figure 2 — geographic sampling and analytical domain | Table S1.1; Figures S1.6–S1.7 | records cohort counts and visualizes filtering and assessability |
| Figure 3 — nested visible variance | Table S1.4; Figure S1.1 | gives the complete established-endpoint decomposition and sampling sensitivities |
| Figure 4 — expanded taxon morphospace | v2 trait-geography report + module PCA supplement | establishes the 18-endpoint multidimensional phenotype map and shows module structure |
| Figure 5 — environmental evidence atlas | Tables S1.3, S1.5, S1.11, S1.13–S1.14 + v2 candidate-quality/evidence-atlas exports | separates A/B established evidence from C-grade expanded signals and makes the relevant sensitivity context visible |
| Table 1 | Table S1.1 | complete cohort audit |
| Table 2 | 27-endpoint contract validation | complete trait, tier and execution-status audit |
| Among-taxon environmental sorting | Table S1.7; Figure S1.8 | reports permutation support and observed environmental separation |
| Spatial robustness | Table S1.9; Figures S1.4–S1.5 | reports residual autocorrelation and region-omission diagnostics |
| Taxonomic robustness | Table S1.10 | reports WCVP synonym-collapse sensitivity |
| Legacy involucre contour screen | Table S1.6 | preserves the failed resolution/sharpness audit without confusing it with the expanded v2 candidate definitions |
| Expanded candidate geometry | v2 candidate-quality and evidence-atlas exports | distinguishes two quality-robust exploratory rows from one image-sensitive row |
| Repeat-photo preflight | Table S1.12 | records the available repeated-image cohort without implying a completed variance estimate |
| Native-range sensitivity | Table S1.14 | reports which established primary rows survive the versioned native-only restriction |
| Historical placement sensitivity | Table S1.8; Figure S1.9 | retains historical-placement provenance outside the main figure sequence |

Withdrawn lability analyses and the old quadrant figure are retained only as repository provenance and are not submitted as biological results. EAzami mechanism-reduction outputs are outside the Chapter 1 input set.

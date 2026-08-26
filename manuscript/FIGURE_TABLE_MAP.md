# Submission figure, table and supplement map

This file records the submission-facing relationship among the main manuscript, Appendix S1 and the repository assets from which the figures and tables are built. The current main manuscript contains five figures and two tables. The historical-placement panel has been demoted to Appendix S1.

## Main figures

1. **Figure 1 — image to multiscale analysis:** open-licensed photographs → capitulum localization → continuous orientation, colour and outline measurements → repeated numerical phenotypes → nested variance and within- and among-taxon environmental analyses.
2. **Figure 2 — geographic sampling and analytical domain:** (a) 1° density of detector-positive observations before spatial filtering; (b) 1° display-cell density of the 46,276-observation spatially thinned primary cohort; (c) balanced-atlas and exhaustive-spatial analysis streams; (d) CHELSA BIO1 × BIO12 environmental domain of the primary cohort. Map panels use the Mollweide equal-area projection and a 5,000-km equatorial scale bar.
3. **Figure 3 — nested visible variance:** source-assigned taxon → photograph → head decomposition across all nine primary endpoints.
4. **Figure 4 — taxon-level multivariate trait architecture:** PC1–PC2 positions of the 148 complete taxa with nine trait-loading vectors, plus PC3 loadings.
5. **Figure 5 — environmental effect sizes and bias controls:** (a) eight frozen BH-supported primary components; (b) grouped SPDE-INLA patterns; (c) native-only coefficients and frozen retention decisions after the seasonal and dominant-taxon audits.

## Main tables

1. **Table 1:** canonical cohorts and analysis roles.
2. **Table 2:** trait modules, endpoints and bounded measurement scope.

## Appendix S1 tables

Submission labels are listed first; repository filenames retain their stable source identifiers.

- **Table S1.1:** `S01_canonical_cohorts.csv` — canonical cohorts and inferential roles.
- **Table S1.2:** `S02_trait_scope_and_assessability.csv` — trait scope and assessability boundaries.
- **Table S1.3:** `S03_primary_within_taxon_BH_supported_rows.csv` — eight BH-supported primary within-taxon component rows.
- **Table S1.4:** `S04_nested_visible_variance.csv` — nine-endpoint nested visible-variance decomposition and sampling sensitivities.
- **Table S1.5:** `S05_grouped_SPDE_stable_supported_patterns.csv` — stable grouped SPDE-INLA supported directions.
- **Table S1.6:** `S06_auxiliary_involucre_resolution_audit.csv` — original and image-quality-adjusted high-resolution involucre rows; all are withdrawn from the headline.
- **Table S1.7:** `S07_niche_permutation_summary.csv` — niche-permutation support.
- **Table S1.8:** `S08_randomized_PGLS_retained_rows.csv` — randomized-tree PGLS retained rows.
- **Table S1.9:** `S09_spatial_robustness.csv` — residual spatial and broad-region-omission robustness.
- **Table S1.10:** `S10_wcvp_synonym_candidates.csv` — WCVP synonym-collapse sensitivity.
- **Table S1.11:** `S13_primary_bias_control_audit.csv` — cyclic-season and dominant-taxon audit of the four non-circular primary rows.
- **Table S1.12:** `S14_repeat_photo_cohort_preflight.csv` — outcome-blind repeat-photo cohort counts and interpretation boundary.
- **Table S1.13:** `S15_hemisphere_season_sensitivity.csv` — phase-aligned and hemisphere-specific cyclic-season sensitivity for the four non-circular primary rows.
- **Table S1.14:** `S16_native_range_sensitivity.csv` — native-only coefficients and frozen retention decisions for the four non-circular primary rows.

Repository QA files `S11_precision_confounding_audit_summary.csv` and `S12_submission_gates.csv` are not included in the submitted Appendix S1.

## Appendix S1 figures

The filenames keep their stable asset numbers; the submission labels follow the journal's Appendix S1 convention.

- **Figure S1.1:** `Figure_S1_nested_variance_sensitivities.*` — full, one-head-per-photo and balanced-10-photo below-taxon variance sensitivities.
- **Figure S1.2:** `Figure_S2_primary_coefficient_map.*` — complete 36-component within-taxon CHELSA coefficient map.
- **Figure S1.3:** `Figure_S3_spde_delta_waic.*` — endpoint-specific ΔWAIC across four grouped SPDE-INLA predictor sets.
- **Figure S1.4:** `Figure_S4_residual_morans_I.*` — residual Moran's I after spatial modelling.
- **Figure S1.5:** `Figure_S5_leave_one_region_out_stability.*` — minimum taxon-rank stability under broad-region omission.
- **Figure S1.6:** `Figure_S6_sampling_geography_across_filters.*` — Mollweide equal-area sampling density across detector-positive, coordinate-usable, ≤10 km and primary-thinned stages.
- **Figure S1.7:** `Figure_S7_geographic_trait_assessability.*` — Mollweide equal-area cell-level image-trait assessability.
- **Figure S1.8:** `Figure_S8_environmental_sorting.*` — environmental separation and overlap between low- and high-trait taxa.
- **Figure S1.9:** `Figure_S1_9_historical_placement_sensitivity.*` — within/among-taxon comparison and randomized-tree PGLS placement sensitivity, demoted from the main manuscript.

Figure provenance and claim boundaries are recorded in `manuscript/supplement/FIGURE_MANIFEST.md`, `manuscript/FIGURE2_GEOGRAPHY_MANIFEST.md` and `manuscript/INTERPRETIVE_FIGURES_MANIFEST.md`.

## Main manuscript ↔ Appendix S1 crosswalk

| Main item / result layer | Appendix S1 support | Purpose |
|---|---|---|
| Figure 1 — image to multiscale analysis | Table S1.2 | defines measurement scope and assessability boundaries |
| Figure 2 — geographic sampling and analytical domain | Table S1.1; Figures S1.6–S1.7 | records cohort counts and visualizes filtering and assessability |
| Figure 3 — nested visible variance | Table S1.4; Figure S1.1 | gives the complete endpoint table and sampling sensitivities |
| Figure 4 — taxon-level PCA | Table S1.2 | keeps trait definitions aligned with the displayed multivariate structure |
| Figure 5 — within-taxon environmental structure | Tables S1.3, S1.5, S1.11, S1.13–S1.14; Figures S1.2–S1.3 | separates frozen coefficients, grouped-spatial context and the final predeclared bias-control decisions |
| Table 1 | Table S1.1 | complete cohort audit |
| Table 2 | Table S1.2 | complete trait and assessability audit |
| Among-taxon environmental sorting | Table S1.7; Figure S1.8 | reports permutation support and observed environmental separation |
| Spatial robustness | Table S1.9; Figures S1.4–S1.5 | reports residual autocorrelation and region-omission diagnostics |
| Taxonomic robustness | Table S1.10 | reports WCVP synonym-collapse sensitivity |
| Withdrawn involucre rows | Table S1.6 | preserves the failed resolution/sharpness audit without retaining the claim |
| Repeat-photo preflight | Table S1.12 | records the available repeated-image cohort without implying a completed variance estimate |
| Native-range sensitivity | Table S1.14 | reports which primary rows survive the versioned native-only restriction |
| Historical placement sensitivity | Table S1.8; Figure S1.9 | retains provenance outside the main figure sequence |

Withdrawn lability analyses and the old quadrant figure are retained only as repository provenance and are not submitted as biological results. EAzami mechanism-reduction outputs are outside the Chapter 1 input set.

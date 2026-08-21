# Chapter 1 figure, table and supplement map

This file fixes the submission-facing relationship among the main manuscript, frozen Supplement tables and release figures. Items are mapped by analysis role so cohorts or inferential families cannot be silently interchanged.

## Main figures

1. **Figure 1 — image to multiscale analysis:** actual open-licensed photographs → YOLO capitulum localization → deterministic continuous orientation/colour/outline measurements → repeated numerical phenotypes → nested variance, within-taxon environment and among-taxon niche/historical questions.
2. **Figure 2 — geographic sampling and analytical domain:** (A) 1° density of detector-positive observations before spatial filtering; (B) the 46,276-observation spatially thinned primary cohort; (C) balanced-atlas and exhaustive-spatial analysis streams; (D) CHELSA BIO1 × BIO12 environmental domain of the primary cohort.
3. **Figure 3 — nested visible variance:** assigned taxon → photograph → head decomposition across all nine primary endpoints.
4. **Figure 4 — taxon-level multivariate trait architecture:** actual PC1–PC2 positions of the 148 complete taxa with nine trait-loading vectors, plus PC3 loadings to expose the third independent trait dimension.
5. **Figure 5 — environmental effect sizes across executed analyses:** (A) eight BH-supported exhaustive primary within-taxon coefficients with 95% CIs; (B) eight stable grouped-SPDE trait × predictor patterns across the four model groups with posterior intervals; (C) three BH-supported high-resolution involucral coefficients with 95% CIs.
6. **Figure 6 — ecological scale specificity and historical-placement sensitivity:** (A) endpoint × CHELSA matrix juxtaposing supported within-taxon primary (`W`) and among-taxon randomized-tree (`H`) directions; (B) the six standardized PGLS coefficients retained in all 50 randomized trees with their across-tree ranges.

## Main tables

1. **Table 1:** canonical cohorts and analysis roles.
2. **Table 2:** trait modules, endpoints and bounded measurement scope.
3. **Table 3:** multiplicity-supported exhaustive primary and high-resolution auxiliary associations.
4. **Table 4:** unique grouped SPDE-INLA associations with global BH q < 0.05.
5. **Table 5:** permutation-supported among-taxon niche contrasts.
6. **Table 6:** randomized-tree PGLS retained associations.

## Supplement tables

- **S01** `S01_canonical_cohorts.csv` — canonical cohorts and inferential roles.
- **S02** `S02_trait_scope_and_assessability.csv` — trait scope and assessability boundaries.
- **S03** `S03_primary_within_taxon_BH_supported_rows.csv` — eight BH-supported primary within-taxon component rows.
- **S04** `S04_nested_visible_variance.csv` — full nine-endpoint nested visible-variance decomposition and sampling sensitivities.
- **S05** `S05_grouped_SPDE_stable_supported_patterns.csv` — stable grouped SPDE-INLA supported directions.
- **S06** `S06_auxiliary_involucre_BH_supported_rows.csv` — BH-supported high-resolution involucre rows.
- **S07** `S07_niche_permutation_summary.csv` — niche-permutation support.
- **S08** `S08_randomized_PGLS_retained_rows.csv` — randomized-tree PGLS retained rows.
- **S09** `S09_spatial_robustness.csv` — residual spatial and broad-region omission robustness.
- **S10** `S10_wcvp_synonym_candidates.csv` — WCVP synonym-collapse sensitivity.
- **S11** `S11_precision_confounding_audit_summary.csv` — withdrawn-lability precision audit.
- **S12** `S12_submission_gates.csv` — completion gates and release boundary.

## Supplement figures

Three independent release routes are retained so descriptive figure assembly cannot silently alter frozen analyses.

- **Figures S1–S5** are the hash-verified numeric Supplement bundle generated from committed figure-data tables.
- **Figures S6–S7** are geographic sampling/measurement diagnostics generated with main Figure 2 directly from frozen observation artifacts.
- **Figure S8** is an interpretive among-taxon environmental-sorting figure generated with main Figures 1/4/5/6 from frozen result tables and the enriched-environment artifact.

Canonical products:

- **Figure S1** `Figure_S1_nested_variance_sensitivities.*` — full, one-head-per-photo and balanced-10-photo below-taxon variance sensitivities.
- **Figure S2** `Figure_S2_primary_coefficient_map.*` — complete 36-component within-taxon CHELSA coefficient map; the eight BH-supported cells are marked and hue sine/cosine remain separate.
- **Figure S3** `Figure_S3_spde_delta_waic.*` — endpoint-specific ΔWAIC across four grouped SPDE-INLA predictor sets; model-fit comparison rather than significance coding.
- **Figure S4** `Figure_S4_residual_morans_I.*` — residual Moran's I after frozen SPDE control.
- **Figure S5** `Figure_S5_leave_one_region_out_stability.*` — minimum taxon-rank stability under broad-region omission.
- **Figure S6** `Figure_S6_sampling_geography_across_filters.*` — common-scale 1° sampling density across detector-positive, coordinate-usable, ≤10 km and primary-thinned stages.
- **Figure S7** `Figure_S7_geographic_trait_assessability.*` — 2° cell-level image-trait assessability for orientation, colour, outline and the mean across all nine primary endpoints; only cells with ≥20 coordinate-usable observations are shown.
- **Figure S8** `Figure_S8_environmental_sorting.*` — observed environmental centroid distance versus Gaussian Bhattacharyya overlap for low/high trait quartiles across all nine endpoints; marker emphasis is taken from the frozen permutation-support table rather than from a new inferential test.

Figure provenance and claim boundaries are documented in `manuscript/supplement/FIGURE_MANIFEST.md`, `manuscript/FIGURE2_GEOGRAPHY_MANIFEST.md` and `manuscript/INTERPRETIVE_FIGURES_MANIFEST.md`.

## Main manuscript ↔ Supplement crosswalk

| Main item / result layer | Frozen Supplement support | Purpose |
|---|---|---|
| Figure 1 — image to multiscale analysis | S02, S12 | shows how measured numerical phenotypes enter distinct ecological questions while keeping measurement-validation gates explicit |
| Figure 2 — geographic sampling and analytical domain | S01, Figures S6–S7 | fixes cohort counts and visualizes sampling geography, filtering and trait assessability |
| Figure 3 — nested visible variance | S04, Figure S1 | provides the complete endpoint table and both pseudoreplication/equal-sampling sensitivities |
| Figure 4 — taxon-level PCA | S02 plus frozen enriched-environment source | shows actual taxon dispersion and all nine trait directions rather than loadings alone |
| Figure 5 — within-taxon environmental structure | S03, S05, S06, Figures S2–S3 | makes coefficient magnitude/uncertainty visible while separating primary OLS, grouped spatial and auxiliary families |
| Figure 6 — scale specificity + randomized-tree PGLS | S03, S08 | exposes which climate axes agree or differ between within-taxon and historical among-taxon analyses, then shows the retained PGLS effect sizes |
| Table 1 | S01 | complete cohort audit |
| Table 2 | S02 | complete trait/assessability audit |
| Table 3 | S03, S06 | complete submission-facing supported rows for the two coefficient families |
| Table 4 | S05, Figure S3 | stable SPDE directions plus full four-group model-fit comparison |
| Table 5 | S07, Figure S8 | permutation-null support plus an interpretable map of observed environmental separation/overlap |
| Table 6 | S08 | randomized-tree retained rows |
| Spatial robustness Results | S09, Figures S4–S5 | residual autocorrelation and region-omission diagnostics |
| Sampling / assessability robustness | Figures S6–S7 | shows where observations survive filtering and where image-derived traits are assessable |
| Taxonomic robustness Results | S10 | WCVP synonym-collapse sensitivity |
| Withdrawn lability audit | S11 only | statistical provenance/QA; not a biological result |
| External completion status | S12 | detector and continuous-measurement validation gates |

The old rho = -0.333 quadrant figure is statistical provenance only and is not a main or supplementary biological result. EAzami mechanism-reduction outputs are downstream context and are not inputs to any Chapter 1 figure or table.

# Chapter 1 interpretive figure manifest

This manifest defines presentation-only figures that make the frozen Chapter 1 results easier to read. None of these figures creates a new inferential family, changes a model, or changes a supported/unsupported decision.

## Main Figure 1 — image to multiscale analysis

`Figure_1_image_to_multiscale_analysis.*`

The existing open-licensed actual-photo / YOLO / deterministic-measurement Figure 1 is retained intact and augmented with one analysis-scale strip. The strip connects repeated numerical head phenotypes to (i) nested visible variance, (ii) within-taxon environmental analyses, and (iii) among-taxon niche/historical analyses. It is a design diagram, not evidence for detector or biological measurement accuracy.

Source presentation image: `manuscript/figures/supervisor_review/figure1_real_photo_yolo_pipeline.png`.

## Main Figure 4 — taxon-level trait architecture

`Figure_4_taxon_trait_architecture.*`

- Panel A shows the actual PC1–PC2 scores for the 148 taxa complete for all nine primary endpoints and overlays the nine loading directions.
- The main PC1–PC2 panel uses a central view so the dense taxon cloud and loading directions remain legible; a labelled inset preserves the complete range, including *Cirsium kawakamii* at PC1 = −15.6 rather than silently dropping the extreme point.
- Panel B shows PC3 loadings so the third independent dimension is not hidden by a PC1–PC2 plot.
- The same frozen standardized taxon-median PCA is used as in the manuscript result: PC1 = 32.9335%, PC2 = 23.1590%, PC3 = 13.2429%, cumulative PC1–PC3 = 69.3354%.

Frozen source: `strict_spatial_thinned_with_climate_topography_soil.csv` from workflow run `29306454759`, artifact `8314947270`, digest `sha256:b1252258fff85b397de8a58f322e24e90731d55f07f53ff7de16a169c5a4dfe3`.

## Main Figure 5 — environmental effect sizes

`Figure_5_environmental_effect_sizes.*`

- Panel A replaces a symbol-only support matrix with the eight frozen BH-supported primary within-taxon coefficients and 95% confidence intervals.
- Panel B shows the eight frozen stable grouped-SPDE trait × predictor patterns across all model groups in which the predictor occurs. Filled markers identify model-group coefficients with global posterior-tail BH q < 0.05; open markers retain the full model-group stability context.
- Panel C shows the three frozen BH-supported high-resolution involucre coefficients and 95% confidence intervals.

Sources:

- `manuscript/supplement/tables/S03_primary_within_taxon_BH_supported_rows.csv`;
- grouped SPDE run `29339137476`, artifact `8313574392`, digest `sha256:b877713456360144106d7531fddf95ac8c3c82fa1597d0cd3cbc4c9149a53cf5`;
- `manuscript/supplement/tables/S06_auxiliary_involucre_BH_supported_rows.csv`.

The panels retain their own coefficient definitions; their numerical scales are not pooled into one effect-size estimator.

## Main Figure 6 — scale specificity and historical sensitivity

`Figure_6_scale_specificity_and_historical_sensitivity.*`

- Panel A places supported primary within-taxon (`W`) and randomized-tree historical (`H`) climate associations in the same endpoint × CHELSA matrix. It is designed to make agreement, scale-specific support and sign reversal across ecological scales visible.
- Panel B retains the six frozen standardized PGLS coefficients supported in all 50 randomized scenario-2 trees with their 2.5–97.5% ranges.

Sources:

- `manuscript/supplement/tables/S03_primary_within_taxon_BH_supported_rows.csv`;
- `manuscript/supplement/tables/S08_randomized_PGLS_retained_rows.csv`.

The historical layer remains a placement sensitivity because only 54 of 216 atlas taxa are direct dated-backbone tips.

## Supplement Figure S8 — environmental sorting

`Figure_S8_environmental_sorting.*`

All nine primary endpoints are plotted by observed low-versus-high taxon environmental centroid distance and Gaussian Bhattacharyya overlap. The observed metrics are recomputed deterministically from the exact frozen enriched-environment table using the same complete-taxon medians, fixed environmental PCA and quartile definitions as the permutation analysis. Marker size is taken from the frozen all-taxa BH support and black rings from the direct-backbone sensitivity in `S07`; no unused colour coding is introduced and no permutation test is rerun by the figure builder.

This figure visualizes among-taxon environmental sorting only. It does not imply within-taxon response or causal adaptation.

## Reproduction

- core generator and validation logic: `analysis/build_ch1_interpretive_figures.py`;
- PGLS column-safe compatibility layer: `analysis/build_ch1_interpretive_figures_v2.py`;
- final presentation entrypoint: `analysis/build_ch1_interpretive_figures_v3.py`;
- Figure 1 augmenter: `analysis/augment_figure1_with_analysis_scales.py`;
- workflow: `.github/workflows/ch1-interpretive-figures-ci.yml`.

The CI downloads frozen source artifacts, validates the 148-taxon PCA cohort and exact PC1–PC3 explained variance, validates primary/auxiliary supported-row counts and the signs of all eight stable SPDE pairs, regenerates every SVG/PNG/PDF, records source/release SHA-256 hashes, and uploads a release artifact. Large observation-level data are not duplicated in Git.

## Verified release

The final presentation code was independently executed in GitHub Actions on head `63ce7a58889bfa3e0b351999c66db85d51fbd064`:

- workflow run: `32448807380`;
- artifact ID: `9434958333`;
- artifact name: `ch1-interpretive-figures-32448807380`;
- artifact digest: `sha256:f9365eddba862d4703e86d4fd38f03d308fd0ff035f3f7d2cf6640ce52c4ea1c`;
- 15 rendered products: five figures × SVG/PNG/PDF.

Subsequent manuscript-map/caption/contract edits do not alter the figure generator or frozen input data.

## Visual QA record

The final presentation pass was explicitly checked for the failure modes that were difficult to see in tables: duplicated Figure 1 headings, PC1 compression by the *C. kawakamii* extreme, overlapping biplot labels, misleading unused colour categories in Figure S8, and unreadable coefficient/support matrices. Figure 5 and Figure 6 required no further data/layout changes after the quantitative release was inspected.

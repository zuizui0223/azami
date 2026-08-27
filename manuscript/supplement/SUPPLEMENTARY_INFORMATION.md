# Supporting Information

## Appendix S1. Full-27 continuous-trait and full-environment analysis

### S1.1 Frozen cohort and endpoint status

The phenotype-blind spatial cohort retained one observation per source-assigned taxon by 0.25-degree cell after a positional-accuracy threshold of 10 km. It contained 46,276 observations from 259 taxa. The registered universe contained 27 endpoints: 22 were measured and five remained explicit unexecuted rows. Hue sine and cosine formed one joint circular inferential unit. Endpoint status, measurement counts and validation boundaries are given in `v2_full27_endpoint_inventory.csv`.

### S1.2 Geographic support and taxon-mean information loss

For each measured endpoint, `v2_full27_trait_geography.csv` records geographic limits and occupied 5-degree cells. `v2_full27_variance_decomposition.csv` reports the fraction of observation-level sums of squares below taxon means and 500 equal-replication resamples with two observations per eligible taxon. Missing measurements were not treated as zeros. The partitions describe visible image phenotypes rather than genetic variance or plasticity.

### S1.3 Environmental predictors

Nine predictors were frozen in six blocks: thermal (BIO1, BIO4), hydric (BIO12, BIO15), radiative/atmospheric (shortwave radiation, VPD), mechanical (wind), growing-season water input (GSP) and resource/productivity (potential NPP). Eight had complete coverage and wind had 99.989% coverage, exceeding the 98% gate. GSP was a climatic layer rather than species-specific flowering-season rain; NPP was modelled potential productivity rather than measured resources.

### S1.4 Within- and among-taxon models

Within-taxon models used taxon-demeaned standardized traits and predictors with taxon-clustered uncertainty. Among-taxon models used standardized taxon trait and environment medians with 9,999 taxon-label permutations. The primary among-taxon scope required five observations per taxon; minimum two formed a separately corrected sensitivity. Benjamini-Hochberg correction was applied once across all successful unit by predictor rows within each scope. Complete output, including unsupported and unexecuted rows, is in `v2_full27_environment_within.csv`, `v2_full27_environment_among.csv` and `v2_full27_environment_cross_scale.csv`.

### S1.5 Sampling-composition sensitivity

All 26 within-taxon and ten among-taxon globally supported rows entered a retrospective direction-only audit. The declared perturbations were equal total weight per taxon for within-taxon rows, omission of each top-ten taxon and the top two jointly, omission of each of six broad regions, and a native-only restriction. All 674 scenarios were evaluable. Sixteen within-taxon and all ten among-taxon pairs retained direction in every applicable scenario. The audit did not calculate post-selection P values and does not correct the public-image data into a probability sample.

### S1.6 Broad-space sensitivity

Globally supported rows were refitted with a second-order spherical-coordinate basis. Passing required spatial P below 0.05, sign preservation for a linear endpoint and residual Moran P of at least 0.05. Five within-taxon and two among-taxon rows passed. Joint hue rows failed because residual spatial autocorrelation remained detectable. Full results are in `spatial/v2_full27_spatial_within.csv` and `spatial/v2_full27_spatial_among.csv`.

### S1.7 Historical-placement sensitivity

The two spatially retained among-taxon rows were tested on 52 audited dated-backbone placement trees. Both retained direction and P below 0.05 on every tree. Only 54 historical-atlas taxa were direct backbone tips; most were within-genus placements. These results quantify stability to the tested placement uncertainty and do not constitute a resolved species-tree or network correction. Full model and summary tables are in `historical/v2_full27_historical_placement_models.csv` and `historical/v2_full27_historical_placement_summary.csv`.

### S1.8 Interpretation limits

The analysis is retrospective exploratory. Directional stability is not causal evidence. Visible CIELAB chroma is not anthocyanin concentration or UV reflectance; EXIF image vertical is not gravity; and image-derived involucral or surface variables are not botanical spine, hair, gland or secretion measurements. The results do not establish phenotypic plasticity, local adaptation, selection, fitness benefit, independent origin or convergence.

## Supplementary data-file index

All files below are under `analysis_outputs/v2_full27_environment_atlas_2026-08-27/` in the repository and in the submission bundle.

- `v2_full27_endpoint_inventory.csv`: all 27 registered endpoint statuses.
- `v2_full27_trait_geography.csv`: endpoint-specific geographic coverage.
- `v2_full27_variance_decomposition.csv`: raw and equal-replication information-loss partitions.
- `v2_full27_environment_within.csv`: complete within-taxon atlas.
- `v2_full27_environment_among.csv`: complete among-taxon minimum-five and minimum-two atlases.
- `v2_full27_environment_cross_scale.csv`: primary matched scale classification.
- `sampling/v2_full27_sampling_composition_scenarios.csv`: all 674 perturbation rows.
- `sampling/v2_full27_sampling_composition_summary.csv`: endpoint-level sampling stability.
- `spatial/v2_full27_spatial_within.csv` and `spatial/v2_full27_spatial_among.csv`: broad-space sensitivities.
- `historical/v2_full27_historical_placement_models.csv` and `historical/v2_full27_historical_placement_summary.csv`: 52-tree sensitivity.

## Supporting figure captions

**Figure S1. Taxon-mean information loss across all measured endpoints.** Raw observation-level and equal-replication fractions below source-assigned taxon means. Values describe visible image phenotypes.

**Figure S2. Endpoint by environmental-gradient support at within- and among-taxon scales.** Every registered inferential unit is retained; unsupported and unexecuted combinations remain visible.

**Figure S3. Sequential robustness of the two final candidates.** Global among-taxon, broad-space, sampling-composition and 52-tree placement evidence for chroma-radiation and orientation-annual-precipitation. The display summarizes robustness, not causal adaptation.

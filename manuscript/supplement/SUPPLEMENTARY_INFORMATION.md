# Supporting Information

The main manuscript contains only the inferential narrative needed to evaluate the general conclusions: continuous-trait recovery, information loss under taxon means, within-versus-among scale structure and the two final candidate patterns. This Supporting Information contains the exhaustive measurement inventory, cohort definitions, multiplicity families, sampling diagnostics, spatial diagnostics and historical-placement sensitivity needed to audit those conclusions. Main-text figures are not repeated here.

## Appendix S1. Full-27 continuous-trait and full-environment analysis

### S1.1 Frozen cohort and endpoint status

The phenotype-blind spatial cohort retained one observation per source-assigned taxon by 0.25-degree cell after a positional-accuracy threshold of 10 km. It contained 46,276 observations from 259 taxa. The registered universe contained 27 endpoints: 22 were measured and five remained explicit unexecuted rows. Hue sine and cosine formed one joint circular inferential unit. Endpoint status, module, analysis tier, validation boundary, measurement counts and taxonomic coverage are given in `v2_full27_endpoint_inventory.csv`.

The five unexecuted endpoints remain in the registry rather than disappearing from the analysis universe. Missing or unassessable measurements were never encoded as biological absence or zero. Candidate involucral, armature-like and surface-image endpoints remain image phenotypes pending botanical calibration and are not renamed as spine length, stiffness, hair density, glands or secretion.

### S1.2 Geographic support and taxon-mean information loss

For each measured endpoint, `v2_full27_trait_geography.csv` records latitude and longitude limits, observation counts, taxon counts and occupancy in 5-degree geographic cells. `v2_full27_variance_decomposition.csv` reports the observation-level fraction of sums of squares below source-assigned taxon means and 500 equal-replication resamples with two observations per eligible taxon. The equal-replication distribution, not only its median, is retained in the frozen output. These partitions describe visible image phenotypes rather than genetic variance, heritability or plasticity.

Figure S1 shows endpoint-specific measurement support rather than repeating the main information-loss figure. It is intended to expose which modules are carried by tens of thousands of photographs and which candidate image geometries are supported by substantially smaller subsets.

### S1.3 Environmental predictors and complete model families

Nine predictors were frozen in six blocks: thermal (BIO1, BIO4), hydric (BIO12, BIO15), radiative/atmospheric (shortwave radiation, VPD), mechanical (wind), growing-season water input (GSP) and resource/productivity (potential NPP). Eight had complete coverage and wind had 99.989% coverage, exceeding the 98% gate. GSP was a climatic layer rather than species-specific flowering-season rain; NPP was modelled potential productivity rather than measured resources.

Within-taxon models used taxon-demeaned standardized traits and predictors with taxon-clustered uncertainty. Among-taxon models used standardized taxon trait and environment medians with 9,999 taxon-label permutations. The primary among-taxon scope required five observations per taxon; minimum two formed a separately corrected sensitivity. Benjamini-Hochberg correction was applied once across all successful inferential-unit by predictor rows within each scope.

Complete output is retained rather than filtering to supported rows. `v2_full27_environment_within.csv`, `v2_full27_environment_among.csv` and `v2_full27_environment_cross_scale.csv` therefore include supported, unsupported, unexecuted and non-comparable combinations. This complete grid is the audit surface underlying the condensed cross-scale figure in the main manuscript.

### S1.4 Sampling-composition sensitivity

All 26 within-taxon and ten among-taxon globally supported rows entered a retrospective direction-only sampling audit. Perturbations were declared before inspecting the sensitivity outcomes: equal total weight per taxon for within-taxon rows, omission of each top-ten taxon separately, joint omission of the two most observed taxa, omission of each of six broad regions, and an explicitly native-only restriction. All 674 scenarios were evaluable.

Sixteen of 26 selected within-taxon pairs and all ten selected among-taxon pairs retained direction in every applicable scenario. The audit did not calculate post-selection P values and does not convert the public-image sample into a probability sample. `sampling/v2_full27_sampling_composition_scenarios.csv` preserves every scenario-level coefficient and `sampling/v2_full27_sampling_composition_summary.csv` records the minimum retained effect-magnitude ratios and direction-stability classifications.

Figure S2 shows every selected row, separating stable and direction-unstable rows and exposing the weakest retained effect magnitude across the declared perturbations. This is intentionally more diagnostic than the two-candidate summary in the main manuscript.

### S1.5 Broad-space sensitivity and residual spatial structure

Every globally FDR-supported row was refitted with a second-order spherical-coordinate basis. Passing required a spatial permutation P below 0.05, sign preservation for a linear endpoint and residual eight-nearest-neighbour Moran P of at least 0.05. Five within-taxon and two among-taxon rows passed all three requirements. Joint hue rows failed the full gate because residual spatial autocorrelation remained detectable even when the spatial-basis association itself was non-zero.

Full results are in `spatial/v2_full27_spatial_within.csv` and `spatial/v2_full27_spatial_among.csv`. Figure S3 displays the complete broad-space diagnostic set, including rows that failed because the spatial effect weakened, reversed or retained residual autocorrelation. The main manuscript shows only the sequentially retained candidates.

### S1.6 Historical-placement sensitivity

The two spatially retained among-taxon rows were tested on 52 audited dated-backbone placement trees: two deterministic scenarios and 50 randomized within-genus placement trees. Both retained direction and P below 0.05 on every tree. Only 54 historical-atlas taxa were direct backbone tips; most were within-genus placements. The analysis is therefore a sensitivity to the tested placement uncertainty, not a resolved species-tree or network correction.

All tree-level coefficients, P values and Pagel-lambda estimates are retained in `historical/v2_full27_historical_placement_models.csv`; endpoint-level ranges are in `historical/v2_full27_historical_placement_summary.csv`. These tables, rather than a duplicate main figure, provide the audit trail for the 52-tree statement.

### S1.7 Whole-capitulum secondary synthesis

The complete-18 secondary analysis was performed only after endpoint-level inference was frozen. It evaluates whether taxon morphospace and within/among association matrices show partial organization among registered measurement modules. It does not define the primary endpoint set, alter multiplicity families or replace endpoint-level inference. Full PCA scores/loadings, module tests and matrix summaries are retained in the frozen multilevel/process artifact and indexed in the submission bundle.

### S1.8 Interpretation limits and independent validation gates

The analysis is retrospective exploratory. Directional stability is not causal evidence. Visible CIELAB chroma is not anthocyanin concentration or UV reflectance; EXIF image vertical is not gravity; and image-derived involucral or surface variables are not botanical spine, hair, gland or secretion measurements. The results do not establish phenotypic plasticity, local adaptation, selection, fitness benefit, independent origin, adaptive radiation or convergence.

The submission remains on scientific HOLD for independent detector-box validation, independent primary-trait reference measurements, botanical calibration of candidate fine geometry, developmental-stage control, paired colour controls and repeat-photo remeasurement. These gates are external validation tasks, not invitations to add further post hoc environmental screens.

## Supplementary data-file index

All files below are under `analysis_outputs/v2_full27_environment_atlas_2026-08-27/` in the repository and in the submission bundle.

- `v2_full27_endpoint_inventory.csv`: all 27 registered endpoint statuses, validation boundaries and coverage.
- `v2_full27_trait_geography.csv`: endpoint-specific geographic support.
- `v2_full27_variance_decomposition.csv`: raw and equal-replication information-loss partitions.
- `v2_full27_environment_within.csv`: complete within-taxon atlas.
- `v2_full27_environment_among.csv`: complete among-taxon minimum-five and minimum-two atlases.
- `v2_full27_environment_cross_scale.csv`: primary matched scale classification.
- `sampling/v2_full27_sampling_composition_scenarios.csv`: all 674 perturbation rows.
- `sampling/v2_full27_sampling_composition_summary.csv`: selected-row sampling stability and effect-magnitude ratios.
- `spatial/v2_full27_spatial_within.csv`: complete within-taxon broad-space diagnostics for globally supported rows.
- `spatial/v2_full27_spatial_among.csv`: complete among-taxon broad-space diagnostics for globally supported rows.
- `historical/v2_full27_historical_placement_models.csv`: all tree-level historical-placement models.
- `historical/v2_full27_historical_placement_summary.csv`: 52-tree endpoint summaries.

## Supporting figure captions

**Figure S1. Measurement support across the registered continuous-trait universe.** Observation and taxon coverage for every measured endpoint, with explicit unexecuted endpoints retained in the inventory. This figure diagnoses data support and does not repeat the main taxon-mean information-loss result.

**Figure S2. Sampling-composition sensitivity of all globally supported rows entering the audit.** Minimum retained effect-magnitude ratio across declared sampling perturbations for each selected within- and among-taxon association, with direction-instability flagged. This figure shows the full selected-row audit rather than only the two final candidates.

**Figure S3. Broad-space and residual-spatial diagnostic surface for all globally supported rows entering spatial sensitivity.** Spatial permutation support and residual Moran diagnostics are shown for within- and among-taxon associations, with full-gate passage distinguished from failure. This figure exposes why many initially supported rows do not enter the final candidate set.

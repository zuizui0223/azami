# Chapter 1 submission manuscript

This is the canonical submission-facing entry point for Chapter 1. The paper is framed as a **global quantitative geography of continuous thistle capitulum traits**: public photographs are converted into repeated continuous observations so that each trait can be paired with predeclared environmental gradients and located at head, photograph, within-taxon, among-taxon and sampled-global scales. The expanded whole-capitulum analysis is a secondary synthesis, not a return to one composite morphology.

The narrative decision and Azami-to-EAzami handoff are fixed in [`AZAMI_CH1_SPATIAL_SCALE_CANONICAL_2026-08-27.md`](AZAMI_CH1_SPATIAL_SCALE_CANONICAL_2026-08-27.md).

The former lability-axis story is not a headline result. Mechanistic explanations are downstream hypotheses rather than the primary inferential target.

**Current status: SUBMISSION HOLD FOR EXTERNAL VALIDATION, NOT FOR COMPUTATIONAL REANALYSIS.** The final artifact-backed GEB v2 continuous-trait workflow, the all-27 / all-nine endpoint atlas, the matched within/among comparison, the retrospective sampling-composition audit and the sequential spatial/historical candidate gates are complete. Whole-capitulum and core-four environmental-sufficiency analyses remain secondary synthesis. Independent detector/trait validation and the other genuinely external gates listed below remain open.

## Preferred title

**Global continuous phenomics reveals trait- and scale-specific environmental geography of thistle capitula**

## Manuscript order

1. [`00_title_abstract.md`](00_title_abstract.md)
2. [`01_introduction.md`](01_introduction.md)
3. [`02_methods.md`](02_methods.md)
4. [`03_results.md`](03_results.md)
5. [`04_discussion.md`](04_discussion.md)
6. [`05_conclusion_and_declarations.md`](05_conclusion_and_declarations.md)
7. [`06_references.md`](06_references.md)
8. [`07_data_code_availability.md`](07_data_code_availability.md)

`06_references.md` is the submission-wide reference list. Older Introduction-only literature files remain drafting provenance and are not the reference source for the final manuscript.

## Final artifact-backed executions

The completed image-measurement and continuous-trait workflow is recorded in [`results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](results/GEB_V2_FINAL_RESULTS_2026-08-26.md):

- run `32975451732` on head `f4a6fd5e01a2befd4f49174984a99e53856c2330`;
- artifact `9612943217`;
- digest `sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`;
- 27 registered continuous endpoints, 22 measured;
- 19 inferential endpoints, 18 measured: 9 primary + 9 candidate;
- `visible_floret_fraction` is the only unexecuted inferential endpoint.

The final multilevel capitulum-space and process-environment extension is recorded in [`GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md`](GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md):

- successful run `33035785120` on analysis head `227c0e7b8c338894806785b8545c7c77c8724de1`;
- artifact `9632715852`;
- digest `sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6`;
- 10,000 permutations per inferential test and 1,000 taxon bootstraps for multivariate geometry;
- 1,874 observations from 124 taxa complete across all 18 measured inferential endpoints;
- complete CHELSA process-variable coverage for shortwave radiation, VPD, wind, growing-season precipitation and potential NPP.

The numerical results remain tied to these immutable artifacts and analysis SHAs. Later manuscript-only commits do not redefine the frozen estimands.

The canonical all-27 / all-nine reanalysis is recorded in [`results/GEB_V2_FULL27_FULL_ENVIRONMENT_RESULTS_2026-08-27.md`](results/GEB_V2_FULL27_FULL_ENVIRONMENT_RESULTS_2026-08-27.md). It retained all 27 registered endpoints, joined the full 46,276-observation cohort to nine predictors, and passed 24/24 atlas checks, 7/7 sampling-composition checks and 11/11 spatial/historical checks. These local outputs are not yet a new remote immutable artifact.

## Current headline

The study moves from:

> species means / categorical trait states

through:

> public photographs → detected capitula → explicitly defined continuous measurements

to:

> repeated phenotype distributions → trait-specific present-day geography → predeclared environmental gradients → within/among scale comparison → sampling-composition and spatial robustness → bounded historical handoff.

Across all 22 measured v2 endpoints, 0.814–0.986 of raw observation-level variation occurs below source-assigned taxon means. Equalizing replication to two observations per eligible taxon retains endpoint medians of 0.287–0.585. The expanded taxon morphospace contains 127 complete taxa × 18 measured inferential endpoints; PC1, PC2 and PC3 explain 18.5%, 12.0% and 11.8%, respectively. The first three axes therefore capture only 42.3%, arguing against one universal capitulum syndrome.

The canonical full-environment comparison shows that phenotype–environment structure changes with biological scale. Across 234 endpoint × predictor rows, seven are supported at both scales, 18 within taxa only, three among taxa only, 152 at neither scale and 54 are not comparable. Orientation aligns with annual mean temperature within taxa but annual precipitation among taxa, and most supported involucral or surface-proxy rows are within-taxon only.

The sample is geographically extensive but uneven: Europe plus North America contribute 92.26% of observations, the Northern Hemisphere 94.91%, and *Cirsium vulgare* plus *C. arvense* 54.03%. Across 674 retrospective sampling-composition scenarios, 16 of 26 selected within-taxon pairs and all ten selected among-taxon pairs retain direction throughout. Ten within-taxon pairs reverse at least once and remain explicitly annotated as sampling sensitive.

Five within-taxon rows pass the new broad-space gate, but projection roughness-radiation reverses under equal taxon weighting and joint omission of the two most observed taxa; the other four retain direction in every sampling-composition scenario. Only two among-taxon rows—lower chroma with higher shortwave radiation and a larger signed head-axis angle relative to EXIF-image vertical with higher annual precipitation—also pass the broad-space gate and all 52 historical-placement trees. Both retain direction throughout the sampling audit. They are adaptive-pattern candidates under current controls, not demonstrated adaptations, mechanisms or globally representative effects. Earlier A/B/C endpoint grades remain provenance and were not used to select or hide rows in this lane.

The candidate-function defence is explicitly bounded. Radiation-induced petal anthocyanin production in other species makes pigment regulation a plausible hypothesis for chroma-radiation, but photographed CIELAB chroma is not an anthocyanin assay and its negative coefficient does not identify the mediator. The orientation scale runs from 0 degrees upward through 90 degrees horizontal to 180 degrees downward in EXIF-oriented image coordinates. Its positive precipitation coefficient is therefore consistent with more downward presentation, and experimental work in a nodding Asteraceae supports rain shielding as a plausible function; image vertical is not gravity and Azami measured neither rain interception nor fitness. These hypotheses enter EAzami-I as tests, not Chapter 1 conclusions.

The environmental-sufficiency test gives an equally scale-dependent result. Adding shortwave radiation, VPD, wind, growing-season precipitation and potential NPP beyond BIO1/BIO4/BIO12/BIO15 provides no supported incremental whole-capitulum information within taxa. Among taxa, the same extension is supported at both replication thresholds; growing-season precipitation is the only block-specific increment supported in both scopes. The four variables are therefore an adequate within-taxon baseline but do not exhaust among-taxon environmental structure.

The whole-organ analysis provides a secondary check on whether the decomposed endpoints can be collapsed again. They cannot: the first three PCs explain only 42.3%, registered-module association exceeds between-module association at both scales, and the two 17-unit matrices are only partly aligned. This rejects one universal syndrome while retaining partial measurement-module organization; it does not replace the endpoint-level geography.

The central biological contribution is **a quantitative map showing that continuous capitulum traits occupy different present-day spatial scales and align with different environmental gradients**. The central methodological contribution is a scalable route from heterogeneous public photographs to repeated, coordinate-bearing continuous trait distributions and complete trait × gradient × scale status grids.

## Current main-result structure

1. **Continuous capitulum phenotypes can be recovered from global public imagery.**
2. **Taxon means compress substantial variation across all 22 measured endpoints, including under equal replication.**
3. **Coordinates turn repeated phenotype values into endpoint-specific geographic support.**
4. **All nine predictors in six hypothesis blocks produce a complete, trait-specific status atlas.**
5. **Matched endpoint–environment support differs between within- and among-taxon scales.**
6. **Sampling perturbations retain 16/26 selected within-taxon and 10/10 selected among-taxon directions; ten within rows are explicitly sensitive.**
7. **Five within-taxon rows pass the v2-native broad-space gate, but only four also retain direction across all sampling perturbations.**
8. **Two among-taxon rows pass both broad-space and all historical-placement gates and retain direction across all sampling perturbations.**
9. **The four-variable CHELSA core is an adequate multivariate baseline within taxa but does not exhaust among-taxon structure.**
10. **The expanded 18-endpoint morphospace rejects one universal capitulum syndrome.**
11. **Registered measurement modules are partially organized, but within- and among-taxon geometry is only partly shared.**
12. **Azami freezes present trait states and scale classes for EAzami-I function, II repeated history, III origin-class and IV convergence gates.**

## Historical evidence tiers

Earlier frozen analyses use the pre-result interpretation policy in [`../analysis/ch1/evidence_tiering_policy.json`](../analysis/ch1/evidence_tiering_policy.json). The canonical full-27 lane retains those labels only as provenance; it does not use them for entry, multiplicity or row suppression.

- **A — robust:** supported globally and consistent across the relevant sensitivity analyses;
- **B — supported but sensitive:** supported in the global primary/spatial analysis but weakened in a restricted-domain sensitivity while retaining direction, or limited by calibration appropriate to the endpoint;
- **C — exploratory continuous signal:** coherent image-derived candidate association requiring botanical validation even when image-quality robust;
- **image-sensitive C:** candidate signal with explicit sign instability across plausible image-resolution strata;
- **validation-only:** surface-image signals that cannot yet be given a botanical ecological interpretation.

The new whole-capitulum and process-environment results form a separate multivariate pattern layer. They do not promote or demote the frozen endpoint-level A/B/C evidence atlas. Stand-alone environmental-block BH decisions, nested incremental support and cross-scale coefficient geometry must remain distinguishable.

## Figure priorities

The main figures should make trait geography, environmental hypotheses and scale structure visible before the whole-capitulum synthesis:

1. **image → quantitative phenotype → ecological scale:** actual photographs, detection and continuous measurement route, including examples from the expanded contract;
2. **geographic sampling and analytical domain:** global sampling, filtering streams and occupied environmental space;
3. **nested visible variation:** taxon → photograph → head decomposition across established continuous endpoints;
4. **trait-specific present-state atlas:** endpoint coverage, continuous morphospace and geographic distributions, with PCA used as a summary rather than one response axis;
5. **trait × gradient × scale:** complete endpoint status beside matched within- and among-taxon support classes;
6. **boundary and EAzami handoff:** broad-space passes, historical-placement candidates, core-four sufficiency and the state bundle passed to the EAzami I–IV evidence ladder.

The full whole-capitulum association matrices, historical-placement detail and sampling/image-quality diagnostics remain Supporting Information unless a compact synthesis panel can be retained without displacing the trait-gradient result.

## Scientific control files

- [`COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](COHORT_FLOW_AND_ANALYSIS_LEDGER.md) fixes cohort names, counts and permitted analyses.
- [`final_claims.json`](final_claims.json) contains the machine-readable numerical claim registry for the established manuscript layer.
- [`results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](results/GEB_V2_FINAL_RESULTS_2026-08-26.md) records final expanded-run provenance, endpoint coverage, morphospace and evidence grades.
- [`GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md`](GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md) records the final multilevel phenotype-space and environmental-sufficiency results.
- [`CAPITULUM_FUNCTIONAL_SPACE_HYPOTHESES_2026-08-27.md`](CAPITULUM_FUNCTIONAL_SPACE_HYPOTHESES_2026-08-27.md) freezes the whole-capitulum organization hypotheses and claim boundaries.
- [`CAPITULUM_ENVIRONMENT_CORE_SUFFICIENCY_HYPOTHESIS_2026-08-27.md`](CAPITULUM_ENVIRONMENT_CORE_SUFFICIENCY_HYPOTHESIS_2026-08-27.md) freezes the nested core-four test and appends its immutable artifact-backed outcome.
- [`AZAMI_CH1_SPATIAL_SCALE_CANONICAL_2026-08-27.md`](AZAMI_CH1_SPATIAL_SCALE_CANONICAL_2026-08-27.md) fixes the PR #71 spatial-scale mainline and the timing/repetition/convergence EAzami handoff.
- [`../analysis/ch1/evidence_tiering_policy.json`](../analysis/ch1/evidence_tiering_policy.json) fixes how global results and sensitivity analyses are interpreted.
- [`FIGURE_TABLE_MAP.md`](FIGURE_TABLE_MAP.md) fixes the main/Supplement figure–table crosswalk.
- [`MAIN_FIGURE_CAPTIONS.md`](MAIN_FIGURE_CAPTIONS.md) contains authoritative main-figure captions.
- [`FIGURE2_GEOGRAPHY_MANIFEST.md`](FIGURE2_GEOGRAPHY_MANIFEST.md) records geographic figure provenance.
- [`INTERPRETIVE_FIGURES_MANIFEST.md`](INTERPRETIVE_FIGURES_MANIFEST.md) records sources and claim boundaries for interpretive figures.
- [`results/nested_visible_variance_summary.csv`](results/nested_visible_variance_summary.csv) records the revised image hierarchy.
- [`results/NICHE_PERMUTATION_RESULTS_2026-08-07.md`](results/NICHE_PERMUTATION_RESULTS_2026-08-07.md) records the established environmental-sorting null.
- [`results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md`](results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md) records residual Moran and region-omission diagnostics.
- [`results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md`](results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md) records synonym-collapse sensitivity.
- [`06_references.md`](06_references.md) is the submission-wide bibliography.
- [`07_data_code_availability.md`](07_data_code_availability.md) contains the blinded-review repository and durable-release statement.
- [`../analysis/ch1/pipeline.json`](../analysis/ch1/pipeline.json) maps canonical executable stages.
- [`EXTERNAL_COMPLETION_GATES.md`](EXTERNAL_COMPLETION_GATES.md) records completed versus genuinely external validation gates.
- [`EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md`](EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md) fixes the Chapter 1 → EAzami boundary.

## Current submission blockers

Submission remains on hold for scientific validation that cannot be supplied by the completed computational reruns:

1. adjudicated human boxes for independent detector precision/recall;
2. independent reference measurements for orientation, colour and outline;
3. independent botanical calibration for expanded involucral/armature image phenotypes before botanical naming or functional promotion;
4. developmental-stage labels or an anthesis-restricted sensitivity;
5. paired flower-versus-background colour negative controls;
6. repeat-photo trait remeasurement and variance partitioning.

`visible_floret_fraction` remains unexecuted and is not used for a floral-display conclusion. The four unexecuted colour-composition endpoints are descriptive and do not reduce the 18-endpoint inferential analysis to a failed run.

Taxonomic robustness, residual spatial autocorrelation, broad-region omission, raw-calendar and hemisphere-aware cyclic collection timing, dominant-taxon omission, native-range restriction, the full-27 sampling-composition audit, environmental-niche permutation, candidate image-quality analysis, matched within/among comparison, complete-18 module organization and process-environment sufficiency are complete for the frozen operational-unit framework. Administrative metadata, licence selection, nomenclatural notes and final durable DOI release remain human/finalization tasks.

## Claim boundary

The final journal package must not:

- use the withdrawn negative lability relation or median-split quadrants as biological conclusions;
- use *climate tracking* or *environmental responsiveness* for cross-sectional spatial associations;
- equate visible image variance with genetic variance or plasticity;
- call registered image modules validated functional or genetic modules;
- interpret hue sine/cosine separately as named flower colours;
- treat candidate image-geometry proxies as direct botanical spine/phyllary measurements before calibration;
- call surface/specular image signals hair density, gland density, mucilage or secretion without targeted validation;
- infer floral-display structure from the unexecuted `visible_floret_fraction` endpoint;
- claim that growing-season precipitation, VPD, radiation or productivity is a demonstrated causal process;
- claim resolved phylogenetic correction from the grafted historical layer;
- convert macro-scale associations into adaptation, pollinator causation, antagonist defence or evolutionary-rate claims.

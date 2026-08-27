# Chapter 1 submission manuscript

This is the canonical submission-facing entry point for Chapter 1. The paper is framed as a **global quantitative phenotypic geography of thistle capitula**: public photographs are converted into repeated continuous trait observations, allowing visible diversity below taxon means, among-taxon morphospace and environmental organization at multiple biological scales to be examined without reducing morphology to species means or coarse categories.

The former lability-axis story is not a headline result. Mechanistic explanations are downstream hypotheses rather than the primary inferential target.

**Current status: SUBMISSION HOLD FOR EXTERNAL VALIDATION, NOT FOR COMPUTATIONAL REANALYSIS.** The final artifact-backed GEB v2 continuous-trait workflow and the multilevel within/among/process-environment extension are complete. Raw-calendar and hemisphere-aware seasonal, dominant-taxon, native-range, spatial, taxonomic, niche-permutation, candidate image-quality and core-four environmental-sufficiency analyses are complete. Independent detector/trait validation and the other genuinely external gates listed below remain open.

## Preferred title

**Global citizen-science images reveal multidimensional and scale-dependent environmental organization of thistle capitula**

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

## Current headline

The study moves from:

> species means / categorical trait states

through:

> public photographs → detected capitula → explicitly defined continuous measurements

to:

> repeated phenotype distributions → a multidimensional and partially organized capitulum space → trait- and scale-specific environmental alignment.

Across the nine established primary continuous endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means, although biological and imaging variance are not yet fully separable. The expanded taxon morphospace contains 127 complete taxa × 18 measured endpoints; PC1, PC2 and PC3 explain 18.5%, 12.0% and 11.8%, respectively. The first three axes therefore capture only 42.3%, arguing against one universal capitulum syndrome.

The whole-organ analysis adds structure to that multidimensionality. In the main complete-18 scope, registered-module association is stronger than between-module association within taxa (contrast 0.1645) and among taxa (0.0885). The within- and among-taxon 17-unit association matrices are positively but incompletely aligned (Spearman 0.3663). The allowed interpretation is **partial phenotypic / measurement-module organization**, not validated functional or genetic modularity.

In the established endpoint-level analysis, orientation–BIO1 and aspect ratio–BIO4 are the two A-grade primary associations. Chroma–BIO12 and aspect ratio–BIO12 retain their global directions but are range-sensitive under native-only restriction; circular colour rows remain calibration-dependent. The expanded candidate layer adds involucre length/width–BIO15 and apical taper–BIO12 as quality-robust C-grade signals, while `bract_projection_p95`–BIO15 remains image-sensitive. All three candidate rows are within-taxon only in the matched cross-scale analysis.

The matched 17-unit comparison shows that phenotype–environment structure changes with biological scale. Across 68 endpoint × four-predictor rows, three are supported at both scales, eight within taxa only, one among taxa only and 56 at neither scale. Orientation aligns with annual mean temperature within taxa but annual precipitation among taxa.

The environmental-sufficiency test gives an equally scale-dependent result. Adding shortwave radiation, VPD, wind, growing-season precipitation and potential NPP beyond BIO1/BIO4/BIO12/BIO15 provides no supported incremental whole-capitulum information within taxa. Among taxa, the same extension is supported at both replication thresholds; growing-season precipitation is the only block-specific increment supported in both scopes. The four variables are therefore an adequate within-taxon baseline but do not exhaust among-taxon environmental structure.

The central biological contribution is **a quantitative map of a multidimensional, partially organized reproductive phenotype whose environmental alignment rotates across biological scales**. The central methodological contribution is a scalable route from heterogeneous public photographs to repeated continuous trait distributions, multilevel phenotype geometry and evidence-graded environmental tests.

## Current main-result structure

1. **Continuous capitulum phenotypes can be recovered from global public imagery.**
2. **Taxon means conceal substantial repeated visible variation.**
3. **The below-taxon pattern survives one-head-per-photo and equal-replication controls.**
4. **The expanded 18-endpoint morphospace remains strongly multidimensional.**
5. **Registered measurement modules show detectable internal organization within and among taxa.**
6. **Within- and among-taxon phenotype geometry is only partly shared.**
7. **Within-taxon environmental associations are trait specific.**
8. **Expanded involucral geometry separates two quality-robust exploratory signals from one resolution-sensitive signal.**
9. **Matched endpoint–environment support differs between within- and among-taxon scales.**
10. **The four-variable CHELSA core is adequate for the current within-taxon whole-capitulum estimand.**
11. **A broader process representation is required among taxa, with growing-season precipitation the most stable increment.**
12. **Historical/PGLS results remain Supporting Information placement sensitivities rather than resolved phylogenetic correction.**

## Evidence tiers

The submission uses the pre-result interpretation policy in [`../analysis/ch1/evidence_tiering_policy.json`](../analysis/ch1/evidence_tiering_policy.json):

- **A — robust:** supported globally and consistent across the relevant sensitivity analyses;
- **B — supported but sensitive:** supported in the global primary/spatial analysis but weakened in a restricted-domain sensitivity while retaining direction, or limited by calibration appropriate to the endpoint;
- **C — exploratory continuous signal:** coherent image-derived candidate association requiring botanical validation even when image-quality robust;
- **image-sensitive C:** candidate signal with explicit sign instability across plausible image-resolution strata;
- **validation-only:** surface-image signals that cannot yet be given a botanical ecological interpretation.

The new whole-capitulum and process-environment results form a separate multivariate pattern layer. They do not promote or demote the frozen endpoint-level A/B/C evidence atlas. Stand-alone environmental-block BH decisions, nested incremental support and cross-scale coefficient geometry must remain distinguishable.

## Figure priorities

The main figures should make phenotype organization and scale structure visible before presenting detailed sensitivity analyses:

1. **image → quantitative phenotype → ecological scale:** actual photographs, detection and continuous measurement route, including examples from the expanded contract;
2. **geographic sampling and analytical domain:** global sampling, filtering streams and occupied environmental space;
3. **nested visible variation:** taxon → photograph → head decomposition across established continuous endpoints;
4. **expanded multivariate capitulum space:** 18 inferential endpoints, with the 9-endpoint primary PCA retained as the frozen comparison;
5. **multilevel phenotype organization:** within- and among-taxon 17-unit association matrices, registered-module contrast and cross-scale similarity;
6. **environmental evidence and sufficiency:** endpoint-level A/B/C atlas plus the nested core-four versus process-extension result, emphasizing the within/among contrast rather than a list of unadjusted block P values.

Historical-placement sensitivity and detailed sampling/image-quality diagnostics remain Supporting Information where needed.

## Scientific control files

- [`COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](COHORT_FLOW_AND_ANALYSIS_LEDGER.md) fixes cohort names, counts and permitted analyses.
- [`final_claims.json`](final_claims.json) contains the machine-readable numerical claim registry for the established manuscript layer.
- [`results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](results/GEB_V2_FINAL_RESULTS_2026-08-26.md) records final expanded-run provenance, endpoint coverage, morphospace and evidence grades.
- [`GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md`](GEB_V2_CROSS_SCALE_FINAL_RESULTS_2026-08-27.md) records the final multilevel phenotype-space and environmental-sufficiency results.
- [`CAPITULUM_FUNCTIONAL_SPACE_HYPOTHESES_2026-08-27.md`](CAPITULUM_FUNCTIONAL_SPACE_HYPOTHESES_2026-08-27.md) freezes the whole-capitulum organization hypotheses and claim boundaries.
- [`CAPITULUM_ENVIRONMENT_CORE_SUFFICIENCY_HYPOTHESIS_2026-08-27.md`](CAPITULUM_ENVIRONMENT_CORE_SUFFICIENCY_HYPOTHESIS_2026-08-27.md) freezes the nested core-four test and appends its immutable artifact-backed outcome.
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

Taxonomic robustness, residual spatial autocorrelation, broad-region omission, raw-calendar and hemisphere-aware cyclic collection timing, dominant-taxon omission, native-range restriction, environmental-niche permutation, candidate image-quality analysis, matched within/among comparison, complete-18 module organization and process-environment sufficiency are complete for the frozen operational-unit framework. Administrative metadata, licence selection, nomenclatural notes and final durable DOI release remain human/finalization tasks.

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

# Chapter 1 submission manuscript

This is the canonical submission-facing entry point for Chapter 1. The paper is framed as a **global quantitative phenotypic geography of thistle capitula**: public photographs are converted into repeated continuous trait observations, allowing spatial diversity within and among taxa and its association with environment to be examined without reducing morphology to species means or coarse categories.

The former lability-axis story is not a headline result. Mechanistic explanations are downstream hypotheses rather than the primary inferential target.

**Current status: SUBMISSION HOLD FOR EXTERNAL VALIDATION, NOT FOR COMPUTATIONAL REANALYSIS.** The final artifact-backed GEB v2 continuous-trait workflow is complete. Raw-calendar and hemisphere-aware seasonal, dominant-taxon, native-range, spatial, taxonomic, niche-permutation and candidate image-quality audits are complete. Independent detector/trait validation and the other genuinely external gates listed below remain open.

## Preferred title

**Global citizen-science images reveal continuous within-taxon variation and environmental structure in thistle capitulum traits**

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

## Final GEB v2 execution

The completed full workflow is recorded in [`results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](results/GEB_V2_FINAL_RESULTS_2026-08-26.md).

- run `32975451732` on head `f4a6fd5e01a2befd4f49174984a99e53856c2330`;
- artifact `9612943217`;
- digest `sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`;
- 27 registered continuous endpoints, 22 measured;
- 19 inferential endpoints, 18 measured: 9 primary + 9 candidate;
- `visible_floret_fraction` is the only unexecuted inferential endpoint.

The computational v2 integration gate is therefore closed. New commits after the recorded run may rerun ordinary PR checks, but the numerical v2 result above remains tied to the immutable artifact and head SHA.

## Current headline

The study moves from:

> species means / categorical trait states

through:

> public photographs → detected capitula → explicitly defined continuous measurements

to:

> repeated within-taxon phenotype distributions → multidimensional global phenotypic geography → trait- and scale-specific environmental structure.

Across the nine established primary continuous endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means, although biological and imaging variance are not yet fully separable. The expanded all-inferential morphospace contains 127 complete taxa × 18 measured endpoints; PC1, PC2 and PC3 explain 18.5%, 12.0% and 11.8%, respectively. The v2 result therefore strengthens the multidimensional-phenotype interpretation rather than producing one universal capitulum syndrome.

In the established global within-taxon analysis, significant environmental structure involves orientation, visible colour and outline aspect ratio. Grouped SPDE-INLA models concentrate support in orientation and visible colour, while among-taxon permutation tests likewise identify the strongest environmental sorting in orientation and colour. Orientation–BIO1 and aspect ratio–BIO4 are the two A-grade primary associations. Chroma–BIO12 and aspect ratio–BIO12 retain their global directions but are range-sensitive under native-only restriction; circular colour rows remain calibration-dependent.

The expanded candidate layer adds three registry-wide FDR signals without promoting them to primary evidence. Involucre length/width ratio–BIO15 (β = +0.054923, q = 0.000282) and apical taper–BIO12 (β = −0.049721, q = 0.000821) retain direction after resolution/sharpness adjustment and across successful fixed resolution strata, so they are `C_exploratory_quality_robust`. `bract_projection_p95`–BIO15 (β = −0.033673, q = 0.014377) shows resolution-stratum sign instability and is `C_exploratory_image_sensitive`. Botanical calibration remains required for all candidate geometry.

The central biological contribution is therefore **the quantitative mapping of multidimensional capitulum phenotype diversity and its environmental structure**, not the demonstration of a single adaptive mechanism. The central methodological contribution is a scalable route from heterogeneous public photographs to repeated continuous trait distributions under an explicit evidence hierarchy.

## Current main-result structure

1. **Continuous phenotypes can be recovered from global public imagery.**
2. **Taxon means conceal substantial repeated visible variation.**
3. **The below-taxon pattern survives one-head-per-photo and equal-replication controls.**
4. **The expanded 18-endpoint inferential morphospace remains multidimensional rather than collapsing onto one morphology axis.**
5. **Within-taxon environmental associations are trait specific, involving orientation, visible colour and selected outline dimensions.**
6. **Grouped spatial models concentrate the most consistent established support in orientation and visible colour.**
7. **Expanded involucral geometry separates two quality-robust exploratory signals from one resolution-sensitive signal.**
8. **Among-taxon permutation tests concentrate environmental sorting mainly in orientation and visible colour.**
9. **Historical/PGLS results are Supporting Information placement sensitivities rather than resolved phylogenetic correction.**
10. **Seasonal, dominant-taxon, native-range and image-quality analyses grade confidence rather than redefining the global primary question.**

## Evidence tiers

The submission uses the pre-result interpretation policy in [`../analysis/ch1/evidence_tiering_policy.json`](../analysis/ch1/evidence_tiering_policy.json):

- **A — robust:** supported globally and consistent across the relevant sensitivity analyses;
- **B — supported but sensitive:** supported in the global primary/spatial analysis but weakened in a restricted-domain sensitivity while retaining direction, or limited by calibration appropriate to the endpoint;
- **C — exploratory continuous signal:** coherent image-derived candidate association requiring botanical validation even when image-quality robust;
- **image-sensitive C:** candidate signal with explicit sign instability across plausible image-resolution strata;
- **validation-only:** surface-image signals that cannot yet be given a botanical ecological interpretation.

Sensitivity analyses therefore rank evidence strength. A result is not automatically erased solely because a restricted-domain BH value crosses 0.05, and a candidate is not automatically promoted because a quality-adjusted q value falls below 0.05.

## Figure priorities

The main figures should make the phenotype landscape and scale structure visible before presenting sensitivity details:

1. **image → quantitative phenotype → ecological scale:** actual photographs, detection and continuous measurement route, including examples from the expanded contract;
2. **geographic sampling and analytical domain:** global sampling, filtering streams and occupied environmental space;
3. **nested visible variation:** taxon → photograph → head decomposition across established continuous endpoints;
4. **expanded multivariate capitulum morphospace:** 18 inferential endpoints, with the 9-endpoint primary PCA retained as the frozen primary comparison rather than the sole morphospace;
5. **environmental evidence atlas:** established primary effect sizes plus the three final candidate signals, visually annotated by A/B/C grade and quality/range context rather than reduced to significant/non-significant cells.

Historical-placement sensitivity and detailed sampling/image-quality diagnostics remain Supporting Information where needed.

## Scientific control files

- [`COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](COHORT_FLOW_AND_ANALYSIS_LEDGER.md) fixes cohort names, counts and permitted analyses.
- [`final_claims.json`](final_claims.json) contains the machine-readable numerical claim registry for the established manuscript layer.
- [`results/GEB_V2_FINAL_RESULTS_2026-08-26.md`](results/GEB_V2_FINAL_RESULTS_2026-08-26.md) records final expanded-run provenance, endpoint coverage, morphospace and evidence grades.
- [`../analysis/ch1/evidence_tiering_policy.json`](../analysis/ch1/evidence_tiering_policy.json) fixes how global results and sensitivity analyses are interpreted.
- [`FIGURE_TABLE_MAP.md`](FIGURE_TABLE_MAP.md) fixes the main/Supplement figure–table crosswalk.
- [`MAIN_FIGURE_CAPTIONS.md`](MAIN_FIGURE_CAPTIONS.md) contains authoritative main-figure captions.
- [`FIGURE2_GEOGRAPHY_MANIFEST.md`](FIGURE2_GEOGRAPHY_MANIFEST.md) records geographic figure provenance.
- [`INTERPRETIVE_FIGURES_MANIFEST.md`](INTERPRETIVE_FIGURES_MANIFEST.md) records sources and claim boundaries for interpretive figures.
- [`results/nested_visible_variance_summary.csv`](results/nested_visible_variance_summary.csv) records the revised image hierarchy.
- [`results/NICHE_PERMUTATION_RESULTS_2026-08-07.md`](results/NICHE_PERMUTATION_RESULTS_2026-08-07.md) records the 10,000-permutation environmental-sorting null.
- [`results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md`](results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md) records residual Moran and region-omission diagnostics.
- [`results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md`](results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md) records synonym-collapse sensitivity.
- [`06_references.md`](06_references.md) is the submission-wide bibliography.
- [`07_data_code_availability.md`](07_data_code_availability.md) contains the blinded-review repository and durable-release statement.
- [`../analysis/ch1/pipeline.json`](../analysis/ch1/pipeline.json) maps canonical executable stages.
- [`EXTERNAL_COMPLETION_GATES.md`](EXTERNAL_COMPLETION_GATES.md) records completed versus genuinely external validation gates.
- [`EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md`](EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md) fixes the Chapter 1 → EAzami boundary.

## Current submission blockers

Submission remains on hold for scientific validation that cannot be supplied by the completed computational rerun:

1. adjudicated human boxes for independent detector precision/recall;
2. independent reference measurements for orientation, colour and outline;
3. independent botanical calibration for expanded involucral/armature image phenotypes before botanical naming or functional promotion;
4. developmental-stage labels or an anthesis-restricted sensitivity;
5. paired flower-versus-background colour negative controls;
6. repeat-photo trait remeasurement and variance partitioning.

`visible_floret_fraction` remains unexecuted and is not used for a floral-display conclusion. The four unexecuted colour-composition endpoints are descriptive and do not reduce the 18-endpoint inferential analysis to a failed run.

Taxonomic robustness, residual spatial autocorrelation, broad-region omission, raw-calendar and hemisphere-aware cyclic collection timing, dominant-taxon omission, native-range restriction, environmental-niche permutation and candidate image-quality analyses are complete for the frozen operational-unit framework. Administrative metadata, licence selection, nomenclatural notes and final durable DOI release remain human/finalization tasks.

## Claim boundary

The final journal package must not:

- use the withdrawn negative lability relation or median-split quadrants as biological conclusions;
- use *climate tracking* or *environmental responsiveness* for cross-sectional spatial associations;
- equate visible image variance with genetic variance or plasticity;
- interpret hue sine/cosine separately as named flower colours;
- treat candidate image-geometry proxies as direct botanical spine/phyllary measurements before calibration;
- call surface/specular image signals hair density, gland density, mucilage or secretion without targeted validation;
- infer floral-display structure from the unexecuted `visible_floret_fraction` endpoint;
- claim resolved phylogenetic correction from the grafted historical layer;
- convert macro-scale associations into adaptation, pollinator causation, antagonist defence or evolutionary-rate claims.

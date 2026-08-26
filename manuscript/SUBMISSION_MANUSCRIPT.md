# Chapter 1 submission manuscript

This is the canonical submission-facing entry point for Chapter 1. The paper is framed as a **global quantitative phenotypic geography of thistle capitula**: public photographs are converted into repeated continuous trait observations, allowing spatial diversity within and among taxa and its association with environment to be examined without reducing morphology to species means or coarse categories.

The former lability-axis story is no longer a headline result. Mechanistic explanations are downstream hypotheses rather than the primary inferential target.

**Current status: SUBMISSION HOLD.** The computational raw-calendar and hemisphere-aware seasonal, dominant-taxon, native-range and involucre image-quality audits are complete; image remeasurement and independent validation gates listed below remain open.

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

## Current headline

The study moves from:

> species means / categorical trait states

through:

> public photographs → detected capitula → explicitly defined continuous measurements

to:

> repeated within-taxon phenotype distributions → global phenotypic geography → trait- and scale-specific environmental structure.

Across nine frozen primary continuous endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means, although biological and imaging variance are not yet fully separable. In the global within-taxon analysis, significant environmental structure involves orientation, visible colour and outline aspect ratio. Grouped SPDE-INLA models concentrate support in orientation and visible colour, while among-taxon permutation tests likewise identify the strongest environmental sorting in orientation and colour.

Robustness is graded rather than treated as one universal fail-closed filter. Orientation–BIO1 and aspect ratio–BIO4 currently have the strongest support across seasonal, dominant-taxon and native-range sensitivities. Chroma–BIO12 and aspect ratio–BIO12 remain same-signed global associations but are classified as range-sensitive because they lose native-only multiplicity support. High-resolution involucral contour signals are retained as exploratory because their strength depends on image resolution and sharpness; they are not promoted to direct botanical spine/phyllary or defensive-function claims.

The central biological contribution is therefore **the quantitative mapping of capitulum phenotype diversity and its environmental structure**, not the demonstration of a single adaptive mechanism. The central methodological contribution is a scalable route from heterogeneous public photographs to repeated continuous trait distributions.

## Why *Cirsium*

The same capitulum expresses conspicuous diversity in orientation, visible colour, gross outline and involucral architecture, allowing multiple phenotype modules to be measured on one reproductive unit. *Cirsium* also includes recent regional radiations and has a complex evolutionary history involving hybridization, incomplete lineage sorting, polyploidy and phenotypic convergence. This makes it an informative system for asking how visible morphology is distributed across geography and whether different phenotype components align with different environments. The current image analysis does **not** itself demonstrate adaptive radiation.

## Current main-result structure

1. **Continuous phenotypes can be recovered from global public imagery.**
2. **Taxon means conceal substantial repeated visible variation.**
3. **The below-taxon pattern survives one-head-per-photo and equal-replication controls.**
4. **The taxon-level phenotype remains multidimensional rather than collapsing onto one morphology axis.**
5. **Within-taxon environmental associations are trait specific, involving orientation, visible colour and selected outline dimensions.**
6. **Grouped spatial models concentrate the most consistent support in orientation and visible colour.**
7. **High-resolution involucral architecture contains provisional continuous environmental signals but is sensitive to image quality.**
8. **Among-taxon permutation tests concentrate environmental sorting mainly in orientation and visible colour.**
9. **Historical/PGLS results are Supporting Information placement sensitivities rather than resolved phylogenetic correction.**
10. **Seasonal, dominant-taxon, native-range and image-quality analyses grade confidence in individual associations rather than redefining the global primary question.**

## Evidence tiers

The submission uses the pre-result interpretation policy in [`../analysis/ch1/evidence_tiering_policy.json`](../analysis/ch1/evidence_tiering_policy.json):

- **A — robust:** supported globally and consistent across the relevant sensitivity analyses;
- **B — supported but sensitive:** supported in the global primary/spatial analysis but weakened in a restricted-domain sensitivity while retaining direction;
- **C — exploratory continuous signal:** coherent image-derived candidate association requiring stronger image-quality and/or botanical validation;
- **D — unresolved/artifact-prone:** important sign reversal or measurement-integrity problem;
- **validation-only:** surface signals that cannot yet be given a botanical ecological interpretation.

Sensitivity analyses therefore rank evidence strength. A result is not automatically erased solely because a restricted-domain BH value crosses 0.05.

## Figure priorities

The main figures should make the phenotype landscape and scale structure visible before presenting sensitivity gates:

1. **image → quantitative phenotype → ecological scale:** actual photographs, detection and continuous measurement route;
2. **geographic sampling and analytical domain:** global sampling, filtering streams and occupied environmental space;
3. **nested visible variation:** taxon → photograph → head decomposition across continuous endpoints;
4. **multivariate capitulum morphospace:** taxon-level positions and trait directions rather than categorical morphology classes;
5. **environmental structure across traits:** effect sizes and environmental axes, annotated by evidence tier rather than reduced to pass/fail cells.

Historical-placement sensitivity and detailed sampling/image-quality diagnostics remain Supporting Information where needed.

## Scientific control files

- [`COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](COHORT_FLOW_AND_ANALYSIS_LEDGER.md) fixes cohort names, counts and permitted analyses.
- [`final_claims.json`](final_claims.json) contains the machine-readable numerical claim registry.
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

Submission remains on hold. The following scientific items are unresolved:

1. adjudicated human boxes for independent detector precision/recall;
2. independent reference measurements for orientation, colour and outline;
3. developmental-stage labels or an anthesis-restricted sensitivity;
4. paired flower-versus-background colour negative controls;
5. repeat-photo trait remeasurement and variance partitioning.

The expanded continuous-trait run additionally evaluates explicitly defined involucral architecture, armature and surface-image metrics. New endpoints are not promoted to botanical or functional claims without the validation appropriate to their tier.

Taxonomic robustness, residual spatial autocorrelation, broad-region omission, raw-calendar and hemisphere-aware cyclic collection timing, dominant-taxon omission, native-range restriction and environmental-niche permutation analyses are complete for the frozen operational-unit analysis. Administrative metadata, licence selection, nomenclatural notes and final durable DOI release remain human/finalization tasks.

## Claim boundary

The final journal package must not:

- use the withdrawn negative lability relation or median-split quadrants as biological conclusions;
- use *climate tracking* or *environmental responsiveness* for cross-sectional spatial associations;
- equate visible image variance with genetic variance or plasticity;
- interpret hue sine/cosine separately as named flower colours;
- treat auxiliary image-geometry proxies as direct botanical spine/phyllary measurements before calibration;
- call surface/specular image signals hair density, gland density, mucilage or secretion without targeted validation;
- claim resolved phylogenetic correction from the grafted historical layer;
- convert macro-scale associations into adaptation, pollinator causation, antagonist defence or evolutionary-rate claims.

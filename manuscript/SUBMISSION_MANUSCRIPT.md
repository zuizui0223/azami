# Chapter 1 submission manuscript

This is the canonical submission-facing entry point for Chapter 1. The current story is the **global recovery of continuous within-taxon capitulum trait distributions from public photographs**, followed by scale-explicit environmental analysis. The former lability-axis story is no longer a headline result.

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

`06_references.md` is the submission-wide reference list. The older Introduction-only literature files remain drafting provenance and are not the reference source for the final manuscript.

## Current headline

Public biodiversity photographs can be converted into repeated **numerical phenotype observations** rather than one species-level category. Across nine primary continuous endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means. Some of that repeated within-taxon variation is structured along environmental gradients, but the pattern is trait- and scale-specific: orientation and visible colour show the clearest environmental structure, most gross-outline traits are weaker after spatial modelling, and fine involucral contour architecture is associated with temperature seasonality in the high-resolution layer.

The central biological conclusion is that capitulum phenotype is **modular rather than one climate-associated thistle syndrome**: orientation, visible colour, gross outline and fine involucral architecture show different environmental signatures, and within-taxon versus among-taxon analyses emphasize partly different gradients. The literature-backed mechanistic interpretation therefore generates separate hypotheses about reproductive microclimate/exposure, colour-mediated abiotic and biotic interactions, and antagonist/protection functions of involucral architecture rather than assigning one adaptive mechanism to all traits.

The central methodological contribution is therefore:

> global public photographs → detected capitula → continuous head-level measurements → within-taxon trait distributions → multiscale environmental analysis

The manuscript does not reduce the primary inference to upright/nodding, pale/pigmented or globose/elongate categories.

## Why *Cirsium*

The same capitulum expresses conspicuous diversity in orientation, visible colour, gross outline and involucral architecture, allowing multiple phenotype modules to be measured on one reproductive unit. *Cirsium* also includes recent regional radiations and has a complex evolutionary history involving hybridization, incomplete lineage sorting, polyploidy and phenotypic convergence. This makes it an informative system for asking whether trait diversity is environmentally structured within taxa, among taxa, or both. The current image analysis does **not** itself demonstrate adaptive radiation.

## Current main-result structure

1. **Repeated continuous phenotypes are recoverable globally.**
2. **Species/taxon means conceal substantial below-taxon visible variation.**
3. **That conclusion survives one-head-per-photo and equal-replication controls.**
4. **Within-taxon environmental associations are trait specific.**
5. **Grouped SPDE-INLA narrows the strongest robust pattern toward orientation and visible colour.**
6. **High-resolution involucre architecture shows a small coherent temperature-seasonality association.**
7. **Among-taxon permutation tests concentrate environmental sorting mainly in orientation and colour.**
8. **Historical/PGLS results remain conditional on incomplete direct-backbone coverage.**
9. **Spatial and taxonomic robustness diagnostics do not overturn the headline results.**

## Figure priorities

The six main figures are now organized to minimize table-like cognitive load and make the study design and scale-specific results visually explicit:

1. **image → phenotype → ecological questions:** retain actual-photo YOLO/measurement evidence and connect repeated numerical phenotypes to the three executed inferential scales;
2. **geographic sampling and analytical domain:** show detector-positive sampling density, the primary spatial cohort, the two analysis streams and occupied BIO1 × BIO12 space;
3. **nested visible variance:** retain the stacked taxon → photograph → head decomposition across all nine primary endpoints;
4. **actual taxon-level multivariate structure:** show PC1–PC2 taxon scores with all nine trait directions and expose PC3 rather than presenting a loading table alone;
5. **environmental effect sizes:** replace a symbol-only support matrix with primary, grouped-SPDE and high-resolution coefficient intervals while keeping the three inferential families visually separate;
6. **ecological scale specificity + historical sensitivity:** juxtapose supported within-taxon and among-taxon climate axes before showing the six randomized-tree PGLS effect sizes.

Supplement Figures S6–S8 then carry the sampling-filter geography, geographic trait assessability and observed among-taxon environmental-separation geometry that would otherwise make the main figures too dense. The legacy precision-confounding audit and withdrawn lability quadrants remain provenance/supplementary QA, not the main narrative.

## Scientific control files

- [`COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](COHORT_FLOW_AND_ANALYSIS_LEDGER.md) fixes cohort names, counts and permitted analyses.
- [`final_claims.json`](final_claims.json) contains the frozen machine-readable claims.
- [`FIGURE_TABLE_MAP.md`](FIGURE_TABLE_MAP.md) fixes the main/Supplement figure–table crosswalk.
- [`MAIN_FIGURE_CAPTIONS.md`](MAIN_FIGURE_CAPTIONS.md) contains the authoritative captions for the six main figures.
- [`FIGURE2_GEOGRAPHY_MANIFEST.md`](FIGURE2_GEOGRAPHY_MANIFEST.md) records the frozen sources for main Figure 2 and Figures S6–S7.
- [`INTERPRETIVE_FIGURES_MANIFEST.md`](INTERPRETIVE_FIGURES_MANIFEST.md) records the frozen sources and claim boundaries for main Figures 1/4/5/6 and Figure S8.
- [`results/nested_visible_variance_summary.csv`](results/nested_visible_variance_summary.csv) records the revised image hierarchy.
- [`results/NICHE_PERMUTATION_RESULTS_2026-08-07.md`](results/NICHE_PERMUTATION_RESULTS_2026-08-07.md) records the 10,000-permutation niche null.
- [`results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md`](results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md) records residual Moran and region-omission diagnostics.
- [`results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md`](results/WCVP_TAXONOMIC_SENSITIVITY_2026-08-07.md) records synonym-collapse sensitivity.
- [`06_references.md`](06_references.md) is the submission-wide bibliography for citations in the manuscript.
- [`../analysis/ch1/pipeline.json`](../analysis/ch1/pipeline.json) maps canonical executable stages.
- [`EXTERNAL_COMPLETION_GATES.md`](EXTERNAL_COMPLETION_GATES.md) records completed versus genuinely external validation gates.
- [`EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md`](EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md) fixes the Chapter 1 → EAzami boundary.

## Current submission blockers

Only two scientific items remain genuinely external to repository computation:

1. adjudicated human boxes for independent detector precision/recall;
2. independent reference measurements for orientation, colour and outline.

Taxonomic robustness, residual spatial autocorrelation, broad-region omission and environmental-niche permutation gates are complete for the frozen operational-unit analysis. Administrative metadata, nomenclatural notes and final durable DOI release remain human/finalization tasks.

## Claim boundary

The final journal package must not:

- use the withdrawn negative lability relation or median-split quadrants as biological conclusions;
- use *climate tracking* or *environmental responsiveness* for cross-sectional spatial associations;
- equate visible image variance with genetic variance or plasticity;
- interpret hue sine/cosine separately as named flower colours;
- treat auxiliary image-geometry proxies as direct botanical spine/phyllary measurements;
- claim resolved phylogenetic correction from the grafted historical layer;
- convert macro-scale associations into adaptation, pollinator causation or evolutionary-rate claims.

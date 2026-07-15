# Chapter 1 final manuscript strategy

## One-sentence claim

Species means conceal a scale-dependent and modular capitulum phenotype: most visible variation occurs within species, but total variation, within-species climate tracking and among-species environmental sorting are not equivalent.

## Problem

Comparative trait ecology usually represents each species by one mean. That practice assumes that species differences contain the most relevant trait information and often treats “lability” as a single conserved-to-variable property. The executed Chapter 1 results show that three structures must instead be separated:

1. the amount of visible phenotype distributed within a species;
2. the extent to which phenotype tracks environmental gradients within species;
3. the environmental sorting of species with contrasting mean phenotypes.

Public photographs provide enough replication to separate these structures, but only when measurement assessability, viewpoint, circular hue, spatial structure and uncertain historical placement are handled explicitly.

## Why Cirsium is a powerful test rather than merely a convenient case study

*Cirsium* is species rich, morphologically diverse, taxonomically difficult and affected by hybridization, polyploidy and uncertain generic or sectional limits. These features mean that visible capitulum traits need not map neatly onto one bifurcating species tree. The genus therefore provides a strong test of whether species means, single lability rankings and single-tree explanations are adequate.

The paper should not claim a demonstrated adaptive radiation. Instead, it should state that regional endemism, hybridization, polyploid complexes and unresolved historical relationships make *Cirsium* a plausible rapidly and reticulately diversifying system in which repeated or weakly conserved trait combinations are expected.

## Evidence base

- 6,626 detected capitula from 3,725 observations and 216 accepted image-analysis taxa;
- nine primary continuous endpoints spanning orientation, colour and outline;
- 46,276 spatially thinned observations from 259 taxa in the strict within-species layer;
- 102 taxa in the complete two-axis lability cohort;
- 36 completed grouped SPDE-INLA models;
- 50 alternative grafting trees, 522 completed historical PGLS fits and 636 completed phylogenetic-signal fits.

## Main results to lead with

1. Across the nine primary endpoints, approximately 82–99% of visible variance is within species. This quantity may contain viewpoint and measurement components and must be described as visible variation, not pure biological variance.
2. Within the 102 complete taxa, species within-variation and environmental responsiveness are negatively associated (Spearman rho = -0.333). All four high/low combinations are well represented.
3. Colour has the highest median climate responsiveness but the lowest median total variation. Shape has the highest visible variation without the highest responsiveness.
4. Grouped spatial effects are concentrated in colour and orientation; no shape effect passes global BH support across the fitted predictor groups. Effect sizes are small.
5. Species quartiles show the largest environmental niche contrasts for orientation and hue/chroma, whereas aspect ratio, circularity and solidity show high niche overlap.
6. The strict <=10 km primary within-species cohort has zero BH-supported effects. An expanded pooled table has four small supported effects. These are different inferential summaries and must not be conflated.
7. Minimum-sample sensitivity preserves the two-axis rankings.
8. Apparent phylogenetic signal is not robust to direct-backbone restriction; only 54 of 216 taxa are direct dated-backbone tips.

## Do not lose the concrete trait–environment results

The conceptual framing must retain the biological detail.

- **Orientation:** the strongest environmental niche centroid contrast among the measured traits; grouped spatial models favour climate-only structure. Discuss precipitation exposure, radiation, cooling, display and stem mechanics only as candidate mechanisms.
- **Colour:** the highest median environmental responsiveness; hue/chroma show strong niche contrasts; soil improves preferred models for lightness and chroma. Discuss pigment physiology and pollinator perception as alternative mechanisms, not demonstrated causes.
- **Shape:** the highest visible dispersion but weak climate responsiveness, high niche overlap for most outline traits and no globally BH-supported grouped spatial effects. Use this as the clearest example that large variation is not equivalent to climatic tracking.

## Novelty

The paper is not merely a global atlas and not merely a new lability index. Its contribution is a distributional, scale-explicit and module-aware framework for comparative image phenomics:

- it demonstrates how much information a species mean discards;
- it distinguishes total variation from environmental tracking;
- it separates within-species association from among-species environmental sorting;
- it shows that colour, orientation and outline occupy different ecological structures;
- it models assessability instead of converting failed images into biological absence;
- it propagates historical-placement uncertainty rather than claiming a false precise species tree;
- it shows why hybridizing and polyploid lineages require distributional and uncertainty-aware trait analysis.

## Why this is stronger than a conventional trait–environment big-data study

A conventional analysis would join one species-level trait value to one climate value and fit a regression. This study instead retains individual observations, repeated measurements, measurement failure, spatial provenance and species-specific environmental slopes. Public photographs do not replace standardized trait databases; they add geographic and within-species replication that conventional databases rarely contain for visible floral traits.

The methodological novelty lies in the full chain:

1. provenance-bearing public observations;
2. capitulum detection and context-preserving crops;
3. deterministic continuous colour, orientation and outline measurements;
4. trait-specific QC and horizontal-mirror technical checks;
5. observation-level retention before aggregation;
6. spatial thinning and coordinate-quality cohorts;
7. separate among-species, within-species and lability analyses;
8. versioned outputs, software environments and checksums.

## Narrative architecture

### Introduction

1. Species means dominate comparative trait databases, although within-species distributions affect ecology, forecasting and evolutionary inference.
2. “Lability” conflates phenotypic dispersion, environmental tracking and among-species differentiation.
3. Public-image phenomics offers global replication but creates non-random assessability, viewpoint and colour limitations.
4. *Cirsium* is an especially informative test because its diversity is shaped by regional differentiation, hybridization, polyploidy and uncertain historical relationships.
5. Thistle capitula combine orientation, colour and outline modules that need not respond at the same scale.
6. We ask where visible variance lies, whether variation predicts climate tracking, whether modules differ, and whether conclusions survive spatial and historical sensitivity tests.

### Results

1. Sampling, measurement retention and assessability.
2. Partitioning of visible variance within and among species and the multivariate PCA architecture.
3. Species-level within-variation versus responsiveness and four quadrants.
4. Module differences in responsiveness and trait-extreme environmental niches.
5. Strict versus expanded within-species inference and grouped spatial-model sensitivity.
6. Historical-placement and molecular-data limitations.

### Discussion

1. **Species means conceal most visible heterogeneity.** The large within-species fractions establish why replicated distributions are necessary, while image limitations prevent interpreting all dispersion as biological variation.
2. **Variation is not responsiveness.** The negative correlation and occupied quadrants reject a one-dimensional conserved–labile continuum.
3. **Trait architecture is modular and scale dependent.** Colour and orientation show stronger environmental structure than most outline traits, while shape remains highly variable.
4. **Among-species sorting is not contemporary within-species response.** Species turnover and longer-term differentiation are plausible explanations, but no single mechanism is proven.
5. **Reticulate and polyploid history complicates trait conservatism.** Hybridization and genome duplication provide plausible background for weak trait–tree correspondence, but are not tested causes in Chapter 1.
6. **Public-image phenomics is useful when failures are data.** Assessability and provenance are part of the inference rather than preprocessing debris.
7. **Broader implication.** Trait databases should retain distributions, environmental slopes and measurement quality, not only means.

## Figure strategy

1. **Figure 1 — From photographs to distributions.** Workflow, repeated observations, trait definitions, retained and failed measurements.
2. **Figure 2 — Species means are incomplete.** Within- versus between-species visible variance fractions and the species-level PCA.
3. **Figure 3 — Main conceptual result.** Species within-variation versus environmental responsiveness for 102 taxa, with four quadrants.
4. **Figure 4 — Module and scale contrast.** Colour, orientation and shape responsiveness together with environmental niche contrasts of trait extremes and the specific environmental-model results.
5. **Figure 5 — Credibility and boundaries.** Minimum-n, coordinate/spatial cohort and direct-tip/tree-placement sensitivity.

Full coefficient heatmaps, individual species slopes, model-group tables, molecular database coverage and all historical fits belong in the supplement.

## Recommended title

> **Species means conceal scale-dependent and modular variation in thistle capitulum traits**

Alternative:

> **Public photographs reveal decoupled variation and climate tracking in a reticulate thistle lineage**

## Abstract logic

1. Comparative trait ecology commonly reduces species to means, although within-species variation and environmental tracking need not coincide.
2. Hybridization, polyploidy and uncertain historical relationships make this simplification especially problematic in *Cirsium*.
3. We extracted continuous orientation, colour and outline from public photographs while modelling assessability.
4. Most visible variance occurred within species, but species variation and responsiveness were negatively associated and occupied all four combinations.
5. Colour and orientation showed stronger climate tracking and among-species niche differentiation than most outline traits; strict within-species effects were weak and cohort dependent.
6. Historical conclusions were sensitive to grafted tree placement.
7. Distributional, scale-explicit and module-aware phenomics therefore reveals ecological structure hidden by species means.

## Claim boundary

Use “visible within-species variation”, “environmental responsiveness”, “climate tracking”, “among-species environmental sorting”, “module-dependent”, “reticulate or polyploid evolutionary context” and “historical sensitivity”.

Do not claim demonstrated plasticity, local adaptation, adaptive radiation, pollinator causation, rain protection, evolutionary rate, hybrid causation, polyploid causation or a resolved *Cirsium* species tree.

## Submission strategy

Position the manuscript as conceptual macroecology enabled by public-image phenomics, with *Cirsium* providing a biologically compelling system in which species means and single-tree narratives are particularly questionable. The general message concerns how comparative trait data should represent species, biological scale and historical uncertainty.

No new headline analysis should be added. Remaining work is limited to credibility gates: detector evaluation, measurement demonstration, residual spatial/region diagnostics, taxonomic freeze and durable archiving.

The detailed evidence inventory and alternative story comparison are in `manuscript/RESULT_LED_STORY_REVIEW.md`. The genus-specific biological synthesis and candidate Discussion text are in `manuscript/CIRSIUM_EVOLUTIONARY_ECOLOGY_AND_DISCUSSION_STRATEGY.md`.

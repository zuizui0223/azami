# Chapter 1 figure and table map

This file maps manuscript items to frozen outputs. It prevents a draft from silently mixing analyses, cohorts or artifact versions.

## Main figures

### Figure 1 — From public photographs to trait distributions

Status: presentation assembly required.

Required components:

- source-image identifier and licence metadata;
- detected capitulum bounding box;
- tight crop and context crop;
- orientation, colour and outline measurement overlays;
- representative assessable and unassessable cases;
- mirrored technical replicate example.

Do not embed redistributed source images unless their licences permit it.

### Figure 2 — Species means conceal visible heterogeneity

Sources:

- frozen within-/among-species variance partition table from the integrated continuous analysis;
- species-level PCA summary and loading table.

Headline values:

- approximately 0.817–0.988 of visible variance occurs within species across primary endpoints;
- PC1 = 32.9%, PC1–2 = 56.1%, PC1–3 = 69.3%.

Required caveat: visible variance includes biological heterogeneity, viewpoint and measurement components.

### Figure 3 — Variation versus environmental responsiveness

Frozen outputs:

- `species_lability_axes.csv`;
- `figure_species_lability_quadrants.pdf`;
- `figure_species_lability_quadrants.png`.

Headline values:

- 102 primary-complete taxa;
- Spearman rho = -0.333;
- quadrants: 22 high/high, 29 low/high, 29 high/low, 22 low/low.

### Figure 4 — Trait-module and scale contrast

Frozen outputs:

- `module_lability_summary.csv`;
- `species_module_lability.csv`;
- trait-extreme environmental-niche summary;
- grouped SPDE-INLA model status, fixed-effect and WAIC tables.

Required visual elements:

- module medians for within-variation and responsiveness;
- standardized effect magnitude or support by module;
- environmental centroid distance and overlap for trait extremes.

Do not combine incomparable coefficient scales without explicit standardization.

### Figure 5 — Credibility and inference boundaries

Frozen outputs:

- `sensitivity_minimum_sample_size.csv`;
- strict versus expanded within-species coefficient summaries;
- spatial/cohort sensitivity outputs;
- phylogenetic evidence classification;
- direct-tip and alternative-tree summaries.

The figure should distinguish sampling sensitivity, spatial sensitivity and historical-placement sensitivity rather than merging them into one significance score.

## Main tables

### Table 1 — Sampling and assessability

Report observations, detected heads, taxa, retained measurements and failure categories by trait module.

### Table 2 — Scale-specific environmental evidence

Rows should be endpoint–predictor combinations. Separate columns must identify:

- among-species niche contrast;
- strict within-species estimate and BH result;
- expanded pooled estimate and BH result;
- grouped SPDE model support;
- interpretation scope.

Never merge strict zero-FDR and expanded four-FDR summaries into one count.

## Supplementary figures

- Figure S1: `figure_s1_effect_size_heatmap.pdf`;
- detector precision/recall and failure categories;
- measurement low/middle/high examples and mirrored repeats;
- full module-level species distributions;
- environmental PCA niche ellipses by trait extreme;
- alternative-tree coefficient and lambda distributions;
- taxonomic and tree-placement audit.

## Supplementary tables

- Table S1: `table_s1_all_global_coefficients.csv`;
- Table S2: `table_s2_species_specific_coefficients.csv`;
- species lability axes and quadrant membership;
- trait-specific within-variation;
- module lability summary;
- minimum-sample sensitivity;
- grouped SPDE model status, WAIC and coefficients;
- taxonomic accepted-name decision table;
- phylogenetic evidence summary;
- NCBI molecular database coverage;
- file-level provenance and SHA-256 manifest.

## Artifact rule

Before a value enters the manuscript, record the workflow run ID, artifact digest, Git commit and source filename in the drafting notes. GitHub Actions artifacts are temporary and must be copied into the durable submission bundle without renaming the internal source tables.

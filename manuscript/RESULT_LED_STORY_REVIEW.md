# Chapter 1 result-led story review

## Purpose

This document reviews the executed Chapter 1 results before choosing the manuscript story. It separates observed results from interpretation, compares several publishable narratives, and recommends the framing that maximizes novelty, generality and credibility without adding a new analysis scope.

## Executed evidence inventory

### Data and measurement scale

- 6,626 detected capitula from 3,725 public observations and 216 accepted taxa support the continuous between-species trait atlas.
- Nine primary endpoints cover orientation, corolla colour and capitulum outline; three auxiliary endpoints describe high-resolution involucre/spine proxies.
- The strict spatial within-species layer contains 46,276 observations from 259 taxa.
- The final two-axis species-lability cohort contains 102 taxa after requiring at least 10 observations, six eligible traits and three environmental predictors per trait.
- Measurement assessability is treated as an analysed outcome. Unmeasurable images are not recoded as biological absence.

### Result family 1 — most visible variance occurs within species

Across the nine primary endpoints, estimated within-species fractions of visible variance range from approximately 0.817 to 0.988.

- lightness: 0.817;
- chroma: 0.842;
- hue sine: 0.842;
- hue cosine: 0.950;
- orientation: 0.971;
- circularity: 0.976;
- width-profile CV: 0.980;
- solidity: 0.969;
- aspect ratio: 0.988.

This is not a claim that all variance is biological. Public-image viewpoint and measurement error can contribute, especially for orientation and outline. The defensible result is that species means discard a large amount of repeatable visible heterogeneity, so distributional summaries and explicit assessability are necessary.

### Result family 2 — variation and climate tracking are not the same property

For the 102 complete taxa, species within-variation and species environmental responsiveness have Spearman rho = -0.333.

The four quadrants contain:

- 22 high-variation/high-response taxa;
- 29 low-variation/high-response taxa;
- 29 high-variation/low-response taxa;
- 22 low-variation/low-response taxa.

The result is not driven by the minimum sample threshold. Relative to the primary n >= 10 analysis, species rankings remain high at n >= 5 and n >= 20.

### Result family 3 — trait modules occupy different parts of the lability space

Median species-level indices are:

| module | within-variation | environmental responsiveness |
|---|---:|---:|
| colour | 0.553 | 0.175 |
| orientation | 0.630 | 0.129 |
| shape | 0.665 | 0.153 |

Colour therefore has the highest median climate tracking but the lowest median total variation. Shape has the highest visible variation but does not have the highest responsiveness. A single conserved-to-labile ranking cannot represent this contrast.

The grouped spatial models support the same qualitative modularity. Across predictor-set fits, BH-supported effects are concentrated in colour, are fewer for orientation, and are absent for the four shape endpoints. Effect magnitudes remain small.

### Result family 4 — climate structure is scale dependent

The strict <=10 km primary within-species cohort has zero BH-FDR-supported primary endpoint-climate effects. An expanded pooled coefficient table contains four BH-supported effects, but all standardized coefficients are small.

In the grouped SPDE-INLA analysis, 36 models all completed successfully. The principal pattern is not a large universal environmental effect:

- climate-only models are best by WAIC for orientation and all four shape endpoints;
- adding soil improves the preferred model for lightness and chroma;
- the strongest repeated effects are concentrated in colour and orientation;
- shape effects remain weak despite large visible within-species variance.

Thus the evidence supports weak, model- and cohort-dependent within-species climate structure rather than strong contemporary adaptation.

### Result family 5 — trait extremes differ in environmental niche mainly for orientation and colour

Environmental PCA contrasts between low and high species quartiles show the largest centroid distances for:

- orientation: 2.398, overlap 0.746;
- hue cosine: 2.271, overlap 0.641;
- hue sine: 1.807, overlap 0.817;
- width-profile CV: 1.663, overlap 0.849;
- chroma: 1.632, overlap 0.694.

Aspect ratio, circularity and solidity have much higher niche overlap (approximately 0.91–0.96). This descriptive result indicates stronger among-species environmental sorting for orientation and colour than for most outline traits.

### Result family 6 — the multivariate architecture is not one-dimensional

Species-level PCA explains:

- PC1: 32.9%;
- PC1 + PC2: 56.1%;
- PC1–PC3: 69.3%.

No single axis dominates the nine traits. The loading structure mixes colour and outline contrasts rather than reproducing a simple orientation-colour-shape ordering. This supports a modular, multidimensional interpretation rather than one latent floral-lability axis.

### Result family 7 — historical interpretation is limited by tree information

- 54 of 216 taxa are direct dated-backbone tips; 162 are grafted.
- All 522 historical PGLS fits and all 636 K/lambda signal fits completed.
- No non-circular endpoint retains FDR support in the direct-backbone subset.
- circularity and solidity are classified as grafting-sensitive;
- hue requires a joint circular test and cannot be interpreted component by component.

The correct conclusion is not that phylogeny is irrelevant. It is that current data do not support a robust, precise species-level historical explanation for the observed trait patterns.

### Result family 8 — existing molecular data cannot presently solve the history at manuscript scale

Among 216 taxa, NCBI coverage is:

- nucleotide: 138;
- ITS: 133;
- plastid: 129;
- complete plastome: 11;
- SRA: 120.

These records are heterogeneous and do not constitute a uniform nuclear species-tree dataset. In a hybridizing and polyploid genus, a plastid-only tree would not solve reticulation or allopolyploid parentage.

## Competing manuscript stories

### Story A — “Trait lability has two axes”

**Question:** Are species that are highly variable also those that track climate strongly?

**Core evidence:** rho = -0.333, four occupied quadrants, module contrasts, sample-size sensitivity.

**Strengths:** clear conceptual result; easy main figure; general beyond Cirsium; directly supported by the newest analysis.

**Weakness:** can sound like a new index paper unless anchored to the broader scale problem and observed variance structure.

**Verdict:** strong core, but should be embedded in a larger scale-explicit argument.

### Story B — “Species means erase most floral-trait structure”

**Question:** What is lost when comparative ecology represents a species by one mean?

**Core evidence:** 82–99% visible within-species variance, weak correspondence between variation and environmental responsiveness, large and balanced four-quadrant occupancy.

**Strengths:** high generality and importance for trait databases, macroecology and image phenomics; explains why replicated public photographs matter.

**Weakness:** within-species fractions include viewpoint and measurement components. The language must remain “visible variation” and be paired with assessability/QC evidence.

**Verdict:** potentially the most important framing when combined with Story A.

### Story C — “Climate structure changes with biological scale”

**Question:** Why do trait-climate patterns appear among species but weaken within species?

**Core evidence:** species-level niche contrasts, strict within-species zero-FDR result, four small expanded pooled effects, grouped spatial models and historical sensitivities.

**Strengths:** strong macroecological problem; links species turnover, niche differentiation and contemporary variation; makes scale explicit.

**Weakness:** several model families and cohorts must be explained carefully. A careless abstract could appear internally contradictory.

**Verdict:** essential second layer of the final story, not the sole headline.

### Story D — “Floral architecture is modular”

**Question:** Do orientation, colour and outline respond as one integrated floral phenotype?

**Core evidence:** module medians, concentration of spatial effects in colour/orientation, niche contrasts, multidimensional PCA.

**Strengths:** biologically intuitive and visually attractive; connects to future pollination, rain exposure and developmental studies.

**Weakness:** orientation is represented by one primary endpoint while colour and shape have four each; mechanism is not tested.

**Verdict:** strong biological interpretation layer, but not sufficiently general as the only headline.

### Story E — “Public-image phenomics is a new comparative-data source”

**Question:** Can uncontrolled biodiversity photographs support global, replicated continuous-trait ecology?

**Core evidence:** 6,626 heads, 216 taxa, continuous measurements, assessability models, mirror repeats, preserved failures and full provenance.

**Strengths:** broad methodological relevance; distinguishes the study from conventional trait databases.

**Weakness:** detector precision/recall and the final measurement demonstration remain manuscript gates. A methods-led title would invite stricter validation expectations.

**Verdict:** the enabling contribution, not the sole result claim.

### Story F — “Phylogenetic uncertainty can create apparent trait conservatism”

**Question:** How stable is phylogenetic signal when most genus tips are grafted?

**Core evidence:** direct-tip versus all-tip discrepancy and 50-tree sensitivity.

**Strengths:** honest and methodologically useful.

**Weakness:** the study lacks the nuclear history required for a positive evolutionary result; this would turn a limitation into the whole paper.

**Verdict:** supplementary robustness result only.

## Recommended synthesis

The strongest manuscript combines Stories B, A, C and D in that order, with Story E as the enabling method and Story F as a credibility check.

### Central problem

Comparative trait ecology commonly collapses each species to one average and treats lability as one quantity. This can confuse three different structures:

1. how much visible phenotype varies within species;
2. whether variation tracks environment within species;
3. whether species with contrasting mean phenotypes occupy different environments.

### Central result

Most observed capitulum-trait variation is within species, yet high variation does not imply strong climate tracking. Climate structure is stronger and more consistent in colour and orientation and at the among-species/niche-sorting scale than in most shape traits or strict within-species slopes.

### General implication

Species means, single lability rankings and single-scale trait-environment models are insufficient. Comparative trait datasets should retain within-species distributions, environmental slopes and measurement assessability, and should distinguish trait modules.

### Why Cirsium is more than a case study

Cirsium provides broad geographic replication and multiple visible modules in one reproductive structure, while its reticulate/polyploid history prevents a falsely neat single-tree explanation. The genus therefore exposes both the opportunity and the limits of global image-derived trait macroecology.

## Recommended headline and title

### Headline

> Public-image phenomics shows that floral-trait variation, environmental tracking and among-species niche differentiation are non-equivalent, module-dependent structures.

### Preferred title

> **Species means conceal scale-dependent and modular variation in thistle capitulum traits**

### Alternative result-led title

> **Public photographs reveal decoupled variation and climate tracking in thistle capitulum traits**

### Alternative macroecological title

> **Within-species variation and among-species climate sorting structure thistle capitulum traits at different scales**

## Recommended abstract logic

1. Trait macroecology commonly represents species by a mean, although within-species variation and environmental tracking need not coincide.
2. We used public photographs to quantify continuous orientation, colour and outline while explicitly modelling assessability.
3. Visible within-species variance was large across all endpoints, but species variation and responsiveness were negatively associated and occupied all four combinations.
4. Colour and orientation showed stronger climate tracking and among-species environmental niche contrast than most outline traits; strict within-species effects were weak and cohort dependent.
5. Historical-signal conclusions were sensitive to grafted tree placement.
6. Distributional, scale-explicit and module-aware trait data therefore reveal ecological structure hidden by species means.

## Figure strategy under the recommended story

1. **Why means are insufficient:** workflow, repeated observations and distributions for representative species, including assessability failures.
2. **Where variation lies:** endpoint-level within- versus between-species variance fractions plus species-level PCA.
3. **Main conceptual figure:** within-variation versus responsiveness for 102 taxa, with four quadrants.
4. **Scale and module contrast:** colour/orientation/shape responsiveness alongside environmental niche contrasts of species quartiles.
5. **Credibility and boundaries:** minimum-n, coordinate/spatial cohort and direct-tip/tree-placement sensitivity.

The old figure order that leads with individual coefficients should be retired. Full coefficient heatmaps belong in the supplement because effect sizes are small and the conceptual result is multiscale structure, not a list of isolated beta estimates.

## Interpretation boundaries

Use:

- visible within-species variation;
- environmental responsiveness or climate tracking;
- among-species environmental sorting;
- modular and scale-dependent structure;
- historical sensitivity.

Do not use the observational results as proof of:

- plasticity;
- local adaptation;
- pollinator-mediated selection;
- rain protection;
- evolutionary rate;
- a resolved Cirsium species tree.

## Remaining work that can change credibility

1. detector precision, recall and missed-head categories;
2. representative low/middle/high measurement and QC panel;
3. residual spatial autocorrelation and leave-one-region/block sensitivity;
4. accepted-name/synonym decision table;
5. durable archive, licence audit and checksum manifest.

These checks can weaken or qualify claims, but no new headline data layer should be added.

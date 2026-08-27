# Capitulum functional-space hypotheses for Chapter 1 — 2026-08-27

## Status and purpose

This document freezes a literature-motivated extension of the GEB v2 Chapter 1 analysis **before process-environment outcomes are inspected**. The extension does not replace the established 9-primary-trait evidence atlas. It asks a higher-level spatial question using the completed 18-endpoint inferential trait universe:

> **How is a multidimensional capitulum phenotype organized within and among taxa, and do predeclared environmental process blocks align with the same phenotypic organization at both scales?**

The chapter remains observational. “Function” below means a candidate functional context supported by independent literature, not a demonstrated function of an image-derived endpoint.

## Why return from separate traits to the capitulum

Asteraceae capitula are condensed inflorescences that function as flower-like reproductive units (pseudanthia). Their component structures can contribute to attraction/presentation, reproductive protection and the geometry of the reproductive unit. Floral integration theory also predicts that subsets of reproductive structures may be internally integrated while retaining partial independence from other subsets. Therefore the GEB v2 endpoint expansion should not end as a list of unrelated trait-by-climate coefficients. It enables a test of whether the observed capitulum occupies a structured, partly modular phenotypic space.

This does **not** imply that the registry modules are already validated functional modules. They were defined primarily from biological construct and measurement scope, and shared image segmentation can also generate covariance. The analysis therefore uses the terms **measurement-module organization**, **partial phenotypic modularity**, and **functionally annotated phenotype space**.

## Frozen hypotheses

### H1 — Measurement-module organization exists at both biological scales

For the 17 inferential units represented by the 18 measured endpoints (circular hue treated jointly), mean association strength among units in the same registered measurement module will exceed mean association strength between modules.

- within-taxon scale: association matrix from taxon-centred observations with equal total taxon weight;
- among-taxon scale: association matrix from taxon medians;
- test: module-label permutation preserving module sizes;
- interpretation: evidence of organized phenotype covariance, **not** proof of functional or genetic modularity.

### H2 — Within- and among-taxon capitulum organization is only partially shared

The upper triangles of the within- and among-taxon inferential-unit association matrices will be positively related, but identity of the matrices is not assumed.

- statistic: Spearman matrix similarity;
- null: jointly permute trait labels on the among-taxon matrix;
- uncertainty: taxon-cluster bootstrap;
- interpretation: shared but scale-dependent phenotypic organization.

This hypothesis directly extends the univariate PR #72 result, which found only partial overlap between within-taxon environment associations and among-taxon environmental sorting.

### H3 — Environmental structure should be tested as process blocks, not as a larger BIO-variable fishing screen

The four established CHELSA variables remain the frozen baseline. The extension adds a small set of process-proximal blocks selected from the literature before their outcomes are inspected:

1. `core_thermal`: BIO1 + BIO4;
2. `core_precipitation`: BIO12 + BIO15;
3. `radiative_atmospheric_drying`: mean shortwave radiation (`rsds`) + mean VPD;
4. `mechanical_exposure`: mean near-surface wind speed;
5. `growing_season_water_input`: growing-season precipitation;
6. `climatic_productivity`: potential NPP.

The whole 18D endpoint vector is the response. Each block receives one multivariate permutation test at the within-taxon scale and one at the among-taxon scale. This avoids treating 18 traits × many correlated environmental rasters as independent discoveries.

### H4 — Environmental effect geometry may rotate across biological scales

For each environmental block, standardized endpoint coefficients define a vector/matrix in the same 18D phenotype space. Within- and among-taxon coefficient geometry will be compared by cosine similarity with a taxon-cluster bootstrap.

A low or negative cosine is interpreted only as **different spatial association geometry across scales**. It is not evidence of opposing selection, plasticity or adaptation.

### H5 — Chapter 1 should export phenotype-space targets, not causal mechanisms, to EAzami

The new Azami → EAzami handoff adds observational targets that a downstream mechanism family may attempt to reproduce:

- within-taxon measurement-module contrast;
- among-taxon measurement-module contrast;
- within-vs-among association-matrix similarity;
- multivariate R² for each predeclared environment block at each scale;
- cross-scale coefficient cosine for each block;
- endpoint/module coefficient-energy allocation as a descriptive target.

EAzami may use these targets to discriminate common-lability and modular mechanism families, but successful reproduction never upgrades Chapter 1 associations to causal claims.

## Literature-motivated candidate functional contexts

### Capitulum as an integrated reproductive unit

Zhang & Elomaa (2024) review Asteraceae capitula as condensed flower-like inflorescences whose architecture is central to the family. Baczyński & Claßen-Bockhoff (2023) review pseudanthia more generally. These sources motivate analysis at the level of the whole capitulum rather than only independent endpoints.

### Partial modularity

Diggle (2014) reviews floral integration and the possibility that subsets of structures involved in different functions form partly independent modules. This motivates H1/H2 but does not establish that the present image modules correspond to those functional modules.

### Advertisement can expose plants to both mutualists and antagonists

In *Cirsium purpuratum*, flower/head production is associated with predispersal seed-predator damage (Ohashi & Yahara 2000). Independent *Cirsium* work also shows that floral signals can be shared by pollinators and floral herbivores. This motivates an advertisement–enemy trade-off in EAzami, but Chapter 1 does not yet contain a calibrated whole-plant display-size trait.

### Orientation can mediate abiotic exposure as well as pollinator presentation

Experiments in nodding Asteraceae and in sunflower show that orientation can alter reproductive exposure, floral thermal/timing conditions and pollination-related outcomes. These studies motivate radiation, atmospheric drying, growing-season water and wind as process contexts. They do not demonstrate that the global *Cirsium* orientation association has the same mechanism.

### Flower colour can have both biotic and abiotic environmental structure

Global comparative work shows that flower colour covaries with biotic context as well as precipitation and solar radiation. The Chapter 1 JPEG-derived colour measurements are not spectral reflectance, so the process extension tests spatial association only. A flower-versus-background image contrast remains a useful negative-control/visibility extension, while true pollinator visual-space inference belongs to calibrated spectral work.

### Water balance requires more than annual precipitation

Water-deficit syntheses and physiological work on floral water loss motivate VPD/water-balance context in addition to BIO12/BIO15. VPD is atmospheric drying power; it is not a direct measure of floral water status.

## Pre-result environmental design rule

Do **not** add all 19 BIOCLIM variables or select predictors after looking at trait coefficients. New environmental variables must enter through the frozen process blocks above. If a data source cannot be extracted reproducibly, that block is recorded as unavailable rather than replaced by a post hoc proxy.

The primary four-variable evidence atlas remains the established result family. The process-block analysis is a multivariate, hypothesis-driven extension.

## Current pre-process results from the already available four-BIO-variable data

These values were calculated after H1/H2 and the block design were formalized locally, but **before** any outcomes from the new process-extension rasters were available.

### Multilevel capitulum-space organization

Using observations complete for all 18 measured inferential endpoints:

- complete 18D observations: 1,874 from 124 taxa;
- main replication threshold (>=5 complete observations/taxon): 1,734 observations, 42 taxa;
- sensitivity threshold (>=2): 1,825 observations, 75 taxa.

Main >=5 scope:

- within-taxon module contrast = **0.16450**, module-label permutation `p = 0.00180`;
- among-taxon module contrast = **0.08848**, `p = 0.02230`;
- within-vs-among association-matrix Spearman = **0.36630**, trait-label permutation `p = 0.00010`;
- 95% taxon-bootstrap intervals: within contrast **0.13269–0.17740**, among contrast **0.02964–0.13156**, cross-scale Spearman **0.08909–0.40149**.

The >=2 sensitivity gives similar values (within contrast 0.15769; among contrast 0.08366; cross-scale Spearman 0.37727). Thus the registered measurement modules have detectable internal organization at both scales, but the full organization is only partly shared across scales.

### Baseline whole-space climate blocks

With only the frozen four BIO variables available, the main >=5 scope gives:

- thermal block within taxa: multivariate R² = **0.00851**, permutation `p = 0.0158`, two-block BH `q = 0.0316`;
- thermal block among taxa: R² = 0.04845, `q = 0.893`;
- precipitation block within taxa: R² = 0.00544, `q = 0.453`;
- precipitation block among taxa: R² = 0.03122, `q = 0.893`.

At the >=2 sensitivity threshold, the within-taxon thermal block no longer passes BH (`q = 0.1568`). Therefore this new whole-space thermal result is **threshold-sensitive exploratory multivariate structure**, not a replacement for or promotion beyond the established endpoint-level A/B evidence.

The cross-scale coefficient cosines are descriptive and currently uncertain; bootstrap intervals include zero. They must not be interpreted as evidence of opposing selection.

## Claim boundary

Chapter 1 may claim:

- multidimensional capitulum phenotype space;
- measurable organization of registered trait modules;
- partial rather than complete correspondence between within- and among-taxon trait organization;
- multivariate spatial association of the phenotype space with predeclared environmental blocks, when supported;
- scale dependence and uncertainty of those associations.

Chapter 1 may **not** claim from these analyses alone:

- validated functional modularity;
- genetic modularity or evolvability;
- adaptive peaks or many-to-one fitness solutions;
- plasticity or local adaptation;
- pollinator-mediated selection;
- defensive efficacy of involucral/armature image geometry;
- direct rain, UV, thermal or aerodynamic protection;
- a resolved causal mechanism.

## Key references used to freeze the hypotheses

- Zhang T, Elomaa P. 2024. Development and evolution of the Asteraceae capitulum. *New Phytologist* 242:33–48. DOI: 10.1111/nph.19590.
- Baczyński J, Claßen-Bockhoff R. 2023. Pseudanthia in angiosperms: a review. *Annals of Botany* 132:179–202. DOI: 10.1093/aob/mcad103.
- Diggle PK. 2014. Modularity and intra-floral integration in metameric organisms: plants are more than the sum of their parts. *Philosophical Transactions of the Royal Society B* 369:20130253. DOI: 10.1098/rstb.2013.0253.
- Ohashi K, Yahara T. 2000. Effects of flower production and predispersal seed predation on reproduction in *Cirsium purpuratum*. *Canadian Journal of Botany* 78:230–236. DOI: 10.1139/cjb-78-2-230.
- Dalrymple RL et al. 2020. Macroecological patterns in flower colour are shaped by both biotic and abiotic factors. *New Phytologist*. DOI: 10.1111/nph.16737.
- Kuppler J et al. 2020. Global gradients in intraspecific variation in vegetative and floral traits are partially associated with climate and species richness. *Global Ecology and Biogeography* 29:992–1007. DOI: 10.1111/geb.13077.
- CHELSA v2.1 climatologies and BIOCLIM+ documentation: https://www.chelsa-climate.org/ .

The final submission bibliography should replace abbreviated entries here with the canonical metadata in `manuscript/06_references.md`.

# Chapter 1 manuscript plan

Status: fixed-scope manuscript plan  
Decision date: 2026-07-11

## Scope decision

Chapter 1 will be completed from the analysis that already exists. It will not
be expanded with MaxEnt, LGM hindcasting, palaeoclimate stability, multi-regime
OU models, hair, mucilage or a new genomic phylogeny.

Those topics remain possible follow-up studies. They are not required to make
the present chapter publishable and would introduce new occurrence datasets,
new assumptions and a second scientific story.

## Working contribution

> Public biodiversity photographs can be converted into continuous capitulum
> measurements at global scale. Across *Cirsium*, colour, image-referenced head
> orientation and outline show between-species climatic structure, but the same
> associations are not recovered as robust within-species environmental effects.

This is a macroecological pattern paper with a methodological contribution. It
is not a demonstration of local adaptation, phenotypic plasticity or a resolved
adaptive radiation.

## Working title options

1. **Global climatic structure of image-derived capitulum traits in *Cirsium***
2. **Citizen-science photographs reveal global variation in thistle capitulum traits**
3. **From public photographs to macroecology: continuous capitulum traits across *Cirsium***
4. **Between-species, not within-species, climatic structure in thistle capitulum traits**

Option 1 is the safest default. Option 4 is the strongest result-led title if the
within-versus-between contrast remains central after the spatial sensitivity
analysis.

## Primary research questions

### Q1. Can uncontrolled public photographs yield reproducible continuous
capitulum measurements?

Evidence:

- 6,626 detected capitula;
- explicit colour, outline and orientation QC;
- horizontal-mirror technical replicates;
- observation- and species-level medians and uncertainty summaries;
- retained and failed measurements reported separately.

### Q2. Is measurement assessability environmentally or taxonomically biased?

Evidence:

- trait-specific QC-retention models;
- taxon, geography and CHELSA summaries;
- all 216 taxa retained in species tables even when a particular measurement is
  unavailable.

### Q3. Do species-level capitulum traits covary with present climate?

Evidence:

- four predeclared CHELSA predictors;
- nine primary continuous endpoints;
- OLS reference and historical-constraint sensitivities;
- deterministic and 50 random within-genus grafting trees.

### Q4. Are the global associations also visible within species?

Evidence:

- within-species demeaning of traits and predictors;
- species-clustered uncertainty;
- full coordinate-eligible and positional-accuracy <=10 km cohorts;
- zero FDR-significant primary effects in the strict coordinate cohort.

## Analysis hierarchy

### Confirmatory core

- colour lightness and chroma;
- circular hue represented jointly by sine and cosine;
- orientation angle relative to EXIF-oriented image vertical;
- outline aspect ratio, circularity, solidity and width-profile variation;
- four predeclared CHELSA variables;
- strict positional-accuracy cohort;
- between-species historical sensitivity across alternative trees.

### Measurement-bias analysis

QC retention is analysed as an outcome. These results explain which images are
measurable; they are not biological trait effects.

### Exploratory supplement

- high-resolution involucre projection and spine-like contour proxies;
- full-coordinate within-species sensitivity;
- deterministic-tree auxiliary PGLS;
- taxonomic and tree-placement audits.

The auxiliary involucre results must not become the title, abstract conclusion or
main causal claim.

## Manuscript structure

## Introduction

### Paragraph 1 — why floral architecture matters

Introduce floral display as the interface between pollination and abiotic
exposure. Orientation, colour and shape can influence visual display,
temperature, rain exposure and pollen protection, but their relative roles are
system dependent.

### Paragraph 2 — the comparative-data bottleneck

Species descriptions commonly reduce variable morphology to categories and are
uneven across regions. Public biodiversity photographs contain visible
phenotypes at much broader geographic coverage, but uncontrolled viewpoint,
lighting and detectability create non-random measurement error.

### Paragraph 3 — why *Cirsium*

*Cirsium* combines substantial variation in head orientation, corolla colour,
capitulum outline and involucral morphology with broad Northern Hemisphere
climatic distributions. Its reticulate and taxonomically difficult history also
makes it a useful case for separating observed macroecological pattern from an
overconfident single-tree evolutionary explanation.

### Paragraph 4 — study aims and predictions

State four aims:

1. derive continuous capitulum measurements from public images;
2. quantify differential assessability;
3. test present-climate associations at between- and within-species scales;
4. test whether between-species conclusions survive alternative historical-tree
   placements.

Do not predict that every trait must show a specific adaptive direction. The
primary prediction is scale-dependent structure: species-level associations may
be stronger than within-species environmental responses.

## Methods

Recommended order:

1. taxon and public-image sampling frame;
2. one-observation/one-primary-photo and spatial balancing rules;
3. capitulum detector and crop reconstruction;
4. continuous colour, orientation and outline measurements;
5. trait-specific QC and mirrored technical replicates;
6. observation and species aggregation;
7. CHELSA extraction and coordinate-quality tiers;
8. QC-retention models;
9. within-species models;
10. between-species OLS and PGLS sensitivities;
11. alternative-tree audit;
12. exploratory high-resolution involucre proxies;
13. multiple-testing and interpretation rules.

## Results

Recommended order:

1. sampling and trait-retention accounting;
2. visual distributions and representative measurements;
3. QC-retention bias;
4. strict within-species results;
5. between-species climate associations;
6. robustness across random grafting trees;
7. auxiliary involucre results;
8. explicit result boundary.

The result boundary should be stated before the discussion:

> Climatic structure was detectable among species, whereas robust primary-trait
> associations were not detected within species in the strict coordinate cohort.

## Discussion

### Section 1 — main result

Discuss global species turnover and long-term ecological or historical
differentiation as plausible explanations for between-species structure.

### Section 2 — why between- and within-species results differ

Possible, non-exclusive explanations:

- geographic replacement of species;
- historical niche differentiation;
- limited within-species climatic coverage after strict coordinate filtering;
- trait measurement error attenuating within-species slopes;
- different evolutionary and ecological time scales.

Do not select one mechanism as proven.

### Section 3 — methodological contribution

Emphasize continuous rather than forced categorical measurement, explicit
assessability modelling, preservation of failed measurements and separation of
within- from between-species inference.

### Section 4 — historical constraint without a false precise phylogeny

Explain that only 54 of 216 taxa were direct dated-backbone tips. Alternative
placements test sensitivity; they do not resolve the reticulate species history.

### Section 5 — limitations

Mandatory limitations:

- image vertical is a proxy for gravity;
- colour is not camera calibrated;
- outline is viewpoint dependent;
- detector false negatives condition the analysed sample;
- public-image sampling is geographically biased;
- species taxonomy and grafted tree positions are uncertain;
- associations are non-causal.

### Section 6 — next tests

Place pollinator observations, field manipulation, genomic phylogeny,
polyploid/hybrid evidence, mucilage and palaeoniche modelling here as future
work, not as missing components of the present paper.

## Proposed figures

### Figure 1 — workflow and measurement definitions

- public observation;
- detected capitulum and context crop;
- colour mask;
- directed PCA orientation axis;
- outline measurements;
- QC and aggregation flow.

### Figure 2 — global sampling and assessability

- species/observation map;
- trait-specific QC-retention fractions;
- geographic or climatic QC-retention effects.

### Figure 3 — trait atlas

- species-level distributions for colour, orientation and outline;
- representative low, middle and high images;
- no forced botanical categories.

### Figure 4 — within-species climate effects

- standardized coefficients and confidence intervals;
- strict <=10 km cohort as the main panel;
- full-coordinate cohort in supplement.

### Figure 5 — between-species historical sensitivity

- stable coefficient ranges across 50 random trees;
- Pagel-lambda distributions;
- direct-backbone versus grafted-tip coverage annotation.

### Optional Figure 6 — scale contrast

A compact comparison of within-species and between-species coefficient direction
and support for the same endpoint–predictor combinations.

The involucre analysis belongs in a supplementary figure unless it becomes
necessary to explain a main result.

## Proposed tables

### Table 1

Sampling, QC-retention and taxon coverage by trait group.

### Table 2

Primary between-species effects that pass the predeclared random-tree stability
criterion.

### Supplementary tables

- full endpoint registry;
- all QC-retention coefficients;
- all within-species coefficients;
- all OLS/BM/lambda-PGLS coefficients;
- tree and taxonomic audit;
- auxiliary involucre results;
- analysis configuration and provenance.

## Abstract logic

1. Public photographs offer global phenotypic coverage but uncontrolled image
   conditions can bias measurable samples.
2. We measured continuous capitulum colour, orientation and outline from 6,626
   detected heads representing 216 *Cirsium* taxa.
3. We quantified QC retention, compared within- and between-species climate
   associations and tested historical sensitivity across alternative grafting
   trees.
4. Several between-species associations were stable to tree placement, whereas
   the strict within-species analysis yielded no FDR-significant primary effect.
5. Global capitulum–climate structure therefore reflects species turnover or
   longer-term differentiation more clearly than contemporary within-species
   environmental response, while public images can support comparative trait
   ecology when assessability is modelled explicitly.

## Excluded from Chapter 1

The following are explicitly out of scope:

- MaxEnt or other species distribution models;
- LGM hindcasting and palaeoclimate stability;
- BM-versus-OU trait-model comparison as a headline result;
- a new nuclear or plastid phylogeny;
- claims about pollinator preference or protection from rain;
- literature-LLM hair or mucilage states;
- experimental proof of adaptation;
- future-climate projections.

## Remaining work, in order

1. residual spatial-autocorrelation and leave-one-region-out sensitivity;
2. detector precision/recall and failure-category audit;
3. representative measurement/QC figure;
4. taxonomic-name decision table;
5. freeze final coefficient and figure tables;
6. draft Methods and Results directly from frozen outputs;
7. write Introduction and Discussion using the novelty boundary;
8. create validated submission bundle and durable archive.

No new broad data layer should be added before these tasks are complete.

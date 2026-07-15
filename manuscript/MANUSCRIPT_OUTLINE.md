# Chapter 1 manuscript writing outline

Working title:

> **Species means conceal scale-dependent and modular variation in thistle capitulum traits**

Alternative title:

> **Public photographs reveal decoupled variation and climate tracking in thistle capitulum traits**

## Abstract

1. Comparative trait studies usually reduce each species to a mean.
2. Public photographs permit replicated, observation-level trait measurement but introduce non-random assessability and viewpoint effects.
3. We extracted continuous orientation, colour and outline measurements from public *Cirsium* photographs and separated visible within-species variation, within-species climate tracking and among-species environmental sorting.
4. Approximately 82–99% of visible variance occurred within species; variation and responsiveness were negatively associated across 102 complete taxa; colour and orientation retained stronger environmental structure than most outline traits.
5. Strict within-species effects were weak and cohort dependent, and historical conclusions were sensitive to grafted tree placement.
6. Trait databases should retain distributions, environmental slopes and measurement quality rather than species means alone.

## Introduction

### Paragraph 1 — the species-mean problem

Explain why comparative ecology commonly represents a species with one trait value and what is lost when within-species distributions are discarded.

### Paragraph 2 — three structures commonly conflated

Separate:

- visible within-species dispersion;
- within-species environmental tracking;
- among-species environmental sorting.

Introduce the prediction that these structures need not covary.

### Paragraph 3 — opportunity and risk of public-image phenomics

Describe global replication and multi-trait extraction, followed immediately by detectability, viewpoint, colour calibration and spatial sampling limitations. State that assessability is analysed rather than treated as absence.

### Paragraph 4 — why *Cirsium*

Present capitulum orientation, corolla colour and outline as modules within one reproductive structure. Introduce the genus as morphologically diverse and complicated by regional differentiation, hybridization, polyploidy and uncertain species-level history, without claiming demonstrated adaptive radiation.

### Paragraph 5 — questions

1. Where does visible trait variance lie, within or among species?
2. Does high within-species variation predict strong environmental responsiveness?
3. Do orientation, colour and shape show different scale-dependent environmental structures?
4. Are conclusions robust to sampling, spatial and historical-placement uncertainty?

## Methods

### Data acquisition and sampling frame

- public biodiversity observations and accepted taxon list;
- one primary photograph per observation where applicable;
- licence-safe identifiers and provenance;
- distinction between the 216-taxon image atlas and the 259-taxon strict spatial layer.

### Capitulum detection and observation units

- define a visible capitulum;
- describe detector, tight crop and context crop;
- explain why detection, anthesis and trait assessability are separate tasks;
- report detector audit once completed.

### Continuous trait extraction

- orientation relative to EXIF-oriented image vertical;
- Lab lightness and chroma;
- circular hue as sine/cosine and presentation angle;
- outline aspect ratio, circularity, solidity and width-profile variation.

### Trait-specific assessability and technical checks

- unassessable is not biological absence;
- masking and failure categories;
- mirrored technical replicates;
- retained versus failed measurements.

### Environment and spatial data

- four predeclared CHELSA variables;
- coordinate-quality tiers and spatial thinning;
- grouped SPDE-INLA predictor sets where reported.

### Statistical hierarchy

1. partition visible variance within and among species;
2. estimate strict and expanded pooled within-species effects;
3. estimate species-specific standardized slopes;
4. calculate species within-variation and environmental-responsiveness indices;
5. compare trait modules and environmental niches of trait extremes;
6. test minimum-sample, spatial/cohort and historical-placement sensitivities.

### Multiple testing and interpretation

Name the cohort for every FDR count. Treat hue jointly. Define associations as observational and non-causal.

## Results

### 1. Sampling and assessability

Report 6,626 detected capitula, 3,725 observations and 216 accepted image-analysis taxa. Report trait-specific retention and failure patterns without interpreting missing measurements as trait absence.

### 2. Species means conceal visible heterogeneity

Lead with the approximately 82–99% within-species visible-variance fractions. Pair this immediately with the viewpoint and measurement limitation.

### 3. Variation is not responsiveness

Report the 102 complete taxa, Spearman rho = -0.333 and quadrant counts 22/29/29/22. Report minimum-n sensitivity.

### 4. Modules and concrete environmental structure

- colour: lowest median total variation, highest median responsiveness;
- shape: highest visible variation without highest responsiveness;
- orientation and hue/chroma: strongest among-species environmental-niche contrasts;
- grouped spatial effects concentrated in colour and orientation; shape effects weak.

### 5. Scale-dependent inference

Report zero BH-supported primary effects for the strict <=10 km cohort and four small effects for the expanded pooled table as distinct analyses. Summarize WAIC/model-group findings without turning model selection into causal evidence.

### 6. Historical sensitivity

Report 54/216 direct backbone tips, complete alternative-tree fits and the absence of direct-tip-supported non-circular signal. State that current molecular records do not form a uniform nuclear species-tree dataset.

## Discussion

### 1. Species means are inadequate summaries

Explain why observation distributions matter for comparative trait ecology while distinguishing visible heterogeneity from pure biological variance.

### 2. Total variation and climate tracking are distinct

Use the negative correlation, occupied quadrants and shape module as the clearest evidence against a single conserved–labile axis.

### 3. Orientation and colour retain biological specificity

Discuss plausible links with rain exposure, radiation, thermal environment, pigment physiology, nutrient conditions and pollinator display. Keep these as hypotheses, not demonstrated mechanisms.

### 4. Scale-dependent climate structure

Explain why among-species sorting may be clearer than strict within-species response: species turnover, longer-term differentiation, limited environmental coverage, measurement attenuation and different time scales.

### 5. Reticulation, polyploidy and historical uncertainty

Use the *Cirsium* literature to explain why one grafted bifurcating tree cannot be treated as a resolved species history. Do not infer that hybridization or polyploidy caused the observed trait patterns.

### 6. Public-image phenomics complements conventional trait databases

Emphasize repeated observations, continuous multi-trait extraction, explicit assessability, failure retention, geographic coverage and reproducible provenance.

### 7. Limitations and next tests

Detector audit, standardized photography, common gardens, pollinator observations, pigment chemistry, nuclear target capture, flow cytometry and explicit hybrid-network analyses.

## Writing rule

Every numerical claim must match `manuscript/final_claims.json` or a frozen output listed in `manuscript/FIGURE_TABLE_MAP.md`. Mechanistic statements must be labelled as hypotheses unless directly tested.

# azami — global image-derived capitulum traits in *Cirsium*

This repository contains a multi-chapter project on the ecology and evolution of
*Cirsium* floral architecture. Chapter 1 is the submission-focused component: a
global macroecological analysis of continuous capitulum traits measured from
public biodiversity photographs.

## Chapter 1 in one sentence

We quantify continuous corolla colour, capitulum outline and image-referenced
head orientation, model which photographs are measurable, compare climatic
structure among and within species, and test whether conclusions survive
uncertain historical-tree placement.

The analysis is observational and non-causal. Climate association is not treated
as proof of local adaptation, plasticity or selection.

## Current executed dataset

The frozen image-measurement layer contains:

- 6,626 detected capitula;
- 3,725 primary public observations;
- 216 accepted image-analysis taxa;
- nine primary continuous measurement components;
- exploratory high-resolution involucre/spine contour proxies.

The strict within-species climate cohort contains 46,276 spatially thinned
observations representing 259 taxa. It is kept separate from the smaller
head-level measurement accounting because the two tables answer different
questions and have different sampling units.

## Main scientific logic

Chapter 1 separates four levels that are often conflated:

1. **measurement assessability** — whether an uncontrolled photograph supports a
   valid trait estimate;
2. **among-species climatic structure** — whether species-level trait summaries
   covary with present climate;
3. **within-species climate association** — whether demeaned traits track
   environmental gradients within species;
4. **species lability decomposition** — whether a species is broadly variable,
   environmentally responsive, both or neither.

The final lability analysis uses two independent axes:

- a species within-variation index based on trait-standardized robust dispersion;
- a species environmental-responsiveness index based on species-specific
  standardized slopes.

Circular hue remains a joint sine/cosine endpoint. A scalar hue angle is used for
presentation, not as an ordinary linear response.

## Canonical submission path

```text
frozen strict spatial cohort
        |
CHELSA extraction + pooled within-species models (script 75)
        |
species variation and species-specific climate slopes (script 77)
        |
minimum-n sensitivity + four-quadrant plot + Tables S1/S2
        |
validated final artifact and durable submission bundle
```

Authoritative entry points:

- `analysis/README.md`
- `workflows/README.md`
- `manuscript/FINAL_STORY_AND_ANALYSIS.md`
- `manuscript/RUNBOOK.md`
- `REPOSITORY_REVIEW.md`
- `ch1_global/CH1_MANUSCRIPT_PLAN.md`

## Repository structure

```text
analysis/             submission-facing analysis definitions
workflows/            documentation for authoritative Actions workflows
manuscript/           final story, interpretation boundary and runbook
archive/              policy and destination for superseded material
ch1_global/v1/        frozen orientation-only baseline
ch1_global/v2/        current measurement and analysis implementation
ch1_shared/           shared collection and model utilities
ch1_japan/            regional exploratory work outside the Chapter 1 headline
tests/                deterministic unit and invariant tests
.github/workflows/    executable GitHub Actions definitions
```

The numbered v2 scripts remain in place for provenance. Mass renaming is deferred
until the final bundle records original paths, commits and checksums.

## Primary traits

- orientation angle relative to EXIF-oriented image vertical;
- corolla Lab lightness and chroma;
- circular hue sine/cosine;
- capitulum aspect ratio, circularity, solidity and width-profile variation.

High-resolution involucral projection and spine-like contours remain exploratory.
Hair and mucilage remain evidence-backed queues; missing evidence is `unknown`,
not biological absence.

## Reproducibility policy

- The final workflow rebuilds the climate table and pooled coefficients rather
  than silently combining artifacts from different commits.
- Primary species-level lability requires at least 10 observations, 6 eligible
  traits and 3 climate predictors per trait.
- Minimum sample sizes of 5 and 20 are reported as sensitivity analyses.
- Every submission bundle must contain CSVs, figures, metadata, software
  versions, commit SHA, workflow run ID and SHA-256 checksums.
- GitHub Actions artifacts are temporary working products; the frozen submission
  bundle should be archived durably, for example in Zenodo.

## Interpretation limits

Citizen-science photographs are not colour calibrated, image vertical is only a
proxy for gravity, and outline measures remain viewpoint dependent. Species
slopes may reflect geographic structure or unmeasured covariates. The grafted
mega-tree is used only as a historical sensitivity analysis because hybridization,
chloroplast capture and polyploid history are not resolved by a single
bifurcating tree.

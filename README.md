# azami — global image-derived capitulum traits in *Cirsium*

This repository contains a multi-chapter project on the ecology and evolution of
*Cirsium* floral architecture. Chapter 1 is the submission-focused component: a
global macroecological analysis of continuous capitulum traits measured from
public biodiversity photographs.

## Chapter 1 in one sentence

We quantify continuous corolla colour, capitulum outline and image-referenced head
orientation, model which photographs are measurable, compare climatic structure
among and within species, and separate total within-species variation from
species-specific climate tracking.

The analysis is observational and non-causal. Climate association is not proof of
local adaptation, phenotypic plasticity or selection.

## Executed datasets

The image-measurement layer contains 6,626 detected capitula from 3,725 primary
public observations and 216 accepted image-analysis taxa. The strict
within-species climate layer contains 46,276 spatially thinned observations from
259 taxa. These tables have different sampling units and are intentionally kept
separate.

The final two-axis lability analysis requires at least 10 observations, six
eligible traits and three climate predictors per trait. A total of 102 taxa pass
this completeness rule.

## Frozen result

Species within-variation and species environmental responsiveness are distinct
axes. Across the 102 complete taxa, their Spearman correlation is -0.333. Thus,
species with broad visible phenotypic variation are not necessarily those whose
traits most strongly track the four CHELSA gradients.

Median module summaries are:

| Module | Within-variation | Environmental responsiveness | Species |
|---|---:|---:|---:|
| Colour | 0.553 | 0.175 | 104 |
| Orientation | 0.630 | 0.129 | 101 |
| Shape | 0.665 | 0.153 | 103 |

Four pooled linear endpoint–predictor effects pass BH-FDR at 0.05, but all are
small. The active manuscript must not retain the older statement that no
within-species effect passed FDR.

## Scientific logic

Chapter 1 separates four levels that are often conflated:

1. **measurement assessability** — whether an uncontrolled photograph supports a
   valid trait estimate;
2. **among-species climatic structure** — whether species summaries covary with
   present climate;
3. **within-species climate association** — whether traits track environmental
   gradients within species;
4. **species lability decomposition** — whether a species is broadly variable,
   climate responsive, both or neither.

Circular hue remains a joint sine/cosine endpoint. Angle is used for presentation,
not as an ordinary linear response.

## Canonical submission path

```text
frozen strict spatial cohort
        |
pinned CHELSA table + pooled coefficients (script 75 provenance)
        |
species variation and species-specific slopes (script 77)
        |
minimum-n sensitivity + four-quadrant plot + Tables S1/S2
        |
validated Actions artifact + durable submission bundle
```

The validated lability workflow is `.github/workflows/ch1-phase1-lability.yml`.
The frozen successful run is `29382767192` at commit
`ccfb15dbce720b324dd7a121f3eda2b23a49448b`.

## Repository structure

```text
analysis/             submission-facing analysis definitions
workflows/            authoritative workflow documentation
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

## Authoritative documents

- `analysis/README.md`
- `workflows/README.md`
- `manuscript/FINAL_STORY_AND_ANALYSIS.md`
- `manuscript/RUNBOOK.md`
- `REPOSITORY_REVIEW.md`
- `ch1_global/CH1_MANUSCRIPT_PLAN.md`

## Primary traits

- orientation angle relative to EXIF-oriented image vertical;
- corolla Lab lightness and chroma;
- circular hue sine/cosine;
- capitulum aspect ratio, circularity, solidity and width-profile variation.

High-resolution involucral projection and spine-like contours remain exploratory.
Hair and mucilage remain evidence-backed queues; missing evidence is `unknown`,
not biological absence.

## Reproducibility policy

- Primary lability requires n >= 10, at least six traits and at least three
  climate predictors per trait.
- Minimum sample sizes 5 and 20 are reported as sensitivity analyses.
- The final bundle must contain CSVs, figures, metadata, package versions, commit
  SHA, workflow run ID and SHA-256 checksums.
- GitHub Actions artifacts are temporary; the submission bundle should be archived
  durably, for example in Zenodo.

## Interpretation limits

Citizen-science photographs are not colour calibrated, image vertical is only a
proxy for gravity and outline measures remain viewpoint dependent. Species slopes
may contain geographic structure or unmeasured covariates. The grafted mega-tree
is a historical sensitivity analysis, not a resolved evolutionary history.

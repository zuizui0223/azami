# azami — global image-derived capitulum traits in *Cirsium*

This repository contains a multi-chapter project on the ecology and evolution of
*Cirsium* floral architecture. Chapter 1 is the submission-focused component: a
global analysis of continuous capitulum traits measured from public biodiversity
photographs.

## Chapter 1 in one sentence

Species means conceal a scale-dependent and modular capitulum phenotype: most
visible variation occurs within species, but total variation, within-species climate
tracking and among-species environmental sorting are not equivalent.

The analysis is observational and non-causal. Climate association is not proof of
local adaptation, phenotypic plasticity, pollinator selection or evolutionary rate.

## Executed datasets

The image-measurement layer contains 6,626 detected capitula from 3,725 primary
public observations and 216 accepted image-analysis taxa. The strict within-species
climate layer contains 46,276 spatially thinned observations from 259 taxa. These
tables have different sampling units and are intentionally kept separate.

The final two-axis lability analysis requires at least 10 observations, six eligible
traits and three climate predictors per trait. A total of 102 taxa pass this
completeness rule.

## Frozen result

Across the nine primary endpoints, approximately 82–99% of visible variance occurs
within species. This is visible heterogeneity, not a claim that all variance is
biological, because viewpoint and measurement conditions also contribute.

Species within-variation and species environmental responsiveness are distinct axes.
Across the 102 complete taxa, their Spearman correlation is -0.333. Species with
broad visible phenotypic variation are therefore not necessarily those whose traits
most strongly track the four CHELSA gradients.

| Module | Median within-variation | Median environmental responsiveness | Species |
|---|---:|---:|---:|
| Colour | 0.553 | 0.175 | 104 |
| Orientation | 0.630 | 0.129 | 101 |
| Shape | 0.665 | 0.153 | 103 |

Colour has the strongest median climate tracking but the lowest median total
variation. Shape has the largest visible variation without the largest climate
responsiveness. Environmental niche contrasts between low- and high-trait species
are strongest for orientation and colour and weaker for most outline traits.

Two within-species summaries must remain distinct:

- the strict <=10 km primary cohort has zero BH-FDR-significant primary effects;
- the expanded pooled coefficient table has four BH-FDR-significant effects, all
  with small standardized magnitudes.

These results are not contradictory because they use different cohorts and
reporting scopes. The active manuscript must name the cohort whenever reporting an
FDR count.

## Scientific logic

Chapter 1 separates five levels that are often conflated:

1. **measurement assessability** — whether an uncontrolled photograph supports a
   valid trait estimate;
2. **visible within-species variation** — how broad the observed phenotype is;
3. **within-species climate association** — whether traits track environmental
   gradients within species;
4. **among-species environmental sorting** — whether species with contrasting mean
   phenotypes occupy different environments;
5. **historical sensitivity** — whether conclusions survive uncertain tree
   placement.

Circular hue remains a joint sine/cosine endpoint. Angle is used for presentation,
not as an ordinary linear response.

## Canonical manuscript path

Start here:

1. `manuscript/RESULT_LED_STORY_REVIEW.md` — full result inventory and competing story comparison;
2. `manuscript/FINAL_MANUSCRIPT_STRATEGY.md` — selected narrative and figure strategy;
3. `manuscript/final_claims.json` — machine-readable frozen numbers;
4. `manuscript/FINAL_FREEZE_CHECKLIST.md` — remaining submission gates;
5. `manuscript/RUNBOOK.md` — reproducible execution and bundle rules;
6. `ch1_global/SUBMISSION_READINESS.md` — current readiness verdict.

The validator `analysis/validate_final_claims.py` and workflow
`.github/workflows/ch1-final-manuscript-claims.yml` enforce the frozen claim set.

## Canonical analysis path

```text
public images + trait-specific assessability
        |
continuous orientation, colour and outline distributions
        |
within-species variance + species-specific environmental slopes
        |
module and among-species environmental-niche contrasts
        |
spatial/cohort and historical-placement sensitivity
        |
validated durable submission bundle
```

The validated lability workflow is `.github/workflows/ch1-phase1-lability.yml`.
The final manuscript-claim workflow is
`.github/workflows/ch1-final-manuscript-claims.yml`.

## Repository structure

```text
analysis/             submission-facing analysis definitions and validators
workflows/            authoritative workflow documentation
manuscript/           result review, final story, claim registry and runbook
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
Hair and mucilage remain evidence-backed queues; missing evidence is `unknown`, not
biological absence.

## Reproducibility policy

- Primary lability requires n >= 10, at least six traits and at least three climate
  predictors per trait.
- Minimum sample sizes 5 and 20 are reported as sensitivity analyses.
- The final bundle must contain CSVs, figures, metadata, package versions, commit
  SHA, workflow run IDs and SHA-256 checksums.
- GitHub Actions artifacts are temporary; the submission bundle must be archived
  durably before submission.

## Interpretation limits

Citizen-science photographs are not colour calibrated, image vertical is only a
proxy for gravity and outline measures remain viewpoint dependent. Species slopes
may contain geographic structure or unmeasured covariates. The grafted mega-tree is
a historical sensitivity analysis, not a resolved evolutionary history.

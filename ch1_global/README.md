# Chapter 1 — global continuous capitulum-trait macroecology

## Status

The Chapter 1 analysis is implemented and has been executed on the 216-taxon
global dataset. The scope is now frozen for manuscript completion: the paper
will be organized around the existing continuous traits, QC-retention analysis,
present-climate comparison and alternative-tree historical sensitivity.

MaxEnt, LGM hindcasting, palaeoclimate stability, a new genomic phylogeny, hair
and mucilage are not part of Chapter 1. They remain possible follow-up projects.

Manuscript-facing documents:

- [`CH1_MANUSCRIPT_PLAN.md`](CH1_MANUSCRIPT_PLAN.md) — fixed questions,
  structure, figures, exclusions and remaining work;
- [`CH1_NOVELTY_REVIEW.md`](CH1_NOVELTY_REVIEW.md) — provisional literature
  audit and defensible novelty wording;
- [`SUBMISSION_READINESS.md`](SUBMISSION_READINESS.md) — remaining submission
  gates;
- [`CODE_REVIEW_SUBMISSION.md`](CODE_REVIEW_SUBMISSION.md) — code and
  reproducibility review.

The frozen orientation-only baseline remains under `v1/`; the active submission
code is under `v2/`.

## Research question

Across global *Cirsium*, how do image-derived continuous capitulum traits vary
with climate and geography, and which associations remain stable after
measurement-QC, coordinate-precision and historical-tree sensitivity analyses?

The chapter reports macroecological pattern. It does not infer adaptation,
selection, pollinator causation or within-species plasticity from species-level
associations.

## Analysis hierarchy

### Primary endpoints

Nine continuous endpoints are used in the main analysis:

1. orientation angle relative to EXIF-oriented image vertical;
2. corolla Lab lightness;
3. corolla Lab chroma;
4. hue sine;
5. hue cosine;
6. capitulum aspect ratio;
7. capitulum circularity;
8. capitulum solidity;
9. capitulum width-profile coefficient of variation.

Hue is represented by sine and cosine because it is circular. Categories may be
used for visual summaries only and are not the inferential measurement scale.

### Auxiliary endpoints

Three high-resolution involucre/spine image proxies are retained as exploratory
supplements. They measure contour projection and protrusion, not botanical
states such as appressed, recurved, unarmed or long-spined.

### Evidence-only traits

Hair, mucilage, ploidy and hybrid status are stored in source-backed evidence
queues. All values begin as `unknown`; absent evidence is not converted to a
negative state.

## Executed sample and QC

The frozen run contains 6,626 detected heads, 3,725 observations and 216 taxa.
QC-retained head measurements were:

- colour: 5,777 heads (87.2%);
- shape: 5,324 heads (80.4%);
- orientation: 4,585 heads (69.2%).

QC retention is an assessability outcome, not accuracy. Failed measurements are
kept as missing and their environmental, taxonomic and spatial distribution is
audited before ecological modelling.

The high-resolution involucre selection contained 1,819 heads. After resolution,
sharpness, mask-quality and horizontal-flip repeatability checks, 1,443 heads
from 1,292 observations and 210 taxa remained usable.

## Canonical execution order

The current submission pipeline begins with the existing detector/crop metadata
and proceeds as follows:

| Stage | Scripts | Purpose |
|---|---|---|
| Continuous measurement | `52–56` | Measure and aggregate colour, outline and orientation |
| QC and environment | `57` | Audit differential assessability and join CHELSA |
| Involucre supplement | `58–60` | Select high-resolution heads, measure proxies, prepare evidence queues |
| Integration and pre-tree models | `61–64` | Build figures/tables and separate within- from between-species models |
| Historical sensitivity | `65–68` | Audit tree coverage and run bounded-lambda PGLS across alternative trees |
| Submission contract | `69–70` | Validate outputs and write checksummed provenance manifest |

Detailed commands and required inputs are in
[`v2/SUBMISSION_PIPELINE.md`](v2/SUBMISSION_PIPELINE.md).

## Statistical design

### QC-retention audit

For each primary trait, usability is modelled as an outcome against four CHELSA
predictors. The analysis tests whether the retained subset is environmentally
selective. These tests are measurement-bias diagnostics and are never
interpreted as biological trait effects.

### Within-species models

Observation-level outcomes and climate predictors are demeaned within species.
This removes all time-invariant species-level differences before fitting the
climate association. Species-clustered uncertainty is retained, and full versus
coordinate-precision <=10 km cohorts are reported separately.

### Between-species models

Species medians are analysed with:

- OLS as a non-phylogenetic reference;
- Brownian PGLS;
- Pagel-lambda PGLS with lambda constrained to [0, 1];
- deterministic grafting scenarios S1 and S3;
- 50 random within-genus grafting trees for the nine primary endpoints.

The vascular-plant backbone directly includes only 54 of 216 taxa. Therefore
these PGLS analyses are historical-constraint sensitivities, not a claim that
one bifurcating tree represents the true history of *Cirsium*.

## Current result boundary

Six primary endpoint–climate combinations passed the strict random-tree
stability screen. However, the coordinate-precision <=10 km within-species
analysis produced no FDR-significant primary effect. The submission should
therefore describe the stable results as between-species macroecological or
geographic associations, not local environmental response.

No auxiliary involucre/spine coefficient survives the deterministic PGLS FDR
screen. Auxiliary results remain hypothesis-generating.

## Provisional novelty boundary

The paper must not claim that automated plant-trait inference, floral
orientation ecology or PGLS are new methods. The likely contribution is their
integration in *Cirsium*: multiple continuous capitulum traits derived from
public photographs, explicit modelling of differential assessability,
within-versus-between species climate inference and sensitivity to uncertain
historical placement.

Until the structured Web of Science/Scopus/Google Scholar audit is completed,
the manuscript should use `to our knowledge` rather than an unconditional
priority claim.

## Retired paths

The earlier CLIP categorical states, representative human-holdout builders and
categorical forbidden-combination analysis are not part of the submission
pipeline. They remain available in Git history but have been removed from the
active tree to prevent accidental reuse.

The annotation and detector utilities retained under `v2/00–46` document how
images and head crops were assembled. They are provenance infrastructure, not
the current trait-inference endpoint.

## Workflow policy

Research-scale workflows are manual because they download public images and use
large intermediate artifacts. Pull requests use lightweight CI for:

- Python syntax;
- unit/invariant tests;
- configuration validation;
- stale-reference detection;
- submission-output contract tests.

Actions artifacts expire and are not citable data deposits. Final tables,
figures, configuration, file checksums and environment reports must be archived
with the manuscript release.

## Submission checklist

Before journal submission:

- complete residual spatial and leave-one-region-out sensitivity;
- complete detector and representative-measurement audits;
- freeze accepted taxon names and input-data checksums;
- complete the structured novelty search log;
- run the canonical submission validator;
- generate the provenance manifest on the final output bundle;
- deposit output tables and required intermediate metadata in a durable archive;
- record image licences and source identifiers;
- obtain a repository licence decision and add the final manuscript citation.

See [`v2/HISTORICAL_CONSTRAINTS.md`](v2/HISTORICAL_CONSTRAINTS.md) for the tree
interpretation boundary.

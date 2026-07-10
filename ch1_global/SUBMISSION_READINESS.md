# Chapter 1 — submission readiness

Review date: 2026-07-10  
Scope: global continuous capitulum traits in 216 *Cirsium* taxa.

## Current verdict

**The analytical codebase is ready to be frozen for manuscript production after
this cleanup. The paper is not yet ready to submit until the durable data,
measurement demonstration and final manuscript package are completed.**

The decisive earlier blocker — reliance on unvalidated categorical CLIP states —
has been removed from the headline analysis. The current inference uses
continuous deterministic image measurements with explicit QC, missingness and
horizontal-flip repeatability checks.

## Completed

- 6,626-head continuous measurement run completed.
- Observation and species aggregation completed.
- QC-retention bias audited by taxon, space and CHELSA environment.
- Primary traits joined to environmental metadata.
- High-resolution involucre/spine supplement completed.
- Hair, mucilage, ploidy and hybrid evidence queues created with `unknown` as the
  default state.
- Within-species and between-species models separated.
- Bounded Pagel-lambda and Brownian PGLS completed on deterministic and 50 random
  grafting trees.
- All 522 historical model fits succeeded.
- Tree provenance audited: 54/216 direct backbone tips; 162/216 grafted taxa.
- Obsolete categorical holdout and categorical integration implementations
  removed from the active submission tree.
- Submission-focused documentation, CI and output-contract tooling added.

## Claims currently supported

The manuscript may report:

- global between-species associations between continuous capitulum traits and
  climate;
- the subset stable across alternative within-genus grafting hypotheses;
- absence of FDR-significant primary effects in the strict <=10 km
  within-species analysis;
- differential image assessability as a quantified measurement limitation;
- exploratory high-resolution involucre/spine patterns.

The manuscript must not claim:

- local adaptation or phenotypic plasticity from species-level associations;
- pollinator causation;
- a resolved *Cirsium* species tree;
- botanical absence from an unmeasurable image;
- categorical bract, spine, hair or mucilage states from the continuous image
  proxies.

## Remaining gates before submission

### 1. Measurement demonstration

A full categorical human holdout is no longer part of the analysis design.
However, the supplement still needs a compact, independently inspectable panel
showing what each continuous variable measures:

- representative low, middle and high values;
- QC failures;
- horizontal-flip technical replicates;
- known failure modes such as camera tilt, occlusion and background leakage;
- comparison with a small set of published morphological descriptions or
  manually measured examples.

This is a face-validity and interpretation audit, not a training-label exercise.

### 2. Detector performance

Report precision, recall and missed-head/error categories for the frozen
capitulum detector on an independent source-image subset. The trait QC does not
replace detector evaluation because undetected heads never enter the trait
pipeline.

### 3. Taxonomic freeze

Review approximate, synonym and multiple-match records from the Open Tree audit.
Record accepted-name decisions in a committed table. Do not silently accept an
external synonym when it conflicts with the study taxonomy.

### 4. Spatial robustness

Add the final residual spatial diagnostic and region/block sensitivity summary
to the manuscript tables. This is especially important because stable
between-species associations may reflect geographic species turnover.

### 5. Durable archive and licences

GitHub Actions artifacts expire and are not citable. Deposit the final tables,
figures, input manifests, configuration and checksums in Zenodo or another
permanent repository. Include iNaturalist identifiers and licence fields without
redistributing images beyond their licence terms.

### 6. Manuscript package

Complete the manuscript, supplement, figure captions, data-availability
statement and code-availability statement. Add the final DOI and citation to the
repository release.

## Submission freeze rule

The final analysis release must pass `69_validate_submission_outputs.py` and
`70_build_submission_manifest.py` on the exact archived bundle. Any rerun that
changes a primary table must receive a new analysis version and manifest rather
than overwriting the submitted result.

## Bottom line

The repository now has a coherent publication analysis rather than competing
categorical and continuous endpoints. The remaining work is validation display,
detector reporting, spatial/taxonomic review, durable archiving and manuscript
writing — not another redesign of the core trait pipeline.

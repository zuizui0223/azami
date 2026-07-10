# Chapter 1 submission code review

Review date: 2026-07-10  
Reviewed scope: active continuous-trait pipeline, scripts 52–70, production
workflows, tests and submission-facing documentation.

## Executive assessment

The repository had a scientifically usable executed analysis, but the active
code tree still represented several abandoned designs. The largest publication
risk was not a numerical bug: it was that documentation, CI and old categorical
workflows described a different analysis from the one that generated the final
results.

This review removes those competing paths, creates one canonical pipeline and
adds a machine-verifiable output contract. The current numerical results are not
changed by the cleanup.

## High-priority findings resolved

### 1. Submission documentation contradicted the executed analysis

**Previous risk:** the Chapter 1 README described v2 as unimplemented, described
the chapter as phylogeny-free and excluded involucral/spine traits even though
continuous measurements, involucre proxies and PGLS sensitivities had already
been executed.

**Resolution:** root and Chapter 1 READMEs now describe the executed 216-taxon
continuous analysis, its interpretation boundary and the actual scripts used.
`SUBMISSION_PIPELINE.md` contains commands checked against each script's current
CLI.

### 2. Abandoned categorical validation remained in the active tree

**Previous risk:** scripts 47–51, their tests and a dedicated workflow suggested
that discrete CLIP classes and a representative human holdout were still the
submission gate. This conflicted with the adopted deterministic continuous
measurement design.

**Resolution:** the categorical holdout and categorical integration paths were
removed from the active tree. They remain recoverable from Git history. The old
CLIP states may be discussed as diagnostics but cannot enter headline analysis
tables.

### 3. Expensive workflows ran on pull requests and used stale artifact defaults

**Previous risk:** editing a script could trigger public-image reconstruction or
large tree generation, and workflow defaults pointed to temporary historical
Actions artifact IDs. A successful rerun could therefore depend on expired or
unrecorded inputs.

**Resolution:** research-scale workflows are now manual. Required artifact IDs
must be supplied explicitly; downloaded ZIP checksums, Git SHA and software
environments are written into the output. Pull requests use one lightweight
submission CI.

### 4. Output validity was encoded only in ad hoc workflow snippets

**Previous risk:** exact sample sizes, endpoint counts, lambda bounds and
random-tree completeness were checked in different workflows but were not a
portable contract that could validate a Zenodo/manuscript bundle.

**Resolution:** `submission_config.json` and
`69_validate_submission_outputs.py` define the frozen result contract. The
validator rejects old AI-state columns, holdout strata, model-weight files,
incomplete random-tree results, unbounded lambda values and unexpected table
dimensions.

### 5. No durable provenance manifest

**Previous risk:** Actions artifacts can expire and contained no single inventory
of file checksums, table dimensions, package versions and commit identity.

**Resolution:** `70_build_submission_manifest.py` writes SHA-256 hashes for every
file, CSV dimensions and columns, Python/package/platform information, Git SHA
and optional R `sessionInfo()`.

### 6. Canonical runbook commands did not match the current CLIs

**Previous risk:** argument names for scripts 53, 56 and 62 were outdated. A
reader following the documented commands would fail before reproducing the
analysis.

**Resolution:** commands were corrected against current `argparse` definitions,
and submission CI compiles all active scripts and validates the JSON contract.

## Real-output verification

The three executed result bundles were assembled without duplicate working
copies and checked against the submission contract. The candidate bundle passed:

- 6,626 detected heads;
- 3,725 unique observations;
- 216 unique taxa;
- nine primary and three auxiliary endpoints;
- 36 strict within-species endpoint–predictor rows;
- 54 direct GBOTB backbone tips;
- 50 random trees;
- 522 successful and zero failed historical model fits;
- 36/36 comparisons with all 50 tree fits;
- Pagel lambda range 0 to 0.8473;
- zero obsolete AI/holdout columns;
- zero detector/model-weight files.

## Intentionally retained technical debt

### Layered measurement engine (`52` → `55` → `56`)

The corrected v2 entrypoint dynamically loads the earlier measurement engine via
an OpenCV compatibility layer. This is not elegant, but it is the exact path
used for the frozen production result and is covered by tests. Rewriting it into
a clean Python package immediately before submission would change import and
numerical behaviour without adding scientific evidence.

**Decision:** freeze it for the submitted analysis. Refactor only on a new
analysis version with numerical regression tests against the frozen head-level
table.

### Thresholds are code constants rather than an external TOML/YAML file

Confidence, sharpness, mask and flip thresholds are explicit in scripts and
workflows but not yet consolidated in one configuration object.

**Decision:** the submission manifest and workflow record the current values.
Move them into a typed configuration file only if the analysis is rerun; doing so
requires a new analysis version.

### Actions artifact transport remains temporary

The improved workflows record artifact IDs and ZIP checksums, but Actions is
still transport rather than durable storage.

**Decision:** deposit the validated submission bundle and provenance directory in
a permanent repository before submission.

## Statistical and scientific items still required

These are not code-cleanup defects, but they remain manuscript gates.

1. **Residual spatial structure and region sensitivity.** Add residual Moran's I
   or an equivalent spatial diagnostic and leave-one-region/block sensitivity for
   headline between-species results.
2. **Taxonomic decision table.** Review approximate/synonym/multiple Open Tree
   matches and freeze accepted-name decisions.
3. **Detector audit.** Report precision, recall and missed-head categories on an
   independent source-image subset.
4. **Measurement demonstration.** Provide representative low/middle/high values,
   QC failures, mirrored technical replicates and known failure modes.
5. **Circular colour inference.** Hue sine and cosine are valid coordinates, but
   biological interpretation should use a joint circular statement rather than
   treating one component as a standalone colour direction.
6. **Repository metadata.** Add a licence only after ownership and third-party
   model/data conditions are reviewed. Add `CITATION.cff` after manuscript title,
   author order, ORCIDs and DOI are frozen.

## Submission rule

Do not alter the frozen output tables in place. Any change to taxon acceptance,
measurement thresholds, endpoint definitions, environmental data or tree
scenarios requires a new analysis version, full rerun, validation report and
checksummed manifest.

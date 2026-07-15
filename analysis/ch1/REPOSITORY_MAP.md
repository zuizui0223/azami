# Chapter 1 active and legacy code map

## Author-facing active layer

These are the files manuscript authors should use directly:

- `analysis/ch1/run_submission.py` — stable CLI for checks, lability, claims and manifest generation;
- `analysis/ch1/pipeline.json` — machine-readable mapping to exact implementation paths;
- `manuscript/MANUSCRIPT_OUTLINE.md` — drafting skeleton;
- `manuscript/FIGURE_TABLE_MAP.md` — output-to-manuscript mapping;
- `manuscript/final_claims.json` — frozen numerical claims;
- `manuscript/RUNBOOK.md` — release procedure.

## Frozen implementation layer

These paths generated or validate the frozen analysis and must remain stable:

- scripts 52–56: continuous image measurement and integration;
- scripts 59–64: continuous endpoints, within-species and between-species analyses;
- scripts 65–70: historical sensitivity, validation and submission manifest;
- script 75: strict CHELSA observation and pooled coefficient reconstruction;
- script 77: two-axis species lability and supplement;
- script 78: phylogenetic-signal sensitivity;
- script 79: live molecular database audit;
- script 80: manuscript-facing phylogenetic evidence summary.

They remain in `ch1_global/v2/` because workflows and frozen provenance refer to
those exact paths. Consolidation is achieved with wrappers and manifests, not by
renaming executed source immediately before submission.

## Supporting active tools

- detector and measurement audit applications;
- taxonomic and tree-placement audits;
- grouped SPDE-INLA sensitivity;
- workflow definitions required to reproduce frozen artifacts;
- tests enforcing the submission contract.

Supporting tools may be cited in Methods or Supplement but are not separate
headline analyses.

## Legacy or exploratory material

The following should not define manuscript claims:

- v1 orientation-only baselines;
- superseded categorical CLIP classifications;
- abandoned RF/XGBoost comparisons;
- unexecuted multivariate TMB experiments;
- exploratory hair, mucilage and high-resolution spine proxies;
- regional analyses outside the Chapter 1 sampling frame.

Legacy files remain recoverable through Git history or clearly labelled directories.
They must not be imported by the canonical runner or required by the final bundle.

## Cleanup rule

Delete or move an existing executed file only after:

1. the durable submission bundle records its path and SHA-256 checksum;
2. all workflows and documentation have been updated;
3. numerical regression tests show that the replacement is identical;
4. the cleanup is released as a new analysis version.

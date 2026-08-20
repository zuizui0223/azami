# Canonical Chapter 1 analysis entry point

This directory is the stable submission-facing interface for the frozen Chapter 1 analysis.

The scientific source of truth is `pipeline.json`; frozen numerical claims are in `manuscript/final_claims.json`. Executed numbered scripts under `ch1_global/v2/` are preserved at their original paths for provenance and are not renamed before durable release.

## Submission checks

```bash
python analysis/ch1/run_submission.py check
python analysis/ch1/run_submission.py claims
python analysis/ch1/run_submission.py summary
```

`summary` prints the active submission story, frozen cohort sizes, completed robustness gates and remaining external validation gates. It does not run the withdrawn raw lability/quadrant workflow.

## Active science

The submission pipeline consists of:

1. frozen continuous head-level trait extraction;
2. nested taxon → photograph → head visible-variance decomposition;
3. exhaustive spatially thinned within-taxon climate coefficients;
4. grouped SPDE-INLA spatial models;
5. high-resolution involucre contour analysis;
6. among-taxon niche permutation tests;
7. historical/PGLS sensitivity across alternative grafted trees;
8. residual-spatial, broad-region and WCVP synonym-collapse robustness;
9. independent detector and measurement validation gates.

The former negative variation–responsiveness relationship and median-split quadrants are statistical provenance only and must not be used as biological conclusions.

## Source of truth

- `pipeline.json` — exact stage map and rerun policy.
- `manuscript/final_claims.json` — frozen machine-readable claims.
- `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — cohort names, counts and permitted analyses.
- `manuscript/SUBMISSION_MANUSCRIPT.md` — current narrative and figure order.
- `manuscript/FIGURE_TABLE_MAP.md` — main/supplement mapping.
- `manuscript/RUNBOOK.md` — release procedure.

## Change policy

Changes to accepted taxa, trait definitions, thresholds, environmental layers or tree scenarios require a new analysis version and full rerun. File organization and prose may change only if frozen results and provenance remain unchanged.

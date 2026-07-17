# Chapter 1 spatial robustness protocol

## Purpose

The accepted grouped SPDE-INLA models already account for residual spatial structure. This gate does not refit those models or introduce a new headline analysis. It checks whether manuscript-facing conclusions remain credible when residual spatial pattern and broad-region dominance are inspected explicitly.

## Required diagnostic input

Prepare one CSV from the frozen model artifact with one row per observation and at least:

- `latitude`, `longitude`;
- `residual`: a clearly documented posterior or fitted-model residual;
- `endpoint`;
- `taxon_name`;
- `broad_region` from a reviewed region mapping;
- the manuscript-facing quantity used for region-omission ranking, defaulting to `residual` unless a separate `--value` is supplied.

The audit deliberately refuses to infer regions from arbitrary longitude cutoffs. Region labels must be supplied by the frozen artifact or a documented country-to-region mapping.

## Diagnostics

1. **Broad-region coverage**: rows and unique taxa per region, including each region's fraction of the dataset.
2. **Residual Moran's I**: endpoint-specific k-nearest-neighbour Moran statistics with permutation p-values.
3. **Leave-one-region-out rank stability**: Spearman agreement between complete-data taxon means and means after each broad region is omitted.

These are diagnostics, not replacement inferential models. A low permutation p-value indicates remaining spatial structure worth disclosing; it does not invalidate the SPDE model automatically. Low leave-one-region-out rank agreement indicates that a ranking depends strongly on one geographic region and must be qualified in Results or Discussion.

## Reproducible execution

Use `.github/workflows/ch1-spatial-robustness-audit.yml` with the exact source workflow run, artifact name and CSV path. The workflow records source and execution provenance and uploads:

- `broad_region_coverage.csv`;
- `residual_morans_i.csv`;
- `leave_one_region_out_rank_stability.csv`;
- `spatial_robustness_summary.json`;
- `workflow_provenance.json`.

## Submission gate

Before Methods and Figure 5 are frozen:

- residual type and extraction procedure must be named;
- all observations must have reviewed broad-region labels;
- endpoint-specific Moran results must be reported without treating significance as effect magnitude;
- any unstable leave-one-region-out result must be disclosed;
- run ID, artifact digest, commit and files must enter the durable submission bundle.

A diagnostic finding can change caveats and presentation. It should trigger a new core analysis version only if it reveals a material violation that requires changing the accepted cohort, spatial model, predictor set or interpretation boundary.

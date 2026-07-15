# Analysis

This directory is the submission-facing entry point for executed analyses.

## Canonical Chapter 1 analyses

- `ch1_global/v2/75_run_exhaustive_within_species_climate_analysis.py`
  rebuilds the strict observation-level CHELSA table and the global
  within-species coefficient table.
- `ch1_global/v2/77_build_lability_axes_and_supplement.py`
  produces the final two-axis lability analysis.
- `analysis/audit_cirsium_molecular_coverage.py`
  freezes taxon-level NCBI nucleotide, ITS, plastid and SRA record counts for the
  Chapter 1 taxon list. It is a database-coverage audit, not a phylogeny builder.

The scientific decision boundary for phylogenetic signal is documented in
`manuscript/PHYLOGENETIC_SIGNAL_FEASIBILITY.md`.

## Frozen primary definition

The lability analysis separates:

1. **Species within-variation index** — mean trait-standardized robust
   within-species dispersion, with circular dispersion for hue.
2. **Species environmental responsiveness index** — RMS of species-specific
   standardized climate slopes across eligible traits and predictors.

Primary inclusion requires at least 10 observations for a species–trait estimate,
at least 6 eligible traits and at least 3 climate predictors per trait. Minimum
sample sizes of 5 and 20 are reported as sensitivity analyses.

These indices describe image-derived variation and climate association. They are
not direct estimates of plasticity, local adaptation or evolutionary rate.

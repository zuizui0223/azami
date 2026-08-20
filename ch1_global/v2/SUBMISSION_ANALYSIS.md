# Chapter 1 submission analysis

## Scope

This directory contains the frozen executed implementation for the Chapter 1 manuscript on global continuous capitulum phenomics in thistles. New users should navigate through `analysis/ch1/`; exact numbered files in this directory are preserved because workflow and artifact provenance cite them.

## Current inferential structure

The submission is not organized around one lability index. Its current multiscale analysis is:

1. continuous head-level orientation, visible-colour and outline extraction;
2. nested assigned-taxon → photograph → head visible-variance decomposition;
3. exhaustive spatially thinned within-taxon CHELSA coefficients;
4. grouped per-endpoint SPDE-INLA models with climate, topography and soil predictor groups;
5. high-resolution involucre/spine-like contour proxies;
6. permutation-tested among-taxon environmental sorting;
7. alternative-tree PGLS historical sensitivity;
8. residual-spatial, broad-region and taxonomic-lumping robustness.

The former raw absolute-slope variation–responsiveness relationship and median-split quadrants were withdrawn after a precision audit showed strong dependence on slope sample size and standard error. Those products are retained only for statistical provenance.

## Grouped SPDE-INLA

The spatial model implementation remains:

- `79_sample_one_environment_layer.py` — samples environmental rasters for the frozen spatial cohort;
- `80_merge_parallel_environment_layers.py` — merges climate, topography and SoilGrids layers;
- `83_run_lightweight_spde_inla.R` — fits grouped per-trait hierarchical SPDE-INLA models;
- `.github/workflows/ch1-lightweight-spde-inla.yml` — reproduces the frozen grouped analysis.

For each endpoint, climate, climate + topography, climate + soil and full-environment groups use the same complete-case cohort. Predictors are centred within source-assigned taxa and standardized; models include a taxon iid random intercept and Matérn spatial field.

## Frozen provenance

The strict spatial/environment artifacts and grouped SPDE outputs remain pinned through the run IDs and digests recorded in `analysis/ch1/pipeline.json`, `manuscript/final_claims.json` and the result ledgers. Do not resample remote rasters or refit frozen models solely for repository cleanup.

## Interpretation rules

- All trait–environment relationships are observational spatial associations.
- Hue sine/cosine are interpreted jointly as a circular response.
- Image vertical is not gravity and camera-recorded colour is not calibrated reflectance.
- Involucre/spine-like variables are 2-D contour proxies, not direct botanical defence measurements.
- The grafted tree is a historical sensitivity device, not a resolved Cirsium phylogeny.
- No adaptive-radiation, plasticity, selection, pollinator-causation or evolutionary-rate claim follows from Chapter 1.

The active manuscript story and current completion gates are defined by `manuscript/SUBMISSION_MANUSCRIPT.md` and `analysis/ch1/pipeline.json`.

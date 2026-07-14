# Chapter 1 submission pipeline

The manuscript analysis is the grouped hierarchical SPDE-INLA workflow documented in `SUBMISSION_ANALYSIS.md`.

## Canonical files

- `79_sample_one_environment_layer.py`
- `80_merge_parallel_environment_layers.py`
- `83_run_lightweight_spde_inla.R`
- `.github/workflows/ch1-lightweight-spde-inla.yml`

## Canonical model set

Nine image-derived continuous traits are fitted separately under four predeclared predictor groups: climate, climate + topography, climate + soil, and full environment. Predictors are centred within species and standardized. Each model includes a species iid random intercept and a Matérn SPDE spatial field.

## Frozen-input reproduction

The workflow currently reuses the strict spatial cohort and successful environment-layer artifacts from Actions run `29306454759`. This avoids silent changes caused by resampling mutable remote raster services. When a permanently archived data release is created, replace the run dependency with that release and record its checksum here.

## Submission outputs

Use only the grouped SPDE-INLA outputs listed in `SUBMISSION_ANALYSIS.md` for manuscript tables, figures, and interpretation. Old one-off execution workflows and RF/XGBoost/TMB experiments are not part of the submitted analysis.

See `ANALYSIS_MANIFEST.tsv` for the active/removed file inventory.

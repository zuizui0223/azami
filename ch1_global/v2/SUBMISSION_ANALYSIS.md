# Chapter 1 submission analysis

## Scope

This directory contains the reproducible analysis used for the Chapter 1 manuscript on global image-derived capitulum traits in *Cirsium*.

The submission analysis is the lightweight grouped SPDE-INLA workflow. Experimental Random Forest, XGBoost, and multivariate TMB implementations are not part of the manuscript analysis and are retained only in Git history.

## Primary analysis

The primary inferential model is implemented in:

- `79_sample_one_environment_layer.py` — samples one environmental raster for the strict spatial cohort.
- `80_merge_parallel_environment_layers.py` — merges available climate, topography, and SoilGrids layers and builds partial-depth 0–30 cm soil composites.
- `83_run_lightweight_spde_inla.R` — fits the final grouped per-trait hierarchical SPDE-INLA models.
- `.github/workflows/ch1-lightweight-spde-inla.yml` — reproduces the analysis from the frozen environment-layer artifacts.

Nine continuous image-derived traits are analysed separately. Four predeclared predictor groups are compared:

1. climate;
2. climate + topography;
3. climate + soil;
4. full environment.

For each trait, the four groups use the same complete-case cohort. Environmental predictors are centred within species and globally standardized. Models include a species iid random intercept and a Matérn SPDE spatial field.

## Frozen data provenance

The strict spatial cohort and extracted environmental layers currently come from GitHub Actions run `29306454759`. The workflow downloads those frozen artifacts rather than resampling remote rasters. This preserves exact reproducibility while remote data services remain mutable.

The analysis tolerates unavailable raw layers. Missing layers are reported, and 0–30 cm soil composites are calculated from the available depth intervals with depth weighting.

## Main outputs

The workflow artifact `ch1-lightweight-per-trait-spde-inla` contains:

- `spde_model_group_summary.csv` — model status, sample sizes, WAIC, DIC, delta-WAIC, delta-DIC, and CPO diagnostics;
- `spde_fixed_effects_by_group.csv` — posterior fixed-effect summaries, 95% credible intervals, posterior tail-area screening values, and BH-adjusted values;
- `spde_effect_stability_across_groups.csv` — sign and credibility stability across predictor groups;
- `spde_hyperparameters_by_group.csv` — species and spatial hyperparameters;
- `spde_predictor_selection_by_group.csv` — predictor availability and collinearity filtering;
- `spde_run_metadata.csv` — mesh and run settings;
- environment coverage and missing-layer reports.

## Interpretation rules

The manuscript should treat SPDE-INLA as the inferential analysis. Effects should be emphasized only when their direction is stable across relevant model groups and uncertainty is adequately narrow. BH-adjusted posterior tail-area values are screening summaries, not frequentist p-values.

`corolla_hue_sin_median` and `corolla_hue_cos_median` are components of circular hue and should be interpreted jointly rather than as independent colour axes. Shape traits bounded near 0–1 require residual and sensitivity checks before strong biological interpretation.

## Reproduction

Run the workflow manually or open a pull request changing one of its tracked files. The workflow retrieves the frozen cohort and successful environment-layer artifacts, merges them, fits all 36 grouped models sequentially, and uploads compact results.

Do not use the deleted one-off execution workflows for manuscript reproduction. Their history documents exploratory and recovery runs only.

# Chapter 1 submission analysis

## Scope

This directory contains the reproducible analysis used for the Chapter 1 manuscript on global image-derived capitulum traits in *Cirsium*.

The submission analysis is the lightweight grouped SPDE-INLA workflow. Experimental Random Forest, XGBoost, multivariate TMB, and phylogenetic analyses are not part of the manuscript analysis and remain only in Git history.

## Primary analysis

The primary inferential model is implemented in:

- `79_sample_one_environment_layer.py` — samples one environmental raster for the strict spatial cohort.
- `80_merge_parallel_environment_layers.py` — merges available climate, topography, and SoilGrids layers and builds partial-depth 0–30 cm soil composites.
- `83_run_lightweight_spde_inla.R` — fits the final grouped per-trait hierarchical SPDE-INLA models.
- `.github/workflows/ch1-lightweight-spde-inla.yml` — reproduces the analysis from the frozen environment-layer artifacts.

Nine continuous image-derived traits are analysed separately. Four predeclared predictor groups are compared: climate, climate + topography, climate + soil, and full environment. For each trait, the four groups use the same complete-case cohort. Environmental predictors are centred within species and globally standardized. Models include a species iid random intercept and a Matérn SPDE spatial field.

## Frozen data provenance

The strict spatial cohort and extracted environmental layers currently come from GitHub Actions run `29306454759`. The grouped SPDE-INLA outputs reused by the submission extensions come from run `29339137476`. The workflows download those frozen artifacts rather than resampling remote rasters.

## Submission extensions

The lightweight extension workflow adds analyses that do not repeat raster extraction or SPDE fitting:

1. a priori orientation, colour, and shape module summaries;
2. species-level multivariate trait PCA;
3. descriptive environmental niche contrasts between low- and high-trait species quartiles;
4. a candidate phenotypic lability/conservatism screen.

The lability screen combines:

- the fraction of total trait variance occurring within species, with species-bootstrap intervals;
- median absolute standardized environmental effects from grouped SPDE-INLA;
- the fraction of environmental effects retained after global BH screening;
- consistency of effect direction across predictor groups.

This is a comparative screening index, not an evolutionary-rate estimate. Without a defensible species-level phylogeny, manuscript language must use **candidate phenotypic lability**, **environmental responsiveness**, or **candidate conservatism**, not demonstrated evolutionary lability or phylogenetic conservatism.

## Main outputs

The primary workflow artifact contains:

- `spde_model_group_summary.csv`;
- `spde_fixed_effects_by_group.csv`;
- `spde_effect_stability_across_groups.csv`;
- `spde_hyperparameters_by_group.csv`;
- `spde_predictor_selection_by_group.csv`;
- `spde_run_metadata.csv`;
- environment coverage and missing-layer reports.

The extension artifact contains:

- trait-module summaries and species trait PCA;
- environmental niche contrast tables;
- `trait_lability_conservatism_screen.csv`;
- `module_lability_conservatism_summary.csv`;
- `lability_screen_metadata.json`.

## Interpretation rules

SPDE-INLA remains the inferential analysis. Effects should be emphasized only when their direction is stable across relevant model groups and uncertainty is adequately narrow. BH-adjusted posterior tail-area values are screening summaries, not frequentist p-values.

`corolla_hue_sin_median` and `corolla_hue_cos_median` are components of circular hue and should be interpreted jointly. Shape traits bounded near 0–1 require residual and sensitivity checks before strong biological interpretation. The candidate lability score must not be treated as a substitute for phylogenetic comparative analysis.

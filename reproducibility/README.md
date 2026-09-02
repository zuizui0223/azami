# Reproducibility artifacts

This directory contains generated scientific artifacts that are intentionally versioned because they are part of the reproducibility record, not manuscript files.

## `figures/`

Frozen Main and Supporting result figures plus compact source/provenance files needed by the figure builders. These files were moved from the former manuscript-facing layout without changing their scientific content.

Canonical rebuild entry points:

```bash
python analysis/figures/run_image_to_trait_perturbation_audit.py
python analysis/figures/build_v2_figures.py
```

The Figure 2 disclosure-safe source table can be rebuilt from the frozen observation-level cohort when that external cohort is available:

```bash
python analysis/figures/build_v2_figure_inputs.py --source-cohort /path/to/strict_spatial_chelsa.csv
```

The generated outputs are written to `reproducibility/figures/`; they are intentionally tracked.

## `contracts/`

Compact trait-definition material required to reconstruct image-measurement and analysis decisions.

## Canonical numerical outputs

Frozen CSV/JSON result tables and validation reports remain in `analysis_outputs/`. Analysis contracts and predictor/trait registries are retained under `analysis/ch1/` and `ch1_global/v2/ontology/`.

Manuscript prose, reviewer material and journal submission packages are kept outside this repository.

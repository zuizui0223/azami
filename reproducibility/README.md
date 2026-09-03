# Reproducibility artifacts

This directory contains generated scientific artifacts that are intentionally versioned because they are part of the reproducibility record, not manuscript files.

## Numerical analysis rebuild

The canonical durable numerical entry point is:

```bash
python -m pip install -e '.[full,test]'
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip /path/to/artifact-9612943217.zip \
  --multilevel-zip /path/to/artifact-9632715852.zip \
  --out-dir rebuild_check \
  --mode validate_sources
```

The runner verifies the archived ZIP digests, the frozen long-trait table, cohort sizes and the complete-18 reference environment before doing any analysis. It accepts local ZIP files deliberately so the reconstruction path remains usable after GitHub Actions retention expires.

To reconstruct the nine-predictor environment and rerun the complete 27-endpoint atlas with fail-closed validation:

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip /path/to/artifact-9612943217.zip \
  --multilevel-zip /path/to/artifact-9632715852.zip \
  --out-dir rebuild_check \
  --mode rebuild_full27
```

To rerun the full-27 atlas plus the whole-capitulum secondary analyses:

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip /path/to/artifact-9612943217.zip \
  --multilevel-zip /path/to/artifact-9632715852.zip \
  --out-dir rebuild_check \
  --mode rebuild_all
```

`--raster-cache-dir` may point to exact URL-basename copies of the frozen CHELSA rasters. Without it, the process variables are sampled from their frozen remote URLs. The rebuilt 46,276-row nine-predictor environment must match the frozen SHA-256 or the runner stops.

The historical-placement input is also artifact-locked. Artifact `8227254443` contains S1, S3 and the 50 randomized S2 trees. It can be supplied directly as a local archived ZIP:

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip /path/to/artifact-9612943217.zip \
  --multilevel-zip /path/to/artifact-9632715852.zip \
  --historical-zip /path/to/artifact-8227254443.zip \
  --out-dir rebuild_check \
  --mode rebuild_full27 \
  --with-historical
```

The historical ZIP and all three tree resources are checked against frozen SHA-256 values before the 52-tree analysis starts. A pre-extracted exact tree directory may be supplied with `--tree-dir` instead.

Sampling-composition sensitivity still requires the exact frozen region and native-status auxiliary tables (`--with-sampling --regions ... --native-status ...`). Their expected SHA-256 values are retained in the frozen sampling-composition report; the runner does not silently reconstruct or substitute them.

The current stage inventory is `ch1_global/v2/ANALYSIS_MANIFEST.tsv`; CI checks that every repository path marked there actually exists.

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

Source artifact identities, digests, expiry metadata and verified archival state are recorded in `actions_artifact_catalog.json`. The immutable recovery tag `azami-ch1-v2-2026-08-27` retains retired workflow definitions without requiring manuscript or submission material to be restored to `main`.

Manuscript prose, reviewer material and journal submission packages are kept outside this repository.

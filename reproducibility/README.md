# Reproducibility artifacts

This directory contains generated scientific artifacts that are intentionally versioned because they are part of the reproducibility record, not manuscript files.

## Numerical analysis rebuild

For byte-identical numerical reconstruction, first recreate the package environment recovered from the frozen 2026-08-26/27 numerical and spatial artifacts:

```bash
python -m pip install -r reproducibility/frozen-numerical-rebuild-requirements.txt
python -m pip install -e . --no-deps
```

`python -m pip install -e '.[full,test]'` remains the looser developer/CI convenience environment. The frozen requirements file is the preferred route when exact archived output hashes are being reproduced.

The canonical durable numerical entry point is:

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip /path/to/artifact-9612943217.zip \
  --multilevel-zip /path/to/artifact-9632715852.zip \
  --out-dir rebuild_check \
  --mode validate_sources
```

The runner verifies the archived ZIP digests, frozen long-trait table, exact 46,276-row strict-spatial source bytes, cohort sizes and complete-18 reference environment before doing any analysis. It accepts local ZIP files deliberately so the reconstruction path remains usable after GitHub Actions retention expires.

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

## Sampling-composition inputs

The sampling audit has two auxiliary inputs and now fails closed on both before execution:

- broad-region lookup SHA-256 `085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff`;
- native-status table SHA-256 `c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a`.

The exact broad-region lookup was recovered during the 2026-09-03 audit from workflow run `31152400475`, artifact `8983877726`, inner path `spatial_regions/broad_region_lookup.csv`. Its artifact ZIP digest and inner-file digest are recorded in `actions_artifact_catalog.json`.

A fail-closed reconstruction path is also retained. Extract `environment/strict_spatial_chelsa.csv` from artifact `9612943217`, and obtain the Natural Earth 1:50m admin-0 archive from `https://naturalearth.s3.amazonaws.com/50m_cultural/ne_50m_admin_0_countries.zip`. The frozen archive SHA-256 recovered from the spatial artifact is `5fed433373581fa648920435f937d95f2d3c0200e067409c6478dcdf1b853139`. Then run:

```bash
python analysis/rebuild_frozen_broad_region_lookup.py \
  --observations /path/to/strict_spatial_chelsa.csv \
  --naturalearth-zip /path/to/ne_50m_admin_0_countries.zip \
  --output broad_region_lookup.csv \
  --report broad_region_rebuild.json
```

The wrapper verifies the strict-spatial cohort SHA, Natural Earth archive SHA, 46,276 rows, the frozen six-region plus `UNMAPPED` counts, assignment-method counts and final CSV SHA. The lower-level deterministic polygon assignment implementation remains `analysis/build_naturalearth_broad_region_lookup.py`.

The exact native-status source remains in the immutable recovery tag and can be recovered without manuscript files:

```bash
git show azami-ch1-v2-2026-08-27:analysis_outputs/native_range_sensitivity_v1/observation_native_status.csv > observation_native_status.csv
sha256sum observation_native_status.csv
# c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a
```

A manuscript-independent reconstruction path is also retained. Using the same frozen `environment/strict_spatial_chelsa.csv`, run:

```bash
python analysis/rebuild_frozen_native_status.py \
  --observation /path/to/strict_spatial_chelsa.csv \
  --contract analysis/ch1/native_range_sensitivity_contract.json \
  --out-csv observation_native_status.csv \
  --report native_status_rebuild.json
```

The native-status builder checks the frozen WCVP dataset identity, pinned TDWG commit and GeoJSON digest, 46,276-row/259-taxon cohort identity, the five frozen status counts, and the final output SHA. Any mismatch stops the reconstruction rather than silently substituting a newer classification.

With exact auxiliary tables present, sampling can be rerun through the canonical runner:

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip /path/to/artifact-9612943217.zip \
  --multilevel-zip /path/to/artifact-9632715852.zip \
  --out-dir rebuild_check \
  --mode rebuild_full27 \
  --with-sampling \
  --regions /path/to/broad_region_lookup.csv \
  --native-status /path/to/observation_native_status.csv
```

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

Source artifact identities, digests, expiry metadata and verified archival state are recorded in `actions_artifact_catalog.json`. The immutable recovery tag `azami-ch1-v2-2026-08-27` retains retired workflow definitions and the exact frozen native-status table without requiring manuscript or submission material to be restored to `main`.

Manuscript prose, reviewer material and journal submission packages are kept outside this repository.

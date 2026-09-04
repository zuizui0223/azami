# Chapter 1 v2 third-party reproducibility runbook

This is the canonical procedure for an **independent reader, reviewer, or researcher** to reproduce the frozen Chapter 1 v2 numerical analysis. It must not depend on the author's personal computer, OneDrive, private Google Drive, or GitHub Actions retention.

The machine-readable contracts are:

- `reproducibility/public_release_manifest.json` — what must be publicly downloadable for an independent rerun;
- `reproducibility/material_availability.json` — what currently exists in GitHub, owner preservation storage, Git history, and public external sources;
- `ch1_global/v2/ANALYSIS_MANIFEST.tsv` — live analysis stages and entry points.

## Current external reproducibility status

**Code and frozen outputs are public, but the full numerical input package has not yet been published in a public data archive.**

The required analysis ZIPs are safely preserved in owner-controlled storage and SHA-verified, but that private preservation copy is **not** a valid third-party distribution channel. Until a public DOI/URL is filled in `reproducibility/public_release_manifest.json`, an independent reader can audit the code, frozen outputs, validation reports, and figures, but cannot perform the complete numerical rerun without obtaining the unpublished input package from the author.

Do not describe the repository as fully independently reproducible until that public-data gate is closed.

## Reproducibility levels

### Level A — audit the reported results from GitHub

Available now to anyone:

- analysis and validation code;
- frozen trait and environmental contracts;
- exact Python package versions;
- frozen full-27 result tables;
- sampling, spatial, and historical sensitivity outputs;
- figure source tables and builders;
- validation reports.

Expected frozen states are:

```text
environment atlas:                PASS 24/24
sampling composition:             PASS 7/7
spatial + historical sensitivity: PASS 11/11
figures:                          status ok; 5 main + 7 supporting
```

### Level B — independently rerun the reported numerical analysis

This is the target publication-level reproducibility contract. It requires a **public, credential-free data archive** containing the files listed in `reproducibility/public_release_manifest.json` with the recorded SHA-256 values.

At minimum the public archive must provide:

```text
artifact-9612943217-continuous.zip
artifact-9632715852-multilevel.zip
artifact-8227254443-historical.zip
artifact-8983877726-spatial.zip
```

The exact native-status table is already recoverable from the public immutable Git tag `azami-ch1-v2-2026-08-27`.

Five added CHELSA v2.1 process rasters (rsds, VPD, wind, GSP, NPP) are identified by frozen public URLs in `analysis/ch1/chelsa_process_environment_sources.json`. If the public data release also archives exact copies of those rasters, the complete numerical rerun becomes offline-capable; otherwise this one stage requires network access to CHELSA.

### Level C — regenerate measurements from original public photographs

This is a stronger upstream reproduction target than the frozen numerical rerun. It additionally requires the detector package, source-image identities/URLs, crop metadata, and image-measurement pipeline. Those upstream materials are preserved, but the publication-level numerical contract above does not require re-downloading every original citizen-science photograph.

## Procedure once the public data release is published

The steps below become executable when `reproducibility/public_release_manifest.json` contains a public DOI/URL and immutable code ref.

### 1. Clone the public code

```bash
git clone https://github.com/zuizui0223/azami.git
cd azami
```

Checkout the immutable `code_ref` recorded in `reproducibility/public_release_manifest.json` rather than an arbitrary future `main` commit.

### 2. Download the public analysis input package

Download the files from the DOI/URL recorded in `reproducibility/public_release_manifest.json` into a directory such as:

```text
data_release/
├── artifact-9612943217-continuous.zip
├── artifact-9632715852-multilevel.zip
├── artifact-8227254443-historical.zip
└── artifact-8983877726-spatial.zip
```

No author credentials should be required.

### 3. Verify the downloaded bytes

```bash
python reproducibility/verify_local_materials.py \
  --archive-dir data_release \
  --scope numerical \
  --report local_material_check.json
```

Expected result:

```text
"status": "PASS"
```

The canonical ZIP SHA-256 values are:

```text
continuous  101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e
multilevel  51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6
historical  499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc
spatial     151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca
```

If any digest differs, stop. Do not substitute a newer or approximately equivalent dataset.

### 4. Recreate the frozen Python environment

Python 3.11 is the CI reference environment.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r reproducibility/frozen-numerical-rebuild-requirements.txt
python -m pip install -e . --no-deps
```

### 5. Validate the canonical direct inputs

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip data_release/artifact-9612943217-continuous.zip \
  --multilevel-zip data_release/artifact-9632715852-multilevel.zip \
  --out-dir rebuild_check \
  --mode validate_sources
```

This fails closed unless the runner recovers the frozen trait universe, 46,276-row / 259-taxon strict-spatial cohort, and complete-18 reference environment with their exact SHA values.

### 6. Rebuild the primary full-27 atlas

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip data_release/artifact-9612943217-continuous.zip \
  --multilevel-zip data_release/artifact-9632715852-multilevel.zip \
  --out-dir rebuild_full27 \
  --mode rebuild_full27
```

The five frozen CHELSA process layers are sampled from the URLs declared in `analysis/ch1/chelsa_process_environment_sources.json`. A public/local raster cache can instead be supplied with `--raster-cache-dir`.

The reconstructed 46,276-row nine-predictor environment must have SHA-256:

```text
e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7
```

The full-27 validator must finish at `PASS 24/24`.

### 7. Recover the exact sampling-composition auxiliary inputs

Extract the broad-region lookup from the public spatial ZIP:

```bash
mkdir -p rebuild_inputs
python -c "import zipfile,pathlib; z=zipfile.ZipFile('data_release/artifact-8983877726-spatial.zip'); pathlib.Path('rebuild_inputs/broad_region_lookup.csv').write_bytes(z.read('spatial_regions/broad_region_lookup.csv'))"
```

Expected SHA-256:

```text
085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff
```

Recover the native-status table byte-for-byte from the immutable public Git tag:

```bash
python -c "import subprocess,pathlib; pathlib.Path('rebuild_inputs/observation_native_status.csv').write_bytes(subprocess.check_output(['git','show','azami-ch1-v2-2026-08-27:analysis_outputs/native_range_sensitivity_v1/observation_native_status.csv']))"
```

Expected SHA-256:

```text
c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a
```

Both tables also have fail-closed reconstruction code in the repository.

### 8. Rerun the complete frozen analysis family

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip data_release/artifact-9612943217-continuous.zip \
  --multilevel-zip data_release/artifact-9632715852-multilevel.zip \
  --historical-zip data_release/artifact-8227254443-historical.zip \
  --out-dir rebuild_everything \
  --mode rebuild_all \
  --with-sampling \
  --regions rebuild_inputs/broad_region_lookup.csv \
  --native-status rebuild_inputs/observation_native_status.csv \
  --with-spatial \
  --with-historical
```

Expected frozen validation states:

```text
environment atlas:                PASS 24/24
sampling composition:             PASS 7/7
spatial + historical sensitivity: PASS 11/11
```

The historical ZIP is also hash-locked internally and must yield S1, S3, and 50 randomized S2 trees (52 placements total).

### 9. Rebuild the figures

```bash
python analysis/figures/run_image_to_trait_perturbation_audit.py
python analysis/figures/build_v2_figures.py
```

Expected `reproducibility/figures/figure_build_report.json`:

```text
status: ok
main figures: 5
supporting figures: 7
```

### 10. Run repository integrity tests

```bash
python -m compileall -q analysis src ch1_global/v2 reproducibility/verify_local_materials.py
python -m pytest -q \
  tests/test_ch1_canonical_manifest.py \
  tests/test_reproducibility_manifest.py \
  tests/test_material_availability.py \
  tests/test_azami_ch1_provenance.py \
  tests/test_azami_ch1_tabular.py
```

## What counts as a successful independent reproduction?

A third party should be able to start from only:

1. the public Git repository at the immutable code ref;
2. the public DOI/data archive;
3. public CHELSA v2.1 sources, unless the exact raster cache is included in the data archive.

They should not need author credentials, a private cloud folder, a local path from the author's machine, or unexpired GitHub Actions artifacts.

Success means all archived input hashes match and the rebuilt validation reports reproduce `24/24`, `7/7`, and `11/11`, with the frozen figure report returning `status: ok`.

## Canonical references

- public-release gate: `reproducibility/public_release_manifest.json`;
- material audit: `reproducibility/material_availability.json`;
- numerical runner: `analysis/rebuild_frozen_analysis.py`;
- stage inventory: `ch1_global/v2/ANALYSIS_MANIFEST.tsv`;
- exact package environment: `reproducibility/frozen-numerical-rebuild-requirements.txt`;
- frozen outputs: `analysis_outputs/v2_full27_environment_atlas_2026-08-27/`;
- CHELSA process source registry: `analysis/ch1/chelsa_process_environment_sources.json`.

# azami analysis and reproducibility repository

This repository retains the code and frozen machine-readable products needed to reproduce and audit the image-phenomics analyses.

Prepublication manuscript prose and journal-submission packages are intentionally kept outside GitHub.

## Canonical Chapter 1 analysis path

The current paper uses this evidence chain. Superseded categorical head-direction classifiers are historical development only and are not manuscript evidence.

```text
public photographs
  ↓
frozen YOLO11n visible-capitulum detector
  artifact 8076736948
  weight SHA-256 4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed
  ↓
exhaustive detector-positive image/head layer
  workflow run 29216617585
  artifact 8269246732
  ↓
category-free continuous trait reconstruction and measurement
  workflow run 32975451732
  artifact 9612943217
  ↓
within-taxon / among-taxon environmental analyses and whole-capitulum synthesis
  workflow run 33035785120
  artifact 9632715852
  ↓
v3 biological-construct reaggregation and robustness chain
  workflow run 34418904597
  artifact 10130210432
  ↓
v3 common-cohort/module upgrade
  workflow run 34432661967
  artifact 10135139012
  ↓
v3 assessability and technical-stress audit
  workflow run 34434467062
  artifact 10135679053
  ↓
v3 direct scale contrast
  workflow run 34435743710
  artifact 10136229131
```

The YOLO model is used to localize visible capitula. It is not the trait classifier. The current trait evidence is generated from detector-defined image regions by the frozen category-free continuous measurement implementation. In particular, the retired MobileNet `upward`/`nodding` classifier and its species-description-derived labels are not part of the current evidence chain.

The v2 implementation registry is [`ch1_global/v2/ANALYSIS_MANIFEST.tsv`](ch1_global/v2/ANALYSIS_MANIFEST.tsv). The frozen Actions/artifact identities are recorded in [`reproducibility/actions_artifact_catalog.json`](reproducibility/actions_artifact_catalog.json). The current v3 synthesis and result receipts are indexed in [`analysis/v3/README.md`](analysis/v3/README.md).

## Active GitHub Actions boundary

The submission/reproducibility surface should contain only workflows that either reproduce a current v3 result or continuously verify the frozen public analysis.

Current v3 result workflows:

- `.github/workflows/ch1-v3-biological-axes.yml` — biological-construct reanalysis and full frozen-v2 robustness chain.
- `.github/workflows/ch1-v3-scale-integration.yml` — pairwise within/among construct integration.
- `.github/workflows/ch1-v3-construct-upgrade.yml` — exact complete-18 common-cohort, bootstrap and module-cohesion upgrade.
- `.github/workflows/ch1-v3-assessability-measurement-stress.yml` — outcome-blind assessability and orientation technical-stress audit.

Persistent public reproducibility workflows:

- `.github/workflows/reproducibility-integrity.yml` — checks frozen contracts, provenance, repository layout and live entry points.
- `.github/workflows/rebuild-frozen-figures.yml` — rebuilds and verifies the frozen v2 figure set in scratch space.

One-off reviewer-defense workflows that do not feed the current v3 synthesis are intentionally removed from the active branch. Their commits and outputs remain recoverable from Git history. Earlier detector bootstrap/recovery, categorical AI/classifier, queue-building and exploratory submission workflows remain historical provenance only and are preserved through immutable tag `azami-ch1-v2-2026-08-27` and branch `archive/ch1-precleanup-20260827` rather than restored to the active submission surface.

## Quick Chapter 1 v3 reproduction

Chapter 1 v3 is a higher-level synthesis layered on the frozen v2 measurements. It does **not** retrain the detector, redownload photographs, or rerun the retired categorical classifier. The v3 workflows begin from frozen v2 trait/environment artifacts and reproduce the construct-level analyses, scale comparisons and technical/selection audits.

Unlike the frozen v2 public bundle below, v3 is currently the active PR #93 analysis branch rather than a separately frozen public release. For an exact reproduction of the current v3 state, use the branch and the byte-verified source artifacts listed here.

### 1. Check out the current v3 analysis branch

```bash
git clone https://github.com/zuizui0223/azami.git
cd azami
git checkout analysis/ch1-v3-biological-axes-20260910
```

Use Python 3.12. The four active v3 workflows pin the numerical runtime to:

```bash
python -m pip install \
  pandas==3.0.5 \
  numpy==2.4.6 \
  scipy==1.17.1 \
  statsmodels==0.14.6 \
  'biopython>=1.83,<2'
```

Biopython is required only for the full biological-axis historical-sensitivity chain; the other three v3 workflows use the first four packages.

### 2. Recover and verify the frozen v2 inputs

The two common v3 inputs are:

```text
continuous trait universe
  artifact 9612943217
  ZIP SHA-256 101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e
  file universe/continuous_trait_universe_observation_long.csv
  file SHA-256 d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344

nine-predictor process environment
  artifact 9633419268
  ZIP SHA-256 d7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4
  file process_environment/strict_spatial_chelsa_process.csv
  file SHA-256 f86a3418e9b21453026bba1aaf350b061f03048949547e2092602620af98cbf6
```

With GitHub CLI authentication that can read this repository's Actions artifacts:

```bash
mkdir -p work/continuous work/environment

gh api repos/zuizui0223/azami/actions/artifacts/9612943217/zip > work/continuous.zip
echo '101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e  work/continuous.zip' | sha256sum -c -
unzip -q work/continuous.zip 'universe/continuous_trait_universe_observation_long.csv' -d work/continuous
echo 'd775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344  work/continuous/universe/continuous_trait_universe_observation_long.csv' | sha256sum -c -

gh api repos/zuizui0223/azami/actions/artifacts/9633419268/zip > work/environment.zip
echo 'd7c0c466f55b67695d06ae46c21a6452dbe6cfd92a52db8042caa200429e97f4  work/environment.zip' | sha256sum -c -
unzip -q work/environment.zip 'process_environment/strict_spatial_chelsa_process.csv' -d work/environment
echo 'f86a3418e9b21453026bba1aaf350b061f03048949547e2092602620af98cbf6  work/environment/process_environment/strict_spatial_chelsa_process.csv' | sha256sum -c -
```

`9633419268` is the frozen nine-predictor environment **input** used by v3. Do not confuse it with artifact `9632715852`, which is the earlier v2 within/among and complete-18 process-analysis **output** recorded in the canonical evidence chain.

### 3. Reproduce the biological constructs and full robustness chain

This is the main v3 workflow corresponding to run `34418904597`, artifact `10130210432`.

First rebuild the construct-level environmental atlas:

```bash
mkdir -p work/results
python -m analysis.v3.run_biological_axis_reanalysis \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --out-dir work/results \
  --permutations 9999 \
  --seed 20260910
```

For the full sampling → spatial/residual → 52-tree historical chain, also recover the frozen spatial lookup and historical trees:

```text
broad-region lookup
  artifact 8983877726
  ZIP SHA-256 151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca
  broad_region_lookup.csv SHA-256 085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff

historical-placement trees
  artifact 8227254443
  ZIP SHA-256 499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc
```

The frozen native-range status is recovered from immutable tag `azami-ch1-v2-2026-08-27` and must match SHA-256 `c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a` after LF/CRLF transport normalization.

Then run:

```bash
python -m analysis.v3.run_biological_axis_sensitivity_chain \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --axis-among work/results/biological_axes_among_min5.csv \
  --axis-within work/results/biological_axes_within.csv \
  --regions work/sampling/broad_region_lookup.csv \
  --native-status work/sampling/observation_native_status.csv \
  --tree-dir work/historical/historical_trees \
  --out-dir work/results/sensitivity_chain \
  --spatial-permutations 999 \
  --moran-permutations 999 \
  --minimum-taxa-historical 30 \
  --seed 20260910
```

The key manuscript anchors should be recovered as:

```text
floral chroma × shortwave radiation
  beta = -0.345372
  q = 0.00225
  historical placement: 52/52 trees P < 0.05

presentation angle × annual precipitation (BIO12)
  beta = +0.304359
  q = 0.00360
  historical placement: 52/52 trees P < 0.05
```

### 4. Reproduce scale-dependent construct integration

The pairwise construct analysis corresponding to run `34421253852`, artifact `10131007603`, is rebuilt from the same trait and environment inputs:

```bash
mkdir -p work/axis work/integration

python -m analysis.v3.run_biological_axis_reanalysis \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --out-dir work/axis \
  --permutations 9999 \
  --seed 20260910

python -m analysis.v3.run_construct_scale_integration_entrypoint \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --axis-among work/axis/biological_axes_among_min5.csv \
  --axis-within work/axis/biological_axes_within.csv \
  --out-dir work/integration \
  --minimum-paired-observations-per-taxon 5 \
  --minimum-taxa 20 \
  --qap-permutations 9999 \
  --seed 20260910
```

Expected pairwise-cohort result:

```text
within/among integration-matrix Spearman rho = 0.374775
QAP one-sided P = 0.0131
```

### 5. Reproduce the complete-18 common-cohort/module upgrade

The common-cohort upgrade corresponding to run `34432661967`, artifact `10135139012`, forces all nine constructs and all 36 construct pairs onto the same frozen complete-18 cohort.

```bash
mkdir -p work/upgrade
python -m analysis.v3.run_construct_scale_upgrade \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --axis-among work/axis/biological_axes_among_min5.csv \
  --axis-within work/axis/biological_axes_within.csv \
  --out-dir work/upgrade \
  --minimum-complete-observations-per-taxon 5 \
  --bootstrap-replicates 1000 \
  --permutations 9999 \
  --seed 20260910

python -m analysis.v3.run_construct_scale_contrast_summary \
  --pairwise work/upgrade/complete18_construct_pairwise.csv \
  --bootstrap work/upgrade/complete18_taxon_bootstrap.csv \
  --out work/upgrade/construct_scale_contrast_summary.json
```

Expected checks:

```text
complete cohort = 1,734 observations / 42 taxa
within-among matrix rho = 0.439125
QAP P = 0.0041
bootstrap median rho = 0.31210
bootstrap 95% interval = 0.01708–0.53517
positive bootstrap replicates = 97.9%
module-label permutation P: within = 0.0013; among = 0.0365
observed median RV: within = 0.002238; among = 0.043212
33/36 construct relations stronger among taxa
```

The direct scale-contrast summary was rerun as run `34435743710`, artifact `10136229131`; its 1,000-bootstrap among-minus-within median-RV difference should be positive in 100% of replicates, with median `+0.065928` and 95% interval `+0.038197` to `+0.119048`.

### 6. Reproduce assessability and technical-stress audits

The current audit corresponds to run `34434467062`, artifact `10135679053`.

```bash
mkdir -p work/results/assessability work/results/technical_stress

python -m analysis.v3.run_assessability_selection_audit \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --out-dir work/results/assessability

python -m analysis.v3.run_frozen_technical_error_stress \
  --traits work/continuous/universe/continuous_trait_universe_observation_long.csv \
  --environment work/environment/process_environment/strict_spatial_chelsa_process.csv \
  --technical-audit-summary analysis/ch1/image_to_trait_automated_technical_audit_summary.json \
  --out-dir work/results/technical_stress \
  --replicates 2000 \
  --seed 20260910
```

Expected audit checks:

```text
endpoint assessability tests = 198; BH-supported = 17
construct assessability tests = 99; BH-supported = 11
largest absolute availability gradient ≈ 1.024 percentage points per within-taxon environmental SD

orientation reconstruction beta = +0.304359
mirror-p95 stress: 2,000/2,000 replicates retain positive sign
5%-bbox-shift-p95 stress: 2,000/2,000 replicates retain positive sign
```

These technical-stress runs test robustness to explicit symmetric random perturbation calibrated from the frozen audit summary. They do not establish gravity-referenced orientation accuracy, calibrated physical colour accuracy, or absence of environment-dependent measurement error.

### 7. Validate the reproduced v3 outputs

At minimum, the following machine-readable reports must exist and parse as JSON:

```text
work/results/biological_axes_report.json
work/results/sensitivity_chain/biological_axis_sensitivity_chain_report.json
work/integration/construct_scale_integration_report.json
work/upgrade/construct_scale_upgrade_report.json
work/upgrade/construct_scale_contrast_summary.json
work/results/assessability/assessability_selection_report.json
work/results/technical_stress/frozen_technical_error_stress_report.json
```

The corresponding frozen/current result interpretation is documented in [`analysis/v3/README.md`](analysis/v3/README.md) and its linked result receipts. If a source ZIP or extracted source file fails its recorded SHA-256 check, stop rather than substituting another file or silently recomputing from a different cohort.

## Third-party reproducibility status

The frozen Chapter 1 v2 analysis code, contracts, frozen outputs, validation reports and figure sources are public and auditable in this repository. The complete independent-reader procedure is [`reproducibility/README.md`](reproducibility/README.md).

The exact minimum numerical reproduction bundle is publicly archived on Zenodo v1 and has been independently re-downloaded without Zenodo author credentials:

```text
DOI       10.5281/zenodo.22295791
record    https://zenodo.org/records/22295791
bundle    azami_ch1_v2_reproduction_inputs_2026-09-04.zip
size      56,942,044 bytes
SHA-256   50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

The immutable code ref paired to that bundle is:

```text
584af97b050d15701f26ce1facea212d5b648d4d
```

On 2026-09-04, a clean GitHub-hosted Ubuntu runner downloaded the Zenodo bundle without author credentials or an access token. The downloaded bundle matched the frozen size and SHA-256, and all four embedded minimum analysis inputs matched the SHA-256 values in [`reproducibility/public_release_manifest.json`](reproducibility/public_release_manifest.json). The anonymous release-gate evidence is GitHub Actions run `33871722778`.

The public release gate tracked in issue #88 is therefore complete. The package status is **ready for independent third-party numerical reproduction**. This is a readiness statement about credential-free public availability and exact byte identity; an actual independent reproduction must still execute the frozen analysis and satisfy every validation criterion below.

## Quick independent reproduction

### 1. Check out the frozen code

```bash
git clone https://github.com/zuizui0223/azami.git
cd azami
git checkout 584af97b050d15701f26ce1facea212d5b648d4d
```

Use this immutable commit for the frozen Chapter 1 v2 reproduction rather than a later mutable `main` state.

### 2. Download and verify the public Zenodo bundle

Download DOI [`10.5281/zenodo.22295791`](https://doi.org/10.5281/zenodo.22295791) without author credentials. Before extraction, verify the bundle SHA-256:

```bash
# Linux
sha256sum azami_ch1_v2_reproduction_inputs_2026-09-04.zip

# macOS
shasum -a 256 azami_ch1_v2_reproduction_inputs_2026-09-04.zip
```

Expected SHA-256:

```text
50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

Extract the bundle into `reproduction_inputs/`. The four embedded archives must match:

| File | SHA-256 |
|---|---|
| `artifact-9612943217-continuous.zip` | `101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e` |
| `artifact-9632715852-multilevel.zip` | `51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6` |
| `artifact-8227254443-historical.zip` | `499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc` |
| `artifact-8983877726-spatial.zip` | `151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca` |

The repository verifier checks these identities fail-closed:

```bash
python reproducibility/verify_local_materials.py \
  --archive-dir reproduction_inputs \
  --scope numerical \
  --report local_material_check.json
```

Expected result: `"status": "PASS"`. If any hash differs, stop rather than substituting a file.

### 3. Run the canonical numerical rebuild

Create the frozen Python environment and install the pinned numerical dependencies:

```bash
python -m venv .venv
# activate .venv using the command appropriate for your platform
python -m pip install --upgrade pip
python -m pip install -r reproducibility/frozen-numerical-rebuild-requirements.txt
python -m pip install -e . --no-deps
```

Then follow [`reproducibility/README.md`](reproducibility/README.md) for source validation, nine-predictor environment reconstruction, sampling/spatial/historical analyses and figure rebuilding.

A successful independent numerical reproduction must reach all canonical frozen states:

```text
environment atlas:                PASS 24/24
sampling composition:             PASS 7/7
spatial + historical sensitivity: PASS 11/11
reconstructed environment SHA:    e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7
figure build:                     status ok; 5 main + 7 supporting figures
```

Finally run the repository integrity checks described in the runbook. The same invariants are guarded on `main` by [`.github/workflows/reproducibility-integrity.yml`](.github/workflows/reproducibility-integrity.yml).

Prepared public-deposition metadata is retained in [`reproducibility/zenodo_metadata.json`](reproducibility/zenodo_metadata.json), and the author-side publication/clean-session verification procedure is [`reproducibility/PUBLIC_DATA_RELEASE.md`](reproducibility/PUBLIC_DATA_RELEASE.md).

The five frozen CHELSA v2.1 process rasters may optionally be cached for a fully offline rerun. If they are not locally cached, their exact public URLs remain frozen in the source registry and the reconstructed nine-predictor environment must match the SHA-256 shown above.

## Owner-side recovery status

Separately from public distribution, the source artifacts are durably preserved in owner-controlled storage, and the complete pre-cleanup scientific tree remains permanently recoverable from immutable tag `azami-ch1-v2-2026-08-27` and branch `archive/ch1-precleanup-20260827`. These preservation locations are recovery infrastructure, not the public third-party distribution endpoint.

## Tracked scientific material

- `ch1_global/v2/` — image collection, screening, measurement, QC and historical-analysis implementation.
- `analysis/` — statistical analyses, sensitivity analyses and figure-generation code.
- `analysis/ch1/` — frozen analysis contracts and registries needed to define the estimands.
- `analysis/v3/` — biological-construct synthesis, robustness analyses and immutable run receipts layered on the frozen v2 measurements.
- `analysis_outputs/` — frozen machine-readable result tables and validation reports.
- `reproducibility/figures/` — generated result figures together with source/provenance material required to rebuild them.
- `reproducibility/contracts/` — compact trait definitions needed to reproduce measurement/analysis decisions.
- `reproducibility/public_release_manifest.json` — public third-party release gate, minimum input hashes, Zenodo identity and frozen code ref.
- `reproducibility/material_availability.json` — public/private/external material boundary.
- `reproducibility/recovery_inventory.json` — owner recovery completeness and separate public-release status.
- `reproducibility/durable_archive_manifest.json` — owner-side preservation identities and reassembly contract.
- `tests/` — regression and invariant tests.
- `src/azami_ch1/` — reusable utilities.

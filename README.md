# azami analysis and reproducibility repository

This repository retains the code and frozen machine-readable products needed to reproduce and audit the image-phenomics analyses.

Prepublication manuscript prose and journal-submission packages are intentionally kept outside GitHub.

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
- `analysis_outputs/` — frozen machine-readable result tables and validation reports.
- `reproducibility/figures/` — generated result figures together with source/provenance material required to rebuild them.
- `reproducibility/contracts/` — compact trait definitions needed to reproduce measurement/analysis decisions.
- `reproducibility/public_release_manifest.json` — public third-party release gate, minimum input hashes, Zenodo identity and frozen code ref.
- `reproducibility/material_availability.json` — public/private/external material boundary.
- `reproducibility/recovery_inventory.json` — owner recovery completeness and separate public-release status.
- `reproducibility/durable_archive_manifest.json` — owner-side preservation identities and reassembly contract.
- `tests/` — regression and invariant tests.
- `src/azami_ch1/` — reusable utilities.

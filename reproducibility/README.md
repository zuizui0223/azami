# Chapter 1 v2 independent reproducibility runbook

This is the public procedure for an **independent reader, reviewer, or researcher**. It must not depend on the author's PC, private cloud account, OneDrive path, or Google Drive credentials.

The machine-readable public-release contract is `reproducibility/public_release_manifest.json`.

## Current public status

The code, contracts, frozen results, validation reports and figure sources are public in this GitHub repository. A single public-release staging bundle has also been built and SHA-verified:

```text
azami_ch1_v2_reproduction_inputs_2026-09-04.zip
size:    56,942,044 bytes
SHA-256: 50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

The bundle contains the exact four minimum binary inputs for the numerical rerun. The immutable code ref paired to it is:

```text
584af97b050d15701f26ce1facea212d5b648d4d
```

**The bundle is prepared but is not yet the public reproducibility endpoint.** Until `public_release_manifest.json` contains a credential-free DOI or stable URL that has been independently re-downloaded and verified, full third-party numerical reproducibility remains blocked. GitHub issue #88 tracks that final publication gate.

## 1. Obtain the code

```bash
git clone https://github.com/zuizui0223/azami.git
cd azami
git checkout 584af97b050d15701f26ce1facea212d5b648d4d
```

Do not substitute a later mutable `main` state when reproducing the frozen package.

## 2. Obtain the public data bundle

After issue #88 is closed, obtain the bundle from the DOI/stable URL recorded in:

```text
reproducibility/public_release_manifest.json
```

The download must not require author credentials. Verify the **bundle itself** before extraction:

```text
SHA-256 50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677
```

Then extract it into any local directory, for example `reproduction_inputs/`.

## 3. Verify the embedded analysis inputs

The extracted directory must contain these four exact files:

| File | SHA-256 |
|---|---|
| `artifact-9612943217-continuous.zip` | `101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e` |
| `artifact-9632715852-multilevel.zip` | `51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6` |
| `artifact-8227254443-historical.zip` | `499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc` |
| `artifact-8983877726-spatial.zip` | `151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca` |

Run the repository verifier:

```bash
python reproducibility/verify_local_materials.py \
  --archive-dir reproduction_inputs \
  --scope numerical \
  --report local_material_check.json
```

Expected result:

```text
"status": "PASS"
```

If any input hash differs, stop. Do not continue with a substituted file.

## 4. Recreate the frozen Python environment

```bash
python -m venv .venv
```

Activate the environment using the platform-appropriate command, then:

```bash
python -m pip install --upgrade pip
python -m pip install -r reproducibility/frozen-numerical-rebuild-requirements.txt
python -m pip install -e . --no-deps
```

The pinned requirements file is the preferred numerical reproduction environment.

## 5. Validate the source checkpoints before analysis

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip reproduction_inputs/artifact-9612943217-continuous.zip \
  --multilevel-zip reproduction_inputs/artifact-9632715852-multilevel.zip \
  --out-dir rebuild_check \
  --mode validate_sources
```

This is fail-closed. Among other identities it verifies:

- continuous source ZIP SHA-256;
- multilevel source ZIP SHA-256;
- continuous trait table SHA-256 `d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344`;
- strict-spatial core environment SHA-256 `2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0`;
- complete-18 process reference SHA-256 `1ab84254a80493776b4c435152ed3d2a1c1e68dd0e0342da0ea081eeb5cd3d9b`;
- 46,276 strict-spatial observations / 259 taxa;
- 1,874 complete-18 process-reference rows / 124 taxa.

## 6. Reconstruct the nine-predictor environment and full-27 atlas

Five process-extension rasters are frozen in:

```text
analysis/ch1/chelsa_process_environment_sources.json
```

They are CHELSA v2.1 layers for rsds, VPD, surface wind, GSP and NPP. If the public data release includes exact raster copies, place them in a local directory and pass `--raster-cache-dir`. Otherwise the canonical runner retrieves the exact frozen public CHELSA URLs.

Run:

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip reproduction_inputs/artifact-9612943217-continuous.zip \
  --multilevel-zip reproduction_inputs/artifact-9632715852-multilevel.zip \
  --out-dir rebuild_full27 \
  --mode rebuild_full27
```

If using a public raster cache, append:

```text
--raster-cache-dir <CHELSA_CACHE_DIR>
```

The reconstructed 46,276-row nine-predictor environment must match:

```text
e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7
```

Expected full-27 validation:

```text
PASS 24/24
```

## 7. Recover the sampling-composition auxiliary inputs

### Broad-region lookup

Extract the exact frozen table from the public spatial input ZIP:

```python
import zipfile
from pathlib import Path

src = Path("reproduction_inputs/artifact-8983877726-spatial.zip")
out = Path("rebuild_inputs")
out.mkdir(exist_ok=True)
with zipfile.ZipFile(src) as z:
    (out / "broad_region_lookup.csv").write_bytes(
        z.read("spatial_regions/broad_region_lookup.csv")
    )
```

Expected SHA-256:

```text
085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff
```

### Native-status table

The exact frozen table is public in immutable Git history. Preserve raw bytes:

```python
import subprocess
from pathlib import Path

out = Path("rebuild_inputs")
out.mkdir(exist_ok=True)
(out / "observation_native_status.csv").write_bytes(
    subprocess.check_output([
        "git", "show",
        "azami-ch1-v2-2026-08-27:analysis_outputs/native_range_sensitivity_v1/observation_native_status.csv",
    ])
)
```

Expected SHA-256:

```text
c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a
```

## 8. Rebuild sampling, spatial, historical, and secondary analyses

```bash
python analysis/rebuild_frozen_analysis.py \
  --continuous-zip reproduction_inputs/artifact-9612943217-continuous.zip \
  --multilevel-zip reproduction_inputs/artifact-9632715852-multilevel.zip \
  --historical-zip reproduction_inputs/artifact-8227254443-historical.zip \
  --out-dir rebuild_everything \
  --mode rebuild_all \
  --with-sampling \
  --regions rebuild_inputs/broad_region_lookup.csv \
  --native-status rebuild_inputs/observation_native_status.csv \
  --with-spatial \
  --with-historical
```

Append `--raster-cache-dir <CHELSA_CACHE_DIR>` if using a public/local copy of the five frozen rasters.

Expected validation states:

```text
environment atlas:                PASS 24/24
sampling composition:             PASS 7/7
spatial + historical sensitivity: PASS 11/11
```

The historical source contains S1, S3 and 50 randomized S2 placement trees: 52 placement trees in total.

## 9. Rebuild the figures

```bash
python analysis/figures/run_image_to_trait_perturbation_audit.py
python analysis/figures/build_v2_figures.py
```

Expected state:

```text
reproducibility/figures/figure_build_report.json
status: ok
main figures: 5
supporting figures: 7
```

The required figure-source tables are tracked in GitHub under `reproducibility/figures/source/`.

## 10. Final repository integrity check

```bash
python -m compileall -q analysis src ch1_global/v2 reproducibility/verify_local_materials.py
python -m pytest -q \
  tests/test_ch1_canonical_manifest.py \
  tests/test_reproducibility_manifest.py \
  tests/test_material_availability.py \
  tests/test_azami_ch1_provenance.py \
  tests/test_azami_ch1_tabular.py
```

The same invariants are guarded by `.github/workflows/reproducibility-integrity.yml`.

## Reproduction success criterion

A third party has independently reproduced the frozen Chapter 1 v2 numerical package only when all of the following are true:

1. code was checked out at `584af97b050d15701f26ce1facea212d5b648d4d`;
2. the public bundle was downloaded without author credentials and matched SHA-256 `50ec15b1280d4660839ca4bf0d55c970a5f49b4d4feaabb7a073b73500253677`;
3. the four embedded analysis inputs passed their recorded SHA-256 checks;
4. the reconstructed full environment matched `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7`;
5. validations returned 24/24, 7/7 and 11/11 PASS;
6. figure build returned `status: ok`, 5 main and 7 supporting figures.

Until a DOI/stable public URL is filled in `reproducibility/public_release_manifest.json`, step 2 cannot be completed by an independent researcher, so the repository correctly does **not** claim full third-party numerical reproducibility yet.

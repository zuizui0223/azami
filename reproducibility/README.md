# Chapter 1 v2 reproducibility runbook

This directory contains the frozen code, machine-readable results, provenance, and archive contracts needed to audit and reproduce the Chapter 1 v2 image-phenomics analysis.

The material inventory was re-audited against GitHub on 2026-09-04. The canonical machine-readable audit is `reproducibility/material_availability.json`.

## 1. What is actually available

Do not assume that every large binary is stored inside the Git repository.

| Material | Where it is | Status |
|---|---|---|
| analysis code, contracts, predictor/trait registries | GitHub `main` | present |
| frozen Python numerical environment | `reproducibility/frozen-numerical-rebuild-requirements.txt` | present |
| frozen full-27 result tables and validation reports | `analysis_outputs/v2_full27_environment_atlas_2026-08-27/` | present |
| figure builders, figure source tables, frozen figures | `analysis/figures/` and `reproducibility/figures/` | present |
| exact continuous numerical source ZIP, artifact `9612943217` | owner-controlled Google Drive archive | durably archived and SHA-verified |
| exact multilevel/process source ZIP, artifact `9632715852` | owner-controlled Google Drive archive | durably archived and SHA-verified |
| historical tree ZIP, artifact `8227254443` | owner-controlled Google Drive archive | durably archived and SHA-verified |
| spatial/broad-region ZIP, artifact `8983877726` | owner-controlled Google Drive archive | durably archived and SHA-verified |
| exact frozen native-status CSV | immutable Git tag `azami-ch1-v2-2026-08-27` | retained in Git history and independently rebuildable |
| detector / Figure 1 / independent audit / exhaustive upstream binaries | owner-controlled Google Drive archive | durably archived; not required for the frozen full-27 numerical rerun |
| five added CHELSA process rasters: rsds, VPD, wind, GSP, NPP | frozen public CHELSA v2.1 URLs | **not vendored in Git or the owner archive as of this audit** |

The durable Drive archive identity, Drive file IDs, SHA-256 values, chunk hashes, and reassembly instructions are committed in `reproducibility/durable_archive_manifest.json`. The broader Actions provenance is in `reproducibility/actions_artifact_catalog.json`.

### Important boundary: the five CHELSA process rasters

The two archived numerical ZIPs contain the frozen continuous-trait universe, the 46,276-row strict-spatial core environment, and the complete-18 process-environment reference. The canonical full-27 rebuild then samples five additional frozen CHELSA v2.1 rasters identified in `analysis/ch1/chelsa_process_environment_sources.json`.

Therefore:

- auditing the frozen results and validating the direct archived numerical inputs does **not** depend on GitHub Actions retention;
- rebuilding the full nine-predictor environment from the source rasters currently needs either network access to those five frozen CHELSA URLs or a local cache containing exact URL-basename copies;
- the reconstructed 46,276-row nine-predictor table must have SHA-256 `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7`; otherwise the canonical runner stops before analysis.

No replacement environment variable is allowed if a frozen source is unavailable.

## 2. Recommended Windows / OneDrive layout

The recommended local root is:

```text
C:\Users\zuizui\OneDrive - Kyoto University\デスクトップ\azami論文材料
```

A convenient layout is:

```text
azami論文材料\
├─ azami\
│  └─ ... GitHub repository ...
└─ azami_ch1_reproducibility_archive_2026-09-04\
   ├─ artifact-9612943217-continuous.zip
   ├─ artifact-9632715852-multilevel.zip
   ├─ artifact-8227254443-historical.zip
   ├─ artifact-8983877726-spatial.zip
   ├─ artifact-8076736948-detector.zip
   ├─ artifact-8066010557-figure1-metadata.zip
   ├─ artifact-8099953404-figure1-global.zip
   ├─ artifact-8225059018-figure1-continuous.zip
   ├─ artifact-8521925441-qc.zip
   ├─ artifact-8521926057-predictions.zip
   ├─ artifact-8521924881-annotation.part-00
   ├─ artifact-8521924881-annotation.part-01
   ├─ exhaustive45-part-00.zip
   ├─ ...
   ├─ exhaustive45-part-19.zip
   └─ exhaustive45-reassembly-manifest.zip
```

The directory name may differ. The commands below use `$ARCHIVE`, so set that variable to the actual synchronized archive directory.

## 3. Clone and pin the recovered repository state

Open PowerShell:

```powershell
$MATERIALS = "C:\Users\zuizui\OneDrive - Kyoto University\デスクトップ\azami論文材料"
Set-Location $MATERIALS

git clone https://github.com/zuizui0223/azami.git
Set-Location (Join-Path $MATERIALS "azami")
git checkout main
git pull
```

For an exactly pinned checkpoint, use the recovery-complete commit recorded in `reproducibility/recovery_inventory.json` or the relevant later `main` commit. The immutable pre-cleanup scientific tree remains recoverable from:

```text
azami-ch1-v2-2026-08-27
archive/ch1-precleanup-20260827
```

## 4. Point to the synchronized durable archive

For the recommended local layout:

```powershell
$ARCHIVE = Join-Path $MATERIALS "azami_ch1_reproducibility_archive_2026-09-04"
```

Confirm that PowerShell sees it:

```powershell
Test-Path $ARCHIVE
Get-ChildItem $ARCHIVE | Select-Object Name, Length
```

## 5. Verify the local archive bytes before analysis

The repository includes a local verifier that uses the committed SHA ledger rather than trusting filenames.

### Numerical inputs only

```powershell
python reproducibility\verify_local_materials.py `
  --archive-dir "$ARCHIVE" `
  --scope numerical `
  --report local_material_check.json
```

Expected result:

```text
"status": "PASS"
```

This verifies the exact bytes of the continuous source, multilevel source, historical source, and spatial/broad-region source.

### Entire durable archive

To additionally verify the detector/Figure 1 sources, independent-audit packet, annotation reassembly, and all 20 exhaustive chunks:

```powershell
python reproducibility\verify_local_materials.py `
  --archive-dir "$ARCHIVE" `
  --scope all `
  --report local_material_check_all.json
```

For the 923 MB exhaustive checkpoint the verifier opens each Drive ZIP, hashes its raw chunk, streams the chunks in numeric order, and requires the original source SHA-256:

```text
5f18b42d18cfcb81691c38ce0f04bcef754e6a67382025ea90110dbc50ae194b
```

These upstream audit files are preservation material; they are not required to rerun the frozen full-27 numerical atlas.

## 6. Recreate the frozen Python environment

From the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r reproducibility\frozen-numerical-rebuild-requirements.txt
python -m pip install -e . --no-deps
```

The looser developer environment is still available as:

```powershell
python -m pip install -e ".[full,test]"
```

but the pinned requirements file is the preferred route for byte-identical numerical reconstruction.

## 7. Validate the two canonical numerical source ZIPs

Run this before any environmental resampling or statistical analysis:

```powershell
python analysis\rebuild_frozen_analysis.py `
  --continuous-zip "$ARCHIVE\artifact-9612943217-continuous.zip" `
  --multilevel-zip "$ARCHIVE\artifact-9632715852-multilevel.zip" `
  --out-dir rebuild_check `
  --mode validate_sources
```

The runner fails closed unless all of the following identities match the frozen record:

- continuous ZIP SHA-256 `101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`;
- multilevel ZIP SHA-256 `51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6`;
- continuous trait table SHA-256 `d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344`;
- strict-spatial core environment SHA-256 `2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0`;
- complete-18 process reference SHA-256 `1ab84254a80493776b4c435152ed3d2a1c1e68dd0e0342da0ea081eeb5cd3d9b`;
- 46,276 observation IDs / 259 taxa in the strict-spatial cohort;
- 1,874 rows / 124 taxa in the complete-18 process reference.

If this step fails, do not continue to the statistical analysis.

## 8. Optional: make the CHELSA process-raster stage offline

The frozen raster URLs are in:

```text
analysis/ch1/chelsa_process_environment_sources.json
```

The five expected URL-basename filenames are:

```text
CHELSA_rsds_1981-2010_mean_V.2.1.tif
CHELSA_vpd_mean_1981-2010_V.2.1.tif
CHELSA_sfcWind_mean_1981-2010_V.2.1.tif
CHELSA_gsp_1981-2010_V.2.1.tif
CHELSA_npp_1981-2010_V.2.1.tif
```

If you have downloaded those exact files into a directory, first check that all five expected names are present:

```powershell
$RASTER_CACHE = Join-Path $MATERIALS "CHELSA_v2.1_process_rasters"
python reproducibility\verify_local_materials.py `
  --archive-dir "$ARCHIVE" `
  --scope numerical `
  --raster-cache-dir "$RASTER_CACHE" `
  --report local_material_check_with_rasters.json
```

Then append this to the rebuild commands below:

```text
--raster-cache-dir "$RASTER_CACHE"
```

If `--raster-cache-dir` is omitted, the canonical sampler uses the frozen public CHELSA URLs. In either case, the resulting nine-predictor table must match the frozen full-environment SHA before the analysis runs.

## 9. Rebuild the primary full-27 atlas

```powershell
python analysis\rebuild_frozen_analysis.py `
  --continuous-zip "$ARCHIVE\artifact-9612943217-continuous.zip" `
  --multilevel-zip "$ARCHIVE\artifact-9632715852-multilevel.zip" `
  --out-dir rebuild_full27 `
  --mode rebuild_full27
```

If using a local raster cache, add:

```powershell
  --raster-cache-dir "$RASTER_CACHE"
```

The reconstructed full environment must match:

```text
e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7
```

The full-27 validation should finish as:

```text
PASS 24/24
```

The canonical frozen reference is:

```text
analysis_outputs/v2_full27_environment_atlas_2026-08-27/v2_full27_environment_validation.json
```

## 10. Recover the two exact sampling-composition auxiliary inputs

### A. Broad-region lookup

The exact frozen lookup is inside the archived spatial ZIP. Extract it without changing bytes:

```powershell
New-Item -ItemType Directory -Force rebuild_inputs | Out-Null
$env:AZAMI_ARCHIVE = $ARCHIVE
python -c "import os,zipfile,pathlib; a=pathlib.Path(os.environ['AZAMI_ARCHIVE']); o=pathlib.Path('rebuild_inputs'); z=zipfile.ZipFile(a/'artifact-8983877726-spatial.zip'); o.joinpath('broad_region_lookup.csv').write_bytes(z.read('spatial_regions/broad_region_lookup.csv'))"
```

Expected SHA-256:

```text
085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff
```

A fail-closed rebuild from Natural Earth 1:50m is also retained in `analysis/rebuild_frozen_broad_region_lookup.py`.

### B. Native-status table

The exact frozen table remains in the immutable Git tag. Use Python to preserve the raw Git blob bytes rather than PowerShell text redirection:

```powershell
python -c "import subprocess,pathlib; o=pathlib.Path('rebuild_inputs'); o.mkdir(exist_ok=True); o.joinpath('observation_native_status.csv').write_bytes(subprocess.check_output(['git','show','azami-ch1-v2-2026-08-27:analysis_outputs/native_range_sensitivity_v1/observation_native_status.csv']))"
```

Expected SHA-256:

```text
c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a
```

The table is also independently rebuildable through `analysis/rebuild_frozen_native_status.py`, which pins the WCVP dataset identity and TDWG source commit.

## 11. Rebuild sampling, spatial, historical, and secondary analyses

Once `rebuild_inputs/broad_region_lookup.csv` and `rebuild_inputs/observation_native_status.csv` exist:

```powershell
python analysis\rebuild_frozen_analysis.py `
  --continuous-zip "$ARCHIVE\artifact-9612943217-continuous.zip" `
  --multilevel-zip "$ARCHIVE\artifact-9632715852-multilevel.zip" `
  --historical-zip "$ARCHIVE\artifact-8227254443-historical.zip" `
  --out-dir rebuild_everything `
  --mode rebuild_all `
  --with-sampling `
  --regions "rebuild_inputs\broad_region_lookup.csv" `
  --native-status "rebuild_inputs\observation_native_status.csv" `
  --with-spatial `
  --with-historical
```

Again, append `--raster-cache-dir "$RASTER_CACHE"` if using a local CHELSA cache.

Expected frozen validation states are:

```text
environment atlas:              PASS 24/24
sampling composition:           PASS 7/7
spatial + historical sensitivity: PASS 11/11
```

Historical input identity is also fail-closed. Artifact `8227254443` contains deterministic S1/S3 trees and 50 randomized S2 trees, for 52 historical-placement model fits.

## 12. Rebuild the frozen figures

Canonical entry points:

```powershell
python analysis\figures\run_image_to_trait_perturbation_audit.py
python analysis\figures\build_v2_figures.py
```

Expected `reproducibility/figures/figure_build_report.json` state:

```text
status: ok
main figures: 5
supporting figures: 7
```

The figure source tables required by these builders are tracked under `reproducibility/figures/source/`.

## 13. Final repository integrity check

Run the same core checks used by CI:

```powershell
python -m compileall -q analysis src ch1_global\v2 reproducibility\verify_local_materials.py
python -m pytest -q `
  tests\test_ch1_canonical_manifest.py `
  tests\test_reproducibility_manifest.py `
  tests\test_azami_ch1_provenance.py `
  tests\test_azami_ch1_tabular.py
```

`main` is guarded by `.github/workflows/reproducibility-integrity.yml`, which checks the live manifest, frozen validation reports, durable archive metadata, and the code-only repository boundary.

## 14. Which files are canonical?

Use these as the first places to look:

- stage inventory: `ch1_global/v2/ANALYSIS_MANIFEST.tsv`;
- audited material inventory: `reproducibility/material_availability.json`;
- durable off-Actions archive: `reproducibility/durable_archive_manifest.json`;
- Actions provenance: `reproducibility/actions_artifact_catalog.json`;
- overall recovery status: `reproducibility/recovery_inventory.json`;
- numerical rebuild runner: `analysis/rebuild_frozen_analysis.py`;
- frozen numerical environment: `reproducibility/frozen-numerical-rebuild-requirements.txt`;
- active frozen results: `analysis_outputs/v2_full27_environment_atlas_2026-08-27/`;
- frozen figure products and source tables: `reproducibility/figures/`.

## 15. Historical recovery boundary

The active `main` tree intentionally does not restore superseded manuscript/reviewer/submission-era material. The pre-cleanup scientific tree is preserved at:

```text
azami-ch1-v2-2026-08-27
archive/ch1-precleanup-20260827
```

These are recovery points, not the active analysis lane.

Manuscript prose, reviewer material, and journal submission packages are not part of the active reproducibility tree.

# Current Chapter 1 numerical reproduction

## One-command current replay

From a clean checkout of the cleaned `main`, record `git rev-parse HEAD`, create a Python 3.12 environment, then run:

```bash
python -m pip install -r reproducibility/requirements-current.txt
python -m reproducibility.run_current_analysis --download
```

`--download` fetches only missing numerical archives through authenticated GitHub CLI (`gh`). It never fetches photographs. Without it, supply the four exact archives in `work/current/archives/`, named `artifact-ID-ROLE.zip` for roles `continuous`, `environment`, `spatial`, `historical`. The public v2 ZIP contains three of these four archives; the environment input requires the current Actions artifact until the new Zenodo version is published. Native status is recovered from the immutable tag, with only the frozen permitted newline normalization.

The runner verifies ZIPs and extracted files, runs all seven existing numerical steps with unchanged seeds/counts, and compares 15 aggregate outputs with the byte-verified original references. Results and logs go under ignored `work/current/`. A successful run ends in `work/current/results/validation.json` with `status: PASS`. A missing input or discrepancy stops the run. `--prepare-only` verifies inputs without claiming numerical completion.

Two fields in the older archived JSON reports were already removed from scientific main `fe25abd`: the pairwise report's `environment_signature_alignment` and the common-cohort report's `integration_environment_coupling`. The validator explicitly excludes those retired fields and reports their names; all current fields and archived CSV values are compared. The original artifact bytes are preserved, and these retired analyses are not restored or presented as current results.

The upgrade report's known historical `claim_boundary` wording is mapped to its already-edited `fe25abd` wording for comparison. No other prose or numerical value is normalized. This is numerical equivalence for the current scope, not byte-for-byte identity of regenerated reports with older artifacts.

The detailed commands below describe the original merged scientific revision. They remain useful for inspecting individual stages; use the single current runner above for the cleaned layout. Frozen endpoint implementations have moved to `legacy/` and are not alternative current entry points.

## Original merged scientific revision and detailed stages


Chapter 1 v3 is a higher-level synthesis layered on the frozen v2 measurements. It does **not** retrain the detector, redownload photographs, or rerun the retired categorical classifier. The v3 workflows begin from frozen v2 trait/environment artifacts and reproduce the construct-level analyses, scale comparisons and technical/selection audits.

The v3 synthesis was merged from PR #93 into `main`. For an exact reproduction of the merged v3 analysis surface, use merge commit `477fd73fc8ebcbad5603775cc0dc980c418c6268` together with the byte-verified source artifacts listed here.

### 1. Check out the merged v3 analysis surface

```bash
git clone https://github.com/zuizui0223/azami.git
cd azami
git checkout 477fd73fc8ebcbad5603775cc0dc980c418c6268
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

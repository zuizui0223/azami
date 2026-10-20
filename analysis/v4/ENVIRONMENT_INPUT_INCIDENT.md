# GSP input incident — interpretation hold

Detected on 2026-10-20 while visually inspecting the observation-level gradient
plots from Actions run 39420523666.

Environment-only audit run 39421038043 confirms that 1,239 of 313,350 native
observations have `chelsa_gsp = 4294967295`, the uint32 maximum. Its median is
6690 and 99th percentile 38206.45; the extreme value is the most frequent single
value. These values entered species means, scaling, VIF selection and fitted
models as finite observations. This is not simply a plotting-unit issue.

The original raster metadata is being inspected in run 39421269354 against the
checksum-locked GSP cache artifact 12243306531. The exact upstream nodata
representation and physical-unit scaling must be resolved before a correction
is finalized. Do not infer the missingness mechanism solely from the value.

## Affected evidence

Any analysis using the uncorrected GSP column is provisional, including selected
environment decomposed run 39416125522, common-slope run 39418226382, decomposed
replay 39420155614 and expanded head-outline run 39420523666. Their inputs and
outputs remain archived for provenance, but their adjusted coefficients and
P values must not be presented as validated ecological conclusions.

## Required recovery

1. Verify source metadata, dtype, declared nodata and scale/offset against the
   encoded values. Audit the other eight exposures as well.
2. Apply a provenance-labelled environmental decoding correction without
   altering raw source archives or trait measurements. Missing values stay
   missing; never replace them with zero or neighbouring pixels.
3. Add regression tests for the observed failure and guard against its recurrence.
4. Repeat environment-only selection under the same declared rule, then refit
   all four modules with endpoint-specific eligibility and report changed counts.
5. Reinspect gradients and update conclusions from corrected results, regardless
   of whether previous associations persist.

## Verified decoding evidence

Actions 39421799730 verified all nine cache archive checksums and exported original
raster metadata, now stored in environment_source_metadata.json. GSP source hash
342858ce86f8b5c4f088327e4fa68afc7ac9cb64eff46ffaa980577a0f0641b6 declares
nodata=2147483647, uint32, scale=0.1 and offset=0; uint32-max is an additional
invalid encoding, not its declared nodata. native_665115_pipeline.sample stored
sample_tiles output directly without scale/offset, explaining the raw-unit archive.

Recovery is implemented as a separate decoder of the exact immutable raw table.
It masks both the declared GSP nodata and the additional invalid uint32-max,
records their treatment, applies each verified affine transform once, and retains
all trait values, eligibility and native observations. Analysis-specific environment
missingness is applied later. Corrected numerical validation remains pending.

Correction validated in Actions 39422378058 at commit 93e8d87: 23 tests passed, 312001 complete nine-environment observations, and all four modules fitted. Corrected aggregate results are archived under results/20261020_corrected_environment/. The earlier between-species rainfall–orientation association is not supported after correction. Follow-up robustness checks remain outstanding.

# azami v4: native-first ORIGINAL upper-source audit (2026-10-20)

**Exploratory. Frozen v3, GEB paper, and main have not been changed.**

## Source population — correction to the earlier 460,036-observation scope

The original pre-detector 2026-07-03 metadata archive has **665,115 observations and 1,122,854 photo links**. The previous 460,036 / 777,766 historical processing queue was a strict subset. The separately repaired inventory (665,139 observations) has 24 observations not present in this checksum-pinned original snapshot.

WCVP June 2026 authority was queried for **all 350 exact species-rank names**. Another 117 genus/infraspecific/hybrid labels are explicitly unresolved rather than being truncated to species. Original locations were tested against pinned TDWG level-3 polygons; records without a clear native distribution row, or with uncertain identification, geoprivacy, coordinates or cultivated status, are **not** called native.

| Original-source result | Number |
|---|---:|
| Source observations | 665,115 |
| Source photo links | 1,122,854 |
| Native-eligible observations | **313,350** |
| Native-eligible photo links | **534,186** |
| Accepted WCVP taxa in native set | **306** |
| Increase from previous 288-name authority | +114 native observations |
| Native observations outside historical 460,036-observation processing queue | **68,780** |

Native describes region-level WCVP/TDWG agreement, **not the indigenous origin of a specific photographed individual**. Native eligibility never uses flower colour, orientation, detection success or ecological outcomes.

Exact original metadata ZIP SHA256: 2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae.
Extended 350-name WCVP JSON SHA256: a5bf23efffbf8e8f9459d5c3391a0136adcedd8f9dbff0117caa969ab545c2bd.

## Full-source original-resolution CHELSA sampling

The nine existing SHA-authenticated 1981–2010 CHELSA v2.1 raster caches were sampled for the new native roster, *without re-downloading rasters*. Their 5-degree tile identifiers are storage groups, not environmental means.

| Coverage | Observations |
|---|---:|
| Full native roster sampled | 313,350 |
| All nine original-pixel predictors finite | **313,240 (99.965%)** |
| Outside cached original-pixel tiles | 106 across 13 extra cells |
| Additional original wind raster nodata | 4 |
| Any environmental value unavailable (never imputed) | **110** |

On all old native-overlap records, exact values match the previous independent extraction: 27,066/27,066 per variable, except wind 27,064/27,064. The new observation-level 313,350-row table with IDs and coordinates remains private and is **not** added to public GitHub.

Local entry points:
- python -m analysis.v4.native_665115_pipeline classify --metadata-zip SOURCE.zip --authority-json native_sources_350.json --private-dir PRIVATE
- python -m analysis.v4.native_665115_pipeline sample --native-roster PRIVATE/native350_roster.csv.gz --cache-root CHELSA_CACHE_DIR --private-dir PRIVATE

These scripts verify source checksums, all 9 tile archive members and expected missingness. They do not impute raster nodata.

## Native-wide phenotype results and honest boundary

Archived image measurements are available for **217,104** of the 313,350 native observations. Another **96,246 native observations lack archived measurements**. The new 68,780 observations outside the historical processing queue have no new measured phenotype in this replay.

Reusing the existing, explicitly uncalibrated photographic measurements gives the same real fitted ecological sample as the earlier PR: presentation n=169,950 / 178 accepted taxa and visible colour n=186,264 / 177 taxa. Among-taxon geographically blocked cross-validation R² is **+0.01982** for image-referenced presentation, **+0.13810** for joint visible colour (L*, C*, hue sin/cos); within-taxon CV R² is **−0.00243** and **−0.00178**. This is **not** independent replication from newly acquired images.

Same-record presentation–colour integration: n=164,052 / 177 accepted taxa; raw RV is 0.00624 within vs 0.03679 among. The taxon bootstrap plus matched-n median difference is +0.02515, 95% interval −0.01982 to +0.09158 (1000 draws). That interval crosses zero: the full native data **do not establish a universal stronger among-taxon linkage**. A circular joint-hue sensitivity also includes zero. Previously seen marginal chroma–solar and presentation–rain associations remain descriptive and cannot be called uniquely causal after joint space and climate adjustment.

## Real-image queue

The 96,246 missing-measurement observations contain **148,939 original photo candidates**, from **277 accepted taxa and 1,360 taxon-by-5-degree cells**. Of these, 68,780 were never in the old processing queue and 27,466 were queued without archived measurement. The private queue retains all photo candidates and uses taxon/cell-balanced scheduling, not outcome-driven selection. A conservative CC0 / CC BY / CC BY-SA public-diagnostic first-photo option exists for 11,162 observations. Another 85,084 first choices need rights review; this **does not** establish that every noncommercial research use is prohibited.

## Verification and unresolved GitHub input

Local reproduction on the exact original archived metadata and new authority completed: **665,115 classified; 313,350 native; 313,240 complete CHELSA records**; 3 additional synthetic tests passed. The previously measured sample gave the same spatially held-out model R² as before.

The original artifact 8066010557 expired (HTTP 410). Its identical bytes have now been restored through an encrypted asset in the existing private draft Release, with the key held separately in Actions secrets. [Run 39407663676](https://github.com/zuizui0223/azami/actions/runs/39407663676) completed the full-source native and climate replay on GitHub: 665,115 original observations, 313,350 native-eligible observations and 313,240 complete nine-predictor records. Aggregate source and climate receipts were downloaded and inspected. No plaintext observation rows were exported. This verifies classification and environmental extraction, not new image measurement or ecological models.

No raw observation IDs, original image URLs, exact coordinates or rights ledgers may be attached to public PR files or workflow artifacts. Do not mistake 99.965% climate coverage for 99.965% actual flower-trait measurement.

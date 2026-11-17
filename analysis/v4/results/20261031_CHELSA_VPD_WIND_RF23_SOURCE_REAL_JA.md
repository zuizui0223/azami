# 原CHELSA VPD・地表風を加えた21→23変数RF：8形質・地理外挿の実測結果（2026-10-31）

**STATUS:** authorized original S1 source-real exploratory. [Successful run #41187988423](https://github.com/zuizui0223/azami/actions/runs/41187988423), artifact **13037393509**. No edits to 2026-10-24 frozen S1/S2 Holm estimates, GEB manuscript, H3c or source images. This is not a confirmatory P-value family.

## Input and registered scientific difference

- Previous 12→21 result: see [full 8-trait and posthoc sensitivity](20261031_RF_BIO21_EIGHT_TRAITS_AND_POSTHOC_SENSITIVITY_JA.md).
- Starting from the **same 313,350 original S1 observations**, exact SHA-verified existing 7+9 original climate raster extractions produced BIOCLIM21. Two extra fields were **already in the original decrypted S1** (not newly downloaded TerraClimate data): `chelsa_vpd_mean` and `chelsa_sfcwind_mean` with source hashes `35328724237339b2f8eb80f21adb7ee38f58aebd17bd18d1c271d563a18fe8e8` and `d4b88b69c98ef50522163727766137a91bc2ffde1b0cedcf208b5230d1835f08`.
- Their original 1981–2010 CHELSA decoding receipt and scale/offset were checked against committed raster metadata. Finite VPD **313,244/313,350**; wind **313,240/313,350**. This is archive source verification, **not independent local measurements of wind exposure or atmospheric stress at the flower**.
- **21 versus 23 RF use exactly the same 23-complete rows/cells per trait**, not a comparison of differing cohorts or the earlier 21-versus-12 table. 80 trees, same fixed hyperparameters and equal-taxon weighting, no raw geography/taxon identifiers as RF predictors. A positive `paired_skill_gain_upper_minus_lower` merely indicates lower prediction loss than 21 on the paired holdout.

## Strict original 5-degree / ≥3 photographs per cell cohort

Across these three domains, original photographic QC and ≥40 photos, ≥8 5°-cells per taxon retained 26 taxa (same group count; 647 angle, 659 colour, 653 outline cells).

| Construct | Unseen taxa 21 | Unseen taxa 23 | Δ23−21 | Held-out sectors 21 | Held-out sectors 23 | Δ23−21 |
|---|---:|---:|---:|---:|---:|---:|
| Image-frame angle | −0.1204 | −0.1343 | −0.0139 | −0.0388 | −0.0008 | +0.0380 |
| Lightness | −0.2055 | −0.2123 | −0.0067 | −0.0156 | −0.0602 | −0.0446 |
| Chroma | +0.0395 | +0.0304 | −0.0091 | −0.2050 | −0.1788 | +0.0262 |
| Joint circular hue | −0.3532 | −0.3717 | −0.0185 | −0.3410 | −0.3375 | +0.0034 |
| Aspect | −0.2771 | −0.2994 | −0.0223 | −0.0597 | −0.1002 | −0.0405 |
| Circularity | −0.1833 | −0.1752 | +0.0081 | −0.1825 | −0.1873 | −0.0048 |
| Solidity | −0.2074 | −0.1813 | +0.0261 | −0.1648 | −0.2261 | −0.0613 |
| Width profile CV | −0.1707 | −0.1576 | +0.0131 | −0.0822 | −0.1013 | −0.0191 |

Positive Δ on **3/8 unknown-taxon** and **3/8 unknown-sector** endpoints; negative Δ on the other 5/8 in each target. **Absolute 23-variable skill is negative for 7/8 unknown-taxon and 8/8 unknown-sector endpoints.** Near-zero −0.0008 for orientation region ≠ demonstrated positive predictability or actual rain shielding. Chroma is only weakly positive in unseen taxa (+0.0304).

Domain-group holdout permutation importance for VPD/wind is mostly small or negative, and correlated climate fields cannot be given independent causal attribution. The two included fields do not measure instantaneous event weather, floral temperature or physical mechanical response.

## Distinct posthoc support and fold design
Separate sensitivity after inspection of the original negative RF used minimum **one photo** per 5° cell, with taxon-count-balanced five-fold group CV (6–7 taxa per fold). Keep the original 3-photo/cell, 26-taxon estimates unchanged.

| Construct | Support (taxa/cells) | Balanced unknown taxa 21→23 | Unknown sector 21→23 |
|---|---|---|---|
| Angle | 31 / 818 | −0.1444 → −0.1636 | −0.0499 → −0.0124 |
| Chroma | 34 / 849 | −0.0106 → +0.0192 | +0.0587 → +0.0318 |
| Aspect | 32 / 827 | −0.2358 → −0.2544 | −0.0693 → −0.0707 |

Chroma changes sign between split/gate conditions, and near-zero positive skills are not robust demonstrations. Because this is a **posthoc** cohort change, no confirmatory label or new multiplicity correction can make the effect independent evidence.

## Decision before pursuing more predictors

The first 21-variable and next 23-variable original-source RF experiments do NOT establish reliable prediction of multiple Cirsium photographic appearance dimensions across unseen species and regions. They do not refute the original S1 BIO18 slope, S2 within-stratum colour coefficients, or H3c local analogue resemblance because **the estimands differ**. Before extending to soil/terrain/vegetation, prioritise (a) improved taxon coverage, (b) more geographically independent validation, and (c) measurement error / photographic assessability. New mechanistically distinct layers should be treated as prospective exploratory candidates, not predictors selected to rescue weak RF skill.

## Immutable evidence

- Artifact 13037393509; ZIP SHA256 `5523fd420d9fa4acc89c87e7aa1625285da7f0d998bca8f892a943224090b6f7`.
- `rf_21_vs_23_vpd_wind_full8.json` SHA256 `999cf78d1dfd7d42a981477a3d04a0aa16d9115055eabbf710df230836870bba`.
- `rf_21_vs_23_vpd_wind_posthoc_balanced.json` SHA256 `9f65666fe7925197b4f5c7e11e9e11ed498dd74d8d93b63a91f50553abdb98c4`.
- Only aggregate metrics/receipts exported, not source photo IDs, complete rows, locations, or taxon identities. Original encrypted S1 source and rasters were deleted from the runner after evaluation.

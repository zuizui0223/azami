# Supporting Information: original native-taxon image-angle × BIO18, raw and adjusted evidence

**Scientific freeze:** 2026-10-24 original V4, GitHub Actions [run 40062126006](https://github.com/zuizui0223/azami/actions/runs/40062126006), artifact **12547876704**. Source ZIP SHA-256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`; JSON member `V4_S1_S2_four_vs_vif5_vif10_alltraits.json` SHA-256 `3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70`.

## Why unadjusted estimates belong here, not as another main-text discovery

The primary Part 2 question is whether photograph-referenced capitulum angle covaries with CHELSA BIO18 **among native taxa** in the original simultaneously adjusted BIO18/BIO05/RSDS/elevation V4 model. This is one member of the **80 biological tests** corrected jointly within the native-taxon family. In the **Main Figure 3**, we therefore show **only the Holm-supported between-taxon result**, its standardized slope, its pointwise interval and corrected P value. The original within-taxon term is not suppressed or discarded: it appears below and in the **complete 80-row supplementary CSV and an explicitly supplementary paired forest plot**.

| Estimand | Original n | Standardized β | Pointwise Wald 95% CI | Unadjusted P | BH P across 80 | Holm P across 80 | Familywise-supported? |
|---|---:|---:|---|---:|---:|---:|---|
| **Among native taxa** (equal-taxon means) | 168,907 | **+0.045446843** | [+0.022348959, +0.068544727] | 0.000115079 | 0.009206351 | **0.009206351** | **Yes** |
| **Within native taxa** (taxon-centred) | 168,907 | **+0.058700384** | [+0.002718592, +0.114682176] | **0.039865108** | **0.395836573** | **1.000000000** | **No** |

**Interpretation of the two rows:**
- Both **original point estimates are positive**. The within-taxon point estimate is numerically larger, but has a much wider pointwise confidence interval and is not supported after adjustment. It is therefore **incorrect** to claim that the *estimated slope* itself is generally stronger among than within taxa just because only the among-taxon P value survives correction.
- An unadjusted P of 0.0399 and a pointwise 95% CI excluding zero are consistent with each other, but **neither is a 80-test-family confirmation**. Holm P=1.0 means the within-taxon association did not survive the multiple-comparison family, **not** that its true coefficient is zero.
- These are **within-vs-among observational estimands**, not a formal test of their difference. Identical outcome-specific n does not make them independent samples or direct evidence of adaptation.
- BIO18 means **long-term mean monthly precipitation during the warmest quarter**, not rain at the exact photographed flowering event. The head angle is *relative to the photograph/camera vertical*, not to a measured gravitational vertical, so rainfall shielding/nodding **cannot be concluded** from these fits alone.
- The four-variable analysis was frozen **post hoc** after earlier model exploration; Holm controls the 80 tests in this chosen family, **not** the entire prior environmental-model search.

## Reproducible complete SI and visual source

From the original checksum-locked source ZIP, `plot_bio18_angle_part2.py` writes:

- `figure3_BIO18_among_taxa_image_angle.png` and `.svg` to the main-figure directory (**Holm-supported among-taxon coefficient only**).
- `SI_Figure_BIO18_image_angle_both_scales.png` and `.svg` to the SI directory (both source coefficients, pointwise intervals and raw/BH/Holm context).
- `SI_S1_all_80_original_raw_BH_Holm.csv` (the **full family**, not cherry-picked unadjusted positives).
- `SI_S1_BIO18_image_angle_both_scales.csv` (the two source rows above, not newly fitted slopes).

The existing `build_full_supplement.py` additionally regenerates `Table_S1_full_208_BIO18_tests.csv` (the entire S1=80 plus S2=128 frozen test ledger), independently verifies Holm/BH correction, and extracts the 28 original matched RV pairs. Figure and CSV generation must fail if source ZIP digest, rows or frozen coefficients differ.

No unadjusted-only result is upgraded to a new Holm discovery. **The Main scientific claim remains one positive among-taxon angle×BIO18 association.**
## Original-source geographic sensitivity added after the 24 October freeze

A separate 2026-10-29 BIO18 original-source run recovered exactly **168,907 / 178** records/taxa and the original global beta +0.045446843 to numerical precision. Longitude-sector adjustment reduced the association to **+0.020035** (HC3 pointwise CI −0.010833,+0.050902); spherical geography to **−0.003295** (−0.063000,+0.056410); removing the east-sector taxa to **−0.010372** (−0.068465,+0.047721). With equal 5° cell weight within taxa, the spherical adjustment also included zero.

This **does not change** the original S1 80-test Holm results above. The geographic coefficients have exploratory, **pointwise HC3 confidence intervals**, without a separate correction for the family of geography models. Thus the original Holm-supported global association is **not a demonstration that head orientation tracks rainfall independently of regional species composition**. Earlier BIO12 geographic sensitivity uses different predictors and source cohorts and is not substituted for the present BIO18 results.

[Original-source BIO18 geographic run](https://github.com/zuizui0223/azami/actions/runs/40837107465) • [reanalysis report](../../analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED_JA.md) • [aggregate machine-readable result](../../analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json). The source-locked original within-/among-taxon coefficient forest remains part of SI; main Figure 3 shows the original Holm-supported among-taxon coefficient, while the regional-sensitivity forest remains supplementary.

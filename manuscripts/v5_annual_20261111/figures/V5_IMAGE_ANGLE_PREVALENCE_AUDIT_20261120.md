# V5 photographed-angle flag prevalence audit (2026-11-20)

**Source-verified diagnostic, not an estimated misclassification rate.** The question is: how many originally analysed head-image measurements are potentially difficult to interpret as true upward/downward orientation, given the projection and principal-axis algorithm?

## Completed original S2 79,421-measurement-record audit

- **Successful private source GitHub Actions run:** https://github.com/zuizui0223/azami/actions/runs/44666062718
- Original WCVP native/explicit-introduced labels, three original image-measured head traits, exact original five CHELSA v2.1 climate variables and filtering; exactly 79,421 records matched original V5 analysis strata.
- Original measurements, IDs, user identifiers and private photo coordinates remain in Actions runner. Only aggregate counts are printed/uploaded.
- Comparison unit is an **observation/photo-derived measurement record**, not a distinct JPEG or detected head.

| Diagnostic candidate (not an error label) | N / 79,421 | Fraction |
|---|---:|---:|
| Rotated projected outline aspect near 1 (1/1.10 ≤ aspect ≤ 1.10) | 50,254 | 63.28% |
| Image-based angle > 90° | 25,167 | 31.69% |
| Both near-aspect and angle > 90° | 15,702 | 19.77% |
| Both near-aspect and angle ≥ 110° | 8,786 | 11.06% |
| Aspect within 1/1.05–1.05 | 29,908 | 37.66% |
| More liberal near 1/1.20–1.20 | 69,655 | 87.70% |

The fraction of angle >90° was **31.25% among near-aspect records** and **32.45% among others**. A near-one shape index alone is not evidence that the source angle is wrong; these near-aspect measures are not the PCA eigenvalue anisotropy in the original angle estimator.

### Four verified source strata

| Species × WCVP status | Original V5 N | Near aspect 1.10 N | Both near and angle>90 N (% of stratum) |
|---|---:|---:|---:|
| *C. arvense* × native | 18,505 | 11,707 | 3,942 (21.30%) |
| *C. arvense* × introduced | 22,269 | 14,432 | 5,123 (23.01%) |
| *C. vulgare* × native | 13,115 | 8,229 | 2,231 (17.01%) |
| *C. vulgare* × introduced | 25,532 | 15,886 | 4,406 (17.26%) |

Across equal-record annual-rainfall quintiles in these two species, the joint near-aspect+angle>90 diagnostic fraction ranged 19.27–20.36%; no clear monotonic increase appears in this *unadjusted source-record diagnostic*. This does not exclude photographer, taxon, phenology or spatial confounding, and is not the same population as Part I.

## Interpretation and next checks

The source continuous angle is computed from PCA of foreground projected head coordinates, with the axis signed by the visible floral-pixel centroid and compared with **EXIF-oriented image-up**, not gravity. A longitudinal PCA axis may be uncertain if eigenvalues are similar or the floral-end sign is weak, but the registered **rotated bounding-box aspect ratio** in these records is **not** the eigenvalue ratio and cannot establish this instability on its own.

The *C. pannonicum* (0.732°), *C. mohavense* (111.148°) and *C. subcoriaceum* (113.284°) archived poster examples are illustration-only and are not a representative validation sample. Counting cases near 1 as actual mistakes would overstate the observed evidence.

## Completed original Part-I S1 audit: 157,822 native records

- **Successful authenticated original-source recovery:** https://github.com/zuizui0223/azami/actions/runs/44666968101
- The historic predecoded archive's separate receipt was missing, so the independent remaining SHA-authenticated original July-2026 observation/measurement source and CHELSA 2.1 five climate normals were used. The annual five-axis **157,822** exact full-cohort gate passed. Both the 1° **63 taxa / 7,117 taxon–cells / 153,843 records** and 5° **31 taxa / 814 taxon–cells / 132,998 records** gates also passed. No original source IDs, point coordinates or individuals were published.

| Predeclared photo-measurement proxy | Full S1 N=157,822 | 1° retained N=153,843 | 5° retained N=132,998 |
|---|---:|---:|---:|
| Rotated outline ratio near 1 (0.909...–1.10) | 96,677 (61.26%) | 94,159 (61.20%) | 81,742 (61.46%) |
| Signed image angle >90° | 48,975 (31.03%) | 47,609 (30.95%) | 42,039 (31.61%) |
| **Both near ratio and >90°** | **29,851 (18.91%)** | **29,010 (18.86%)** | **25,565 (19.22%)** |
| Both near ratio and ≥110° | 16,574 (10.50%) | 16,128 (10.48%) | 14,368 (10.80%) |

Among near-ratio source records, 30.88% had image angle >90°; among the others, 31.28% did. The joint diagnostic fraction varied from **18.36% to 19.75%** over five annual rainfall quintiles in the original S1 records; this is an unadjusted record-level summary, not proof of absent taxon/climate-linked measurement error. Importantly **the 5° retained cohort did not exhibit a substantially higher flag fraction than the 1° retained cohort** by this diagnostic. The two cohorts overlap and are not independent.

**Machine-readable audited counts**: `manuscripts/v5_annual_20261111/figures/V5_ANGLE_QC_SOURCE_VERIFIED_COUNTS_20261120.csv`. The public source table records the originating Actions run for every aggregate.

### Decisions and mandatory next verification

These are **candidate flags, not error rates**. Near-one **rotated bounding-box aspect** must not be substituted for PCA **eigenvalue anisotropy**, visible floral-end projection confidence or an independent ground-truth angle. No conclusion can yet be made that 18.91% (or 61.26%) of measurements are wrong.

1. Blind-annotate orientation and image roll on a stratified sample that spans floral shapes, regions, taxa, rainfall, photographic stage and near/outside near-aspect groups; quantify signed-angle disagreement and 0/90/180 categorical agreement separately.
2. Export source-verified distributions of **actual head-level PCA eigenvalue ratio, sign confidence, and quality score**, if retained in raw head receipts. If unavailable, rerun these on licensed original photographs.
3. Re-estimate BIO12–image-angle slopes and original five-axis three-trait climate-analogue results with the same predeclared QC rule, documenting resulting eligible taxa/cells and uncertainty. Run a chroma+outline (angle-excluded) test on identical taxon/cell source support.
4. Address geographic proximity and source taxon eligibility separately; a retained effect in both photo-based grains does not automatically demonstrate climate causation or independently repeated evolution.

The S1 and S2 source-recovery audits above answer the **prevalence** question but do not claim to have recomputed V5 inference after QC exclusion.

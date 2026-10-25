# S2 natural experiment-like contrast: frozen V4 four-axis source audit (25 October 2026)

## Design and choice of the two species

The ecological contrast uses non-random introductions as **historical natural experiments in the limited comparative sense**, not randomized manipulations of rainfall or replicates of adaptive evolution. The source-verified all-350 species-name screen at `analysis/v4/results/20261021_unified_abcde/S2_ALL_350_SPECIES_SOURCE_SCREEN_VERIFIED_20261021.json` found **seven** names with both wild native and explicit WCVP-introduced range statuses. Only **two** passed the predeclared stronger data coverage gate in **both** ranges: ≥100 orientation-eligible photographs/observations, ≥15 broad 5°×sine-latitude equal-spherical-area occupied cells and ≥10 original observers per range. The initial focal species were selected before that all-350 re-screen; the re-screen subsequently confirmed no additional eligible cases. This is data-coverage validation **after a focal choice**, not outcome-blind randomized species selection.

- *C. arvense*: native 78,646 wild source records, introduced 70,964; angle-eligible 49,072 vs 43,195; stronger gate passes in both ranges.
- *C. vulgare*: native 62,368, introduced 85,644; angle-eligible 33,637 vs 47,833; stronger gate passes in both.
- The 156,608 wild introduced records represented by both species constitute **99.91%** of all 156,750 original-source WCVP explicitly introduced wild records. The remaining five both-range species have too few or too spatially confined introduced images. No inference for other thistles follows.

Source/status caution: TDWG/WCVP regional native/introduced flags are coarse, and do not demonstrate photographed-individual ancestry, introduction date, reproductive naturalization, or a mutation. The S1 native photos and S2 native photos overlap. Two species are **not independent evolutionary replication units**.

## Exact original source, frozen BIO18 model

- Frozen model: **BIO18 (warmest-quarter mean monthly precipitation)** + BIO05 (warmest-month maximum temperature) + RSDS (shortwave radiation) + elevation.
- Original experiment code commit: `d2c45b00018a0c04e7c69a29857207dc1ac06d0c`.
- Source-verified success: [Actions #40062126006](https://github.com/zuizui0223/azami/actions/runs/40062126006) / artifact **12547876704**, ZIP SHA-256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`.
- Independent re-extraction checks each source member checksum and enumerates **all 16 camera-angle × four exposure** stratum coefficients in original `S2/four_bio18` Holm family (128 S2 tests overall); **none of the 16 passes Holm**.
- **No new full 4D climate-space density comparison has yet been computed from source-level rows**. The original aggregate holds marginal quartiles/means and climate correlations, not joined site climate records. A 4D centroid, analog fraction, niche expansion or background-standardized D cannot be reconstructed uniquely from four separate medians.

## Results: positive association in source-site occupancy, not proof of niche evolution

| Introduced minus native location median | *C. arvense* | *C. vulgare* |
|---|---:|---:|
| BIO18 (mean mm / month in warmest quarter) | **+48.2** | **+49.7** |
| BIO05 (°C) | **+4.5** | **+3.7** |
| RSDS (original CHELSA radiation scale) | **+3.740** | **+3.688** |
| Elevation (m) | **+115** | **+74** |

Source-level counts differ by trait/environment due to predictor coverage; see `s2_four_axis_source_recheck.py` outputs. All BIO18 site rows were present. The four marginal changes mean the colonized sites are **warmer, wetter in the warmest quarter, sunnier and at higher elevations in median occupancy**. The climatic contrast is NOT a single-axis rainfall addition. Environmental availability and observer sampling may jointly explain the geographic redistribution.

Introduced stratum **four-variable adjusted** camera-angle × BIO18 coefficient:
- *C. arvense* n=24,777, standardized β = **+0.001659**, raw P=.95043, Holm P=1.
- *C. vulgare* n=25,107, standardized β = **−0.022218**, raw P=.55263, Holm P=1.

These are **within-range photo-angle environmental slopes**, not difference-of-range mean-angle or genotype contrasts. The 16 camera-angle coefficients all have Holm P≥0.05, but **absence of statistical support is not an equivalence/nonresponse test**. A prior BIO12 source-appearance comparison reported median 68→67° and 64→63°; do not relabel that historical BIO12 five-exposure comparison as the new BIO18 common-support analysis.

## What to test next using private genuine source rows

The source-level source reconstruction exists in the frozen original workflow `.github/workflows/v4-s1-s2-extreme-lego.yml` at the source execution ref. Once the privacy-safe four-variable join is available, the proper multivariate analysis will:
1. evaluate the complete original BIO18/BIO05/RSDS/elevation four-axis per-site overlap, **separately from photographic measurement eligibility**;
2. use equal-space-cell and observer×cell balancing (before comparing range medians and native-reference scaled climate centroids);
3. quantify 4D niche distance, non-analog introduction-site occupation and cell-resampled uncertainty while distinguishing occupancy from accessible climates;
4. evaluate native versus introduced photographic angles **within matched four-axis common climate support**, including range-status × BIO18 interactions and conditional time/observer/site adjustments;
5. preserve null or nonsupported estimates and 95% intervals and refrain from an equivalence claim without a biologically defined angle-margin and independent physical-camera validation.

This precludes declaring a joint **BIO18 four-axis niche expansion** merely from the verified marginal +48.2/+49.7 mm/month differences.

## GEB framing and relevant literature

Invasion as a natural comparison can be informative but is vulnerable to nonrandom geography, source ancestry and climate availability. Use Guisan et al. (2014; doi:10.1016/j.tree.2014.02.009) for stability/expansion/unfilling and climate analog/background, and Colautti & Lau (2015; doi:10.1111/mec.13162) for differentiation vs actual adaptive evolution. The multi-region native-to-introduced comparison should anchor a **scale-dependent observational claim**, not infer rainfall-driven downward floral evolution or fundamental niche shifts.

**Strong defensible statement:** In native *Cirsium*, warmer-quarter precipitation corresponds to larger camera-frame capitulum angles **among taxa**, while the same two thistle species occupy introduced locations with higher four-axis climate medians, yet exhibit **no common Holm-supported photo-angle gradient** in the fixed strata. This is suggestive of differing ecological scaling or covariance structure, not proof that the same genotypes did not change.

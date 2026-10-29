# Supporting Information: S2 native–introduced floral-colour spatial-block sensitivity (29 October 2026)

**Purpose:** Assess whether the three original introduced-range photograph-derived floral-colour associations, and their three native counterparts, are sensitive to non-independent geography **without adding latitude or longitude as fixed covariates**. This is a separate spatial sensitivity of the original **S2 four-variable BIO18/BIO05/RSDS/elevation model**, not a refit of the original 128-test Holm family or a causal separation of climate and space.

## Source and exact-coefficient gate

- Source: original WCVP-qualified **297,622 wild photographs/observations** from *C. arvense* and *C. vulgare*, in four native/introduced strata. The original confidential observation metadata, WCVP region roster, original photographic colour measurements and five archived CHELSA layers were independently checksum-verified inside a secure Actions runner. The seven additional original CHELSA/WorldClim layers were resampled at original source pixels; missing raster values were **retained, not imputed**. Corresponding individual observations, observer IDs and coordinates were **not published**.
- All **six** original colour–environment coefficients were recomputed on exactly the original **22,936 / 26,171 / 17,335 / 27,200 / 17,335 / 27,200** outcome-eligible observation rows before any spatial sensitivity. Every standardised β reproduced the frozen original to absolute tolerance **2 × 10⁻⁷**.
- Same original month/year nuisance effects where identifiable, original observer × 5° geographic-cell deduplication, inverse-cell WLS weighting, and 5° cell-cluster inference. Fixed environment predictors were **BIO18/BIO05/RSDS/elevation**, without any longitude or latitude fixed effects.
- Additional post-freeze **20° latitude–longitude blocks used only as resampling groups**: 399 nonparametric *whole-block* bootstrap draws; empirical 2.5–97.5 percentile intervals, and source-only **five highest-contributing blocks** removed one at a time. These 20° tiles have varying physical area, are not biogeographic provinces or phylogenetic clades, and are not independent experimental replicates. The original 5° cell-cluster and new 20° block methods represent different approximations to spatial dependence.
- [Original successfully verified GitHub Actions run #40855508122](https://github.com/zuizui0223/azami/actions/runs/40855508122); aggregate artifact **12889487967**; [public aggregate-only machine-readable result](../../analysis/v4/results/20261029_S2_COLOUR_SPATIAL_BLOCK_SOURCE_VERIFIED.json); analysis script `analysis/v4/s2_original_colour_spatial_block_sensitivity.py`. The original 24 October V4 ZIP SHA256 is still `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`.

## Results: exact original coefficients, spatial-block uncertainty

| Species and stratum | Trait × environment | Original n | Frozen β | 20° blocks | Block-bootstrap percentile 95% interval | Percent of draws with original sign |
|---|---|---:|---:|---:|---|---:|
| *C. arvense* native | image lightness × elevation | 22,936 | +0.056267 | 24 | [−0.05208,+0.12724] | 90.0% |
| **C. arvense introduced** | **image lightness × elevation** | **26,171** | **−0.173732** | **10** | **[−0.22949,−0.05578]** | **100.0%** |
| *C. vulgare* native | image chroma × RSDS | 17,335 | +0.136883 | 13 | [+0.06363,+0.22540] | 100.0% |
| **C. vulgare introduced** | **image chroma × RSDS** | **27,200** | **+0.134504** | **14** | **[−0.07581,+0.20119]** | **94.7%** |
| *C. vulgare* native | image lightness × RSDS | 17,335 | −0.048808 | 13 | [−0.11396,+0.12466] | 81.0% |
| **C. vulgare introduced** | **image lightness × RSDS** | **27,200** | **−0.142801** | **14** | **[−0.28093,−0.01521]** | **97.7%** |

All six series completed 399 valid source-block samples. The five largest source-contributing block holdouts retained the original coefficient **sign in 5/5 single-block exclusions for each series**; this is a source-blind influence check, *not* five independent climate replications.

## Interpretation relative to GEB 3.2

- The introduced *C. arvense* elevation–lightness negative and introduced *C. vulgare* radiation–lightness negative were retained under this additional spatial block interval (zero excluded). Their original **Holm P=0.00131233** and **0.01925414** remain valid *within their original fixed regression family*.
- Introduced *C. vulgare* radiation–chroma remained positive in **94.7%** of 20° block bootstrap draws and in all five largest block deletions, but its exploratory percentile interval includes zero (**−0.07581,+0.20119**). Do **not** call this independently spatially confirmed; the original **Holm P=0.01456813** does **not** disappear or become a spatial-bootstrap P-value. Some 20° blocks comprise many of the source's 5° cells and limited geographical support can inflate block uncertainty.
- The native *C. vulgare* radiation–chroma relation is itself positive with a block CI excluding zero (**+0.06363,+0.22540**), consistent with native β+0.136883 being similar to the introduced β+0.134504. It directly reinforces the original warning: **introduced-only Holm detection is not evidence that a new colour–radiation relationship evolved after introduction**.
- Direct native–introduced range × predictor interaction, common-garden genotype contrasts, spectral reflectance, an evolutionary process, realised camera illumination and spatially independent climate contrasts were **not** tested.
- **Neither spatial resampling nor spatially structured predictor adjustment isolates a causal climate component**. The aim is to show sample-geometry dependence without erasing the observed original positive associations. The three introduced Holm-positive observations remain the GEB **Results 3.2** ecological pattern, with the above honest sensitivity qualifications.

This sensitivity is **SI only**; the 24 October original S1 80- and S2 128-test Holm families, Fig. 3/4 original coefficients, and new independent H3c positive local climate-alignment result in Results 3.3 are unchanged.

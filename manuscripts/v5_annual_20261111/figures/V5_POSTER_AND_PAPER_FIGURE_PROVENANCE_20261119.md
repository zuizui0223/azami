# V5 GEB / 58th Species Biology Symposium — figure provenance (19 Nov 2026)

This note accompanies the scientific poster and manuscript working draft in PR #173. It does **not** change the adopted V5 scientific analyses, V4 frozen manuscript or closed exploratory PR #174.

## Shared result order — unchanged

1. **Native climatic analogues and joint capitulum resemblance.** Main-paper Fig. 2.
2. **Individual photographed head traits along different climatic gradients.** Main-paper Fig. 3.
3. **Native–introduced contrasts in two focal species.** Main-paper Fig. 4 (PCA density) and Fig. 5 (trait overlap / chroma).

## Reproducible graphic sources

Run from the repository root:

```bash
python -m pip install numpy matplotlib
python manuscripts/v5_annual_20261111/scripts/build_v5_results_figures.py \
  --repo-root . \
  --output manuscripts/v5_annual_20261111/figures/generated
```

The script reads only **aggregate source-locked CSV files** and writes SVG, PDF and 300-dpi PNG, plus `figure_source_receipt.json` with input SHA256 digests. It has numerical gates preventing accidental substitution of another result family.

| Manuscript use | Frozen source | Figure caution |
| --- | --- | --- |
| Fig. 2: five-axis analogues | `analysis/v5/results/native_seasonality/v5_within_analogue_seasonal.csv`; `v5_cross_analogue_seasonal.csv` | Excess standardised *distance reduction*, **not** regression beta. Exact original Holm P. No original CIs available, so no error bars. |
| Fig. 3: BIO12/angle and BIO15/C* | `analysis/v5/results/native_seasonality/v5_signed_pairs_seasonal.csv` with `sample=all_cells`, `scale=within_pair` | Coefficients 1°/+0.117, 5°/+0.183 and 1°/+0.206, 5°/+0.226; intervals are **pointwise taxon-bootstrap** and exploratory. 1° and 5° reuse the same photos. |
| Fig. 4: climate occupancy | `analysis/v5/results/niche_overlap/V5_S2_niche_overlap.pdf` plus archived density source | Existing plot, must be inspected separately against original climate-density surfaces. Occupancy density D 0.474 / 0.442, outside-native-95%-core mass 18.8% / 44.0%. |
| Fig. 5: head space and chroma | `analysis/v5/results/novel_head_overlap/novel_head_overlap.csv`; `analysis/v5/results/invasion_design/climate_balanced_heads.csv` | Geography-matched head D 0.740→0.735 / 0.847→0.789 with two-species Holm P 0.953 / 0.096. C. vulgare introduced-minus-native mean C* +1.741; climate-class-reweighted +2.465: **descriptive only**. |
| Source photo measurement panel | `legacy/v3/reproducibility/figures/source/Figure_1_real_photo_measurement_provenance.png` and `.csv` | Real source photography and real computer-vision overlays. Orientation is relative to **image vertical**, not gravity. |
| **V5 three-trait photo coverage maps (replaces V4 proxy)** | `manuscripts/v5_annual_20261111/figures/v5_photo_coverage_5deg_public.csv` generated from authenticated confidential source in Actions run **44514111389**; `scripts/build_v5_public_photo_coverage_maps.py` | **Exact V5 79,421 measurement-record cohort**, WCVP native/introduced and five CHELSA annual variables complete. 350 public 5° cells represent 79,141 records (n≥5 per cell); 280 records suppressed. Points are **cell centres**, not original coordinates or individual photos. |
| Superseded V4 contextual source map (do not mix with V5 in poster) | `analysis/v4/results/20261030_SAMPLING_COVERAGE_SOURCE_VERIFIED.json` with `manuscripts/v4_s1_s2_frozen_20261025/plot_sampling_coverage_v4.py` | The map uses **real 5° cells** from a V4 photo-measurement cohort, not fictitious dots. **It is not an exact V5 three-trait-complete 79,421-record location map.** Do not label it as such in the paper or poster. |

## Photographic attribution

- *Cirsium pannonicum*: © Bobo-X, source iNaturalist photograph 532420148; CC BY 4.0. Annotation/cropping is a modification.
- *C. mohavense*: Jesse Rorabaugh, iNaturalist photo 4051162, CC0.
- *C. subcoriaceum*: Diego Blanco, iNaturalist photo 475732253, CC0.

Credits and original licences are recorded in `manuscripts/v4_s1_s2_frozen_20261025/FIGURE1_PHOTO_CREDITS_20261030.json`.

## Ecological interpretation (poster, not tested function)

The positive BIO12–photo-angle association is **consistent with** a pollen-wetting or rain-shedding hypothesis (cf. Chen et al., 2013 on *Cremanthodium*), but the V5 dataset does not measure capitulum gravitational inclination, pollen protection or fecundity. The positive BIO15–photographic C* association is **consistent with** seasonal water-related colour or pigmentation hypotheses (cf. Warren & Mackenzie, 2001, including *C. palustre*) but neither flower-stage drought nor pigment concentration was measured. Climatic similarity among taxon samples does not establish phylogenetically independent parallel evolution, genetic fixation or speciation.

## Remaining data gate

**Completed (2026-11-19):** The original-source workflow produced and committed the V5-specific photo-eligible **5° cell count table** for the four native/introduced species×range strata (Actions run **44514111389**, successful). The exact cohort totals are *C. arvense* 18,505 native / 22,269 introduced and *C. vulgare* 13,115 native / 25,532 introduced: **79,421 measurement records**. Privacy filtering leaves 350 publicly displayable cells covering 79,141 records, with 280 records withheld in low-frequency cells. The poster now uses these **V5** cells, not the earlier V4 proxy. Their plotted centres are illustrative symbols of true **coarse cell membership**, never source photo coordinates. This completes a spatial *coverage figure*, not new statistical verification of photo measurement accuracy.

## Existing draft vector artwork

The current manuscript branch also retains `fig2_native_annual_resemblance.svg`, `fig3_annual_trait_gradients.svg`, `fig5_native_introduced_chroma.svg` as earlier source-driven artwork. The new generator provides a repeatable independent rendering route for those same **frozen result values** rather than claiming a new statistical result.

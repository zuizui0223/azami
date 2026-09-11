# Chapter 1 v3 biological-construct sensitivity chain — 2026-09-10

## Status

Complete. GitHub Actions run `34418904597` passed the complete construct-level sequence using the frozen Chapter 1 v2 source cohort and the same v2 sensitivity logic: sampling composition -> broad/residual spatial sensitivity -> historical-placement sensitivity.

This is a biological-construct robustness layer. The frozen v2 endpoint-level atlas, multiplicity family and manuscript conclusions remain canonical.

## Starting point

The 22 already-measured v2 endpoints were reorganized without using environmental outcomes into 11 named constructs: 10 inferential constructs plus descriptive surface specularity. The construct atlas used the same nine frozen environmental gradients. At the primary among-taxon min5 scale, 12 construct-gradient rows passed the reduced construct-family BH correction; 9 within-taxon rows passed their separate construct-family BH correction.

## Sampling-composition gate

The v2 sampling scenarios were reapplied: omission of each of the ten most-observed taxa, joint omission of the top two, leave-one-broad-region-out, native-only restriction, and equal-total-weight-per-taxon for within-taxon rows.

All 12 selected among-taxon construct-gradient rows retained direction/vector alignment in every declared sampling scenario. Six of nine selected within-taxon rows were stable in every declared scenario; three hue rows reversed vector direction in at least one scenario. Overall, 18 of 21 selected rows were stable under all declared sampling scenarios.

For the two frozen v2 headline patterns:

- floral chroma x shortwave radiation: all sampling scenarios stable; minimum effect-magnitude ratio = `0.614876`.
- presentation angle x annual precipitation (BIO12): all sampling scenarios stable; minimum effect-magnitude ratio = `0.780783`.

## Broad-space and residual-spatial gate

The same second-order spherical-coordinate basis, Freedman-Lane permutation logic and residual Moran-I screen used by v2 were applied to the construct-level rows.

Among-taxon: 12 rows entered and 3 passed.

1. **floral chroma x shortwave radiation**
   - baseline beta = `-0.3453720171`
   - spatial beta = `-0.7124111818`
   - spatial permutation P = `0.005`
   - residual Moran I = `-0.047171`
   - residual Moran P = `0.188`
   - pass = `true`

2. **presentation angle x annual precipitation (BIO12)**
   - baseline beta = `+0.3043592859`
   - spatial beta = `+0.2860860885`
   - spatial permutation P = `0.023`
   - residual Moran I = `-0.066578`
   - residual Moran P = `0.072`
   - pass = `true`

3. **floral chroma x NPP**
   - baseline beta = `+0.2463466832`
   - spatial beta = `+0.4136986109`
   - spatial permutation P = `0.001`
   - residual Moran I = `-0.031772`
   - residual Moran P = `0.418`
   - pass = `true`

Presentation angle x growing-season precipitation retained a positive coefficient but did not pass the broad-space permutation gate (`P = 0.171`). The selected hue rows did not pass the residual-spatial screen; several retained spatial associations but residual Moran P values were approximately `0.001–0.002`.

Within-taxon: 9 rows entered and 0 passed the complete broad-space plus residual-Moran gate. This does not retroactively alter the frozen v2 endpoint-level within-taxon results; it shows that no reduced construct-level within-taxon row met the stricter full construct-level spatial gate in this secondary synthesis.

## Historical-placement gate

Only the three among-taxon spatial passes entered the frozen 52-tree placement sensitivity (scenario 1, scenario 3 and 50 randomized scenario-2 placements) using Pagel-lambda GLS/PGLS.

All three passed on all 52 trees:

- **floral chroma x shortwave radiation:** 52/52 P < 0.05; P approximately `2.39e-05` on all trees; lambda = `0` throughout; coefficient direction stable.
- **presentation angle x annual precipitation:** 52/52 P < 0.05; P range `0.000201–0.000231`; lambda range `0–0.053690`; coefficient direction stable.
- **floral chroma x NPP:** 52/52 P < 0.05; P range `0.002080–0.003993`; lambda range `0–0.212192`; coefficient direction stable.

The 52 trees are alternative historical-placement sensitivity scenarios, not 52 independent biological confirmations.

## Interpretation relative to frozen v2

The central result is unchanged and strengthened: **both frozen v2 headline among-taxon patterns survive the full robustness sequence after the 22 image endpoints are reorganized into biologically interpretable constructs.** Thus their survival is not restricted to the original endpoint-level parameterization.

`floral_chroma x NPP` is a new v3-only exploratory pattern. It became FDR-supported only in the smaller construct-level multiplicity family and therefore must not be promoted as a replacement for, or expansion of, the frozen v2 headline conclusions without a separately declared confirmatory analysis.

The endpoint-level bract projection-peak-density x VPD association did not enter the construct-level chain because the reaggregated projection construct was not FDR-supported for VPD. This suggests that the original VPD result was more endpoint-specific than the two headline patterns.

## Reproducibility receipt

- successful run: `34418904597`
- workflow job: `102689810335`
- branch head at run: `3cc65f71b25abfb5bd730649ec36b69a26672f7f`
- Actions artifact: `10130210432`
- artifact digest: `sha256:7568ec98b0709f6cda027b928298d58c7e6aa42af0f6c0947e1e91c3d1abbab2`
- frozen continuous-trait artifact: `9612943217`
- frozen process-environment artifact: `9633419268`
- frozen broad-region artifact: `8983877726`
- frozen historical-tree artifact: `8227254443`
- frozen native-status table recovered from immutable tag `azami-ch1-v2-2026-08-27`; Git transport LF was normalized to CRLF and then matched the frozen SHA-256 `c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a` exactly.

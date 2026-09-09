# Full-cohort VIFstep + multivariable sensitivity (2026-09-09)

## Status

This is a **retrospective post-hoc sensitivity** layered on the frozen GEB-v2 atlas. It does not replace the frozen marginal one-predictor results and it does not identify causal environmental effects.

The analysis uses the full strict-spatial cohort (**46,276 observations / 259 taxa**) and all measured registered endpoints (22 measured endpoints = 21 measured inferential units after treating hue sine/cosine jointly). No native-range filter is used.

## Input validation

The recovered 46,276-row nine-predictor environment was checked by rerunning the repository's frozen `among_taxon_min5` univariate code. All **234 rows** matched the frozen atlas in status; all **180 successful rows** reproduced the frozen p values exactly, with coefficient/q-value differences only at floating-point epsilon. The frozen **10 FDR signals** were recovered exactly before the new sensitivity was fitted.

## Why the final VIFstep is response-cohort specific

An initial common-set VIFstep on all 259 taxon environmental medians removed BIO1 at threshold 10 and BIO1 + GSP at threshold 5. That common set was not strict enough after response-specific taxon availability was applied: endpoint subsets could again have VIF above the requested threshold. The final analysis therefore recomputes VIFstep inside each endpoint's actual taxon cohort before model fitting.

This closes the loop: **failure location -> subset-specific collinearity -> endpoint-cohort VIFstep -> every successful model below its threshold**.

For `among_taxon_min5`, thresholds 10 and 5 converge to the same effective models for the primary colour/orientation and most other endpoints. There are 137 successful unit-predictor tests across 20 measurable units, and the maximum final VIF is **3.746**.

## Strict result

Only one among-taxon min5 association survives simultaneous adjustment and global BH-FDR at both VIF thresholds:

- `corolla_lab_chroma ~ chelsa_rsds_mean`
  - n taxa = **143**
  - adjusted standardized beta = **-0.467430**
  - HC3 95% CI = **[-0.705230, -0.229630]**
  - Freedman-Lane permutation P = **0.0003**
  - BH q = **0.0411**
  - final max VIF = **2.643**

The frozen marginal coefficient for the same association was beta = -0.345372, q = 0.006. Thus the negative chroma-radiation relationship is not explained away by the other retained environmental gradients and becomes stronger in the simultaneous model.

The other frozen min5 signals do not survive the strict sensitivity:

- hue × BIO1 is removed by VIFstep;
- orientation × GSP is removed by VIFstep;
- orientation × BIO12 remains positive (beta = 0.266823) but CI crosses zero and q = 0.229475;
- bract projection peak density × VPD remains positive (beta = 0.169031) but q = 0.793458;
- the retained joint hue associations are not FDR-supported after simultaneous adjustment (e.g. radiation joint P = 0.0048 but q = 0.20824).

At `among_taxon_min2`, chroma × radiation is again the only FDR-supported result at both VIF thresholds. A chroma × GSP coefficient appears under VIF<10 but is not FDR-supported and disappears from the stricter VIF<5 predictor set, so it is not a stable adjusted result.

## Interpretation for Chapter 1

The original atlas should continue to be described as a map of **marginal associations along correlated observational gradients**. This sensitivity sharpens the hierarchy among those patterns: **lower chroma under higher shortwave radiation is the only among-taxon environmental association that remains FDR-supported after explicit multivariable collinearity control across the full measured trait universe.** Orientation/precipitation, hue and VPD/involucre signals remain useful marginal or hypothesis-generating patterns, but they should not be described as environmentally independent effects.

## Reproduction

Run `analysis/v3/run_full27_vif_multivariable_sensitivity.py` with the continuous trait universe and a 46,276-row nine-predictor environment table. Defaults use 9,999 Freedman-Lane permutations, 1,000 taxon bootstraps for the circular joint response, VIF thresholds 10 and 5, and BH correction across every successful tested unit-predictor row separately for min5/min2 and each VIF threshold.

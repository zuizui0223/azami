# V5 photographic image-angle measurement error: frozen analysis sensitivity protocol

**Date:** 2026-11-20. **Purpose:** quantify whether uncertainty in automatically measured photographic capitulum angle changes V5's central five-axis climatic-analogue result and exploratory BIO12–image-angle association. This is methodological stress testing, not validation of gravity-based biological orientation.

## Original analysis and source gates

Input: exact source-recovered, authenticated native V5 **157,822 observation/photo-derived measurements**, five climate axes (BIO01, BIO04, BIO12, BIO15 and RSDS), three co-recorded visible photographic traits (signed image-axis angle, CIELAB chroma, rotated 2D outline aspect). The original 1° and 5° samples are 63 and 31 taxa, 7,117 and 814 taxon–cells, respectively. Across both grains, retain the exact source aggregation (taxon × geographic cell median), taxon-balanced pooled standardisation, original 24 eight-cell rarefaction draws, original 499 permutations of *intact phenotypic vectors within taxon*, and original within- versus between-taxon correction families.

Baseline hard gates before reporting a sensitivity: within-taxon excess = +0.0842122783 at 1°, +0.2383692461 at 5°; cross-taxon excess +0.0309379488 and +0.0653295626; and BIO12–photo-angle all-cell within-pair standardized coefficients approximately +0.1171 and +0.1829. If the gates fail, report a source/algorithm discrepancy rather than treating the new value as a sensitivity effect.

## Distinct errors, distinct estimands

1. **Angle deletion, chroma + outline**: same source cohort, same climate analogue pairs, original within-taxon whole-vector null, but a two-dimensional rather than three-dimensional phenotype distance. Test source-angle dependence. Its effect size need not be numerically on the same scale as the 3D result; compare evidence/uncertainty, not just excess differences.
2. **Small symmetric technical perturbation**: synthetic iid Gaussian angle jitter at the original *photographic observation record* level, sigma = 4.67°/1.95996. Calibrated to the **historical post-QC horizontal-mirror absolute discrepancy p95** across a different technical audit; not an empirical V5 record-specific error model.
3. **Strong symmetric perturbation**: sigma = 54.1°/1.95996 based on **deliberately shifted 5% detector boxes** in a separate 100-head perturbation audit; a strong hypothetical challenge, not an estimate that the typical V5 photo is wrong by 54°.
4. **Differential environmental bias stress**: synthetic photo-angle offset of ±10° per within-taxon BIO12 standard deviation (winsorized z to ±2), holding source observations and climate strata fixed. **No evidence says the actual error follows these values.** This is a failure-mode assessment because an error correlated with climate can create/mask a slope or alter the climatic-analogue test.

Recompute source-record angle after perturbation and then recompute taxon–cell medians, pooled standardization and the original within-pair signed slope (including its original taxon bootstrap pointwise interval), rather than editing the published aggregate coefficient directly. Repeat the two random scenarios with independent seeds; regard replicate ranges as Monte Carlo sensitivity, not biological confidence intervals. Record original 1° and 5° as **linked views of the same photo sources**, not independent external replications.

## Interpretation and manuscript decision

- If small symmetric-error perturbation leaves the source-specific null-tested result and the signed slopes similar, say *insensitive to the stated synthetic random noise*, **not** that image angles are biologically accurate.
- If the strong detector-shift scenario destroys the 3D result or the signed gradient, report that sensitivity; do not reinterpret this as a measured field error prevalence.
- If the 2D chroma + outline analysis retains permutation support, the *climatic resemblance result does not require the angle coordinate*, though it remains subject to geography-preserving null and image validity checks.
- If a ±10°/SD rainfall-linked offset changes the sign or inference of BIO12–image-angle slopes, be explicit: the result is vulnerable to differential measurement bias. If it is stable, other untested camera/taxon/observer biases may still exist.
- **Primary caveat:** within-taxon phenotype permutation is not a geography-preserving null. Spatial autocorrelation and climate-driven photographic sampling remain separate confounds requiring matching/restricted permutation.
- A smaller P value at 5° is not automatically a resolution effect because the 5° cohort contains fewer, more widely sampled taxa; the same-31-taxon analysis is a separate gate.

## Existing results that should NOT be described as V5 error validation

The historical v3 frozen synthetic angle perturbation preserved a positive *different* BIO12–angle regression in 2,000/2,000 iid-error replicates under both 4.67° and 54.1° p95 scenarios. That analysis operated on 142 taxa/37,042 original v2/v3 observations and a different model (unperturbed standardized beta +0.304359); it cannot stand in for V5's annual five-axis within-pair estimates.

The V5 original-source QC prevalence report (see \`ANGLE_MEASUREMENT_ERROR_AUDIT_20261120.md\`) showed 61.3% of S1 observations within a factor 1.10 of rotated **outline bounding-box aspect** 1, and 31.0% angle >90°, but no meaningful enrichment of >90° among near-unity outline ratios. This does not estimate PCA eigenvalue anisotropy, camera roll or angle ground-truth error.

## Additional gold-standard validation that sensitivity analyses cannot replace

Select licensed, original photographic records stratified by angle, outline, taxon, wet/dry climate and photographer; have independent annotators mark the *floral-facing head longitudinal axis in an EXIF-oriented photograph* and, separately where possible, peduncle and actual gravity reference. Record per-image angular absolute/circular deviation, P90/P95 error, near-horizontal disagreement, detector-box variation and loss by assessability class. Enforce photograph- and taxon-cluster uncertainty. Use a held-out annotation set before any final claim of rain-protective nodding.

**Reproducibility code:** \`analysis/v5/v5_angle_error_sensitivity.py\` via \`analysis/v5/run_angle_suspects_s1_recovered.py\` on the encrypted-source Actions workflow \`.github/workflows/v5-angle-suspects-s1-recovered.yml\`. All per-photo identities, measured values and coordinates remain private.

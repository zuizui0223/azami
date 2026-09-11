# Chapter 1 v3 assessability and technical-stress audit — 2026-09-10

## Purpose

This layer addresses two reviewer-facing failure modes without changing the frozen v2 endpoint atlas, its multiplicity families, or its two manuscript-level positive headline associations:

1. whether trait **measurement opportunity itself** changes along the same environmental gradients used in the ecological atlas;
2. whether the orientation–annual-precipitation association is fragile to technical perturbation magnitudes preserved in the frozen image audit.

The first analysis is outcome-blind: trait values are never responses. The second is an explicitly synthetic stress envelope, not an empirical measurement-error model.

## Outcome-blind assessability audit

The frozen v2 universe contains 46,276 observations, 259 source-assigned taxa and 22 measured endpoints. For every endpoint and each of the nine frozen environmental predictors, measurement availability was regressed on a taxon-centred predictor. The coefficient is a change in measurement probability per one within-taxon SD of environment; uncertainty uses a taxon-cluster sandwich estimate.

- endpoint tests: **198**;
- endpoint rows passing BH FDR: **17**;
- construct tests: **99**;
- construct rows passing BH FDR: **11**;
- largest absolute availability gradient in either family: **1.024 percentage points per within-taxon SD**.

Thus measurement opportunity is not perfectly environment-independent, but the detected gradients are small in absolute probability scale. Statistical support should not be confused with a large selection effect because the observation denominator is 46,276.

### Frozen headline 1: floral chroma × shortwave radiation

- frozen trait association: beta = **-0.345372**;
- chroma measurement availability = **87.54%**;
- assessability slope versus radiation = **+0.873 percentage points per within-taxon SD**;
- raw P = **8.18e-09**;
- BH q across the 198 endpoint assessability tests = **3.73e-07**.

The availability gradient is therefore statistically detectable but numerically small and has the **opposite numeric direction** to the frozen trait-value association: chroma becomes slightly more measurable, not less measurable, toward higher radiation. This does not prove absence of selection bias because value-dependent missingness remains possible, but simple environment-dependent availability does not directly reproduce the negative chroma–radiation slope.

### Frozen headline 2: presentation angle × annual precipitation

- frozen trait association: beta = **+0.304359**;
- orientation measurement availability = **80.50%**;
- assessability slope versus annual precipitation = **+0.688 percentage points per within-taxon SD**;
- raw P = **0.01246**;
- BH q across the 198 endpoint assessability tests = **0.11065**.

The availability slope is in the same numeric direction as the trait association, so it is retained as a caution rather than dismissed. However, its magnitude is below one percentage point per within-taxon SD and it does **not** pass the global endpoint-level assessability FDR family.

## Frozen-summary orientation technical stress

The frozen automated audit preserves two value-scale orientation discrepancies:

- post-QC horizontal-mirror discrepancy p95 = **4.67 degrees**;
- discrepancy p95 under an intentional **5% bounding-box shift = 54.1 degrees**.

The audit does not preserve per-image error draws. Therefore each value was used only to calibrate an explicit zero-mean independent Gaussian stress convention whose absolute-error p95 equals the frozen summary. Orientation values were perturbed at observation level, clipped to 0–180 degrees, taxon medians were recomputed, and the frozen among-taxon orientation–BIO12 standardized slope was recalculated in **2,000 replicates** per scenario.

The unperturbed reconstruction exactly reproduced the frozen coefficient:

- **beta = 0.3043592858977514** across **142 taxa**.

### Mirror-p95-calibrated routine stress

- median beta = **0.30294**;
- 95% simulation interval = **0.28848–0.31642**;
- median effect ratio to unperturbed = **0.995**;
- positive/sign-retaining replicates = **2,000/2,000**.

### 5%-bbox-shift-p95-calibrated severe stress

- median beta = **0.25481**;
- 95% simulation interval = **0.15497–0.35180**;
- median effect ratio to unperturbed = **0.837**;
- positive/sign-retaining replicates = **2,000/2,000**.

This shows that the positive orientation–precipitation direction is not fragile to **independent symmetric random perturbation** at either frozen technical scale. The severe scenario is deliberately not presented as the actual field error distribution.

## What remains unresolved

### Visible colour

The frozen perturbation summary records colour as **generally stable**, but it preserves no numeric colour error distribution or per-image perturbation values. Quantitative chroma error propagation was therefore **not run**. Inventing a chroma error distribution would create stronger-looking but unsupported evidence. Independent calibrated visible/UV colour validation remains an external completion gate.

### Orientation

The stress test does not address:

- systematic camera roll;
- environment-dependent measurement error;
- detector bounding-box selection correlated with precipitation or phenotype;
- gravity-referenced botanical accuracy.

Thus it supports technical **sign robustness under explicit random stress**, not accuracy or mechanism.

### Selection

The assessability audit uses availability as the outcome and is therefore independent of measured trait values. It cannot rule out missingness that depends jointly on environment and the unobserved trait value.

## Net effect on the v3 story

**Positive reinforcement**

- the orientation–precipitation sign survives both routine and deliberately severe frozen-summary-calibrated random stress in 100% of 2,000 replicates;
- simple chroma availability changes only slightly with radiation and in the opposite numeric direction to the negative chroma–radiation trait slope;
- the orientation availability gradient is also small and is not FDR-supported across the 198 endpoint assessability tests.

**Caution retained**

- assessability is not perfectly random: 17/198 endpoint rows and 11/99 construct rows are FDR-supported;
- therefore the manuscript should report trait-specific assessability rather than claim an unbiased probability sample.

**Unresolved external validation**

- no frozen numeric colour-error distribution exists for quantitative chroma propagation;
- technical orientation stress is not gravity-referenced biological validation.

None of these diagnostics changes the frozen v2 positive results or their claim level. They narrow alternative technical explanations and make the remaining external validation requirements explicit.

## Reproducibility

- workflow run: `34434467062`
- job: `102736566352`
- artifact: `10135679053`
- artifact digest: `sha256:fa80f6079d676186d48665a574b5e68f906037a466a9e3101a57d87ee02f5ba7`
- orientation stress replicates: `2000` per scenario
- frozen trait artifact: `9612943217`
- frozen environment artifact: `9633419268`
- frozen technical audit: `analysis/ch1/image_to_trait_automated_technical_audit_summary.json`

Claim boundary: post-hoc reviewer-response diagnostics only; frozen v2 inference remains canonical.

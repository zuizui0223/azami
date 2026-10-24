# Chapter 1: multitrait ecological synthesis and covariance across scales (2026-10-24)

Status: **new, post-hoc exploratory extension**. This note does not revise the frozen GEB
numerical reference, establish new headline statistical support, or alter the blinded
manuscript (which is stored outside GitHub). The proposed re-framing below is an
**interpretive option**, not a retrospective claim that these tests were preregistered.

## Central ecological question

**How are visible capitulum presentation, corolla colour, capitulum form, and
involucral projection traits assembled across taxa and within taxa, and to what
degree are their observed associations explained by shared environmental gradients?**

Rather than make precipitation x orientation the whole story, show:

1. **Trait–environment coverage:** all nine biological constructs and nine frozen
   environmental predictors at each scale, not only significant pairs.
2. **Trait–trait geometry:** all 36 pairs on exactly the same 1,734 observations /
   42 taxa, with within vs among strength shown side-by-side.
3. **Environmental conditioning:** retain the original raw RV values, then compare
   RV after partialling out three outcome-blind PCs from the original nine
   environmental variables. This is exploratory; it does not identify mediation.
4. **Estimator sensitivity:** do not equate raw 33/36 to a bias-robust count.
   Existing equal-n sensitivity gives median 23/36 and null-centering 19/36.
   Include equal-n comparisons for climate-conditioned RV as well.
5. **Evidence tiers:** only the two frozen ecological anchors deserve the full
   "survives sampling, spatial/residual and historical placement" status.

## Evidence read directly from frozen current-reference CSVs

Cohort for these pairwise RV comparisons: 1,734 observations, 42 taxa; nine
constructs; 36 pairs. RV coefficients are **unsigned association strengths**.

| Construct pair | within-taxon RV | among-taxon RV | difference | Interpretation priority |
|---|---:|---:|---:|---|
| presentation angle × head elongation | 0.002343 | 0.339617 | +0.337274 | geometric coupling candidate; camera pose is a major alternative |
| floral chroma × projection prominence | 0.000166 | 0.334516 | +0.334349 | colour–armature coupling candidate; no demonstrated defence efficacy |
| floral lightness × floral hue | 0.184336 | 0.465806 | +0.281469 | strong colour-module coherence; colour errors may be shared |
| floral chroma × floral hue | 0.033916 | 0.290331 | +0.256415 | colour-module pattern |
| presentation angle × floral hue | 0.014833 | 0.196294 | +0.181461 | potential presentation–colour contrast |
| projection prominence × projection pattern | 0.350738 | 0.099930 | -0.250808 | exception: stronger raw within-taxon co-variation |

**These rows were selected after looking at all 36 results. No pair-specific
FDR, permutation or bootstrap tests are attached to these numeric comparisons.
Never present the ranking as a confirmatory discovery list.**

Original scale-wide results:
- raw median RV: within **0.002238**, among **0.043212**, 33/36 raw stronger among;
- common-cohort matrix alignment Spearman rho **0.439125**, QAP P **0.0041**;
- 1,000-taxon bootstrap for global contrast: within-vs-among strength difference
  95% interval **+0.038197 to +0.119048** (this is the original estimator);
- equal-n estimator audit: median stronger-among count **23/36**, overall median
  difference 95% interval **+0.00188372 to +0.03035156**;
- null-centered pair count **19/36** (diagnostic only).

Sources:
- reproducibility/current_reference/upgrade/complete18_construct_pairwise.csv
- reproducibility/current_reference/upgrade/construct_scale_contrast_summary.json
- reproducibility/current_reference/estimator_validity/rv_estimator_validity_summary.json

## Trait-environment atlas and inferential boundaries

Frozen v3 among-taxon associations on the broader construct-specific cohorts:
- presentation angle × BIO12: beta +0.304359; BH q 0.00360; passes the full
  sampling -> spatial/residual -> 52-tree placement sensitivity.
- floral chroma × shortwave radiation: beta -0.345372; BH q 0.00225; likewise
  passes the full robustness sequence.
- lightness × GSP, chroma × NPP, projection pattern × BIO4: v3 exploratory
  BH associations, not new frozen-v2 anchors.
- hue: several BH associations in the ordinary atlas, but the hue rows fail
  the current residual-spatial robustness gate.
- head elongation and head compactness: no among-taxon FDR-supported
  single-predictor associations in the frozen nine-predictor atlas.

**Do not interpret an unsigned RV as an environmentally directed trait change.**
The existing construct environmental atlas is univariate/marginal; even FDR
support there is not an independent predictor effect from a joint model.

Sources:
- reproducibility/current_reference/axes/biological_axes_among_min5.csv
- analysis/v3/BIOLOGICAL_AXIS_RESULT_SUMMARY_20260910.md
- reproducibility/CURRENT_FIGURE_SURFACE_20260912.md

## Added executable sensitivity: first stage

Run the post-hoc sensitivity against the **same frozen inputs** and confirm it
reconstructs all 36 raw RV edges before calculating anything new:

```bash
python -m analysis.v3.run_multitrait_environment_conditional \
  --traits work/current/inputs/traits.csv \
  --environment work/current/inputs/environment.csv \
  --out-dir work/multitrait/results \
  --environment-pcs 3 \
  --matched-replicates 500 \
  --seed 20261024
```

The named GitHub Actions workflow recovers the original frozen numerical archives,
verifies checksums, runs synthetic unit tests, then writes a conditional-pair
CSV, matched-n replicates and a machine-readable report. Until that run succeeds
the climate-conditioned results are **not reported**.

The adjustment uses 3 PCs fitted **to environment, without trait outcomes**.
These PCs are orthogonal but this does not solve all climate correlation or
geographic confounding. Conditioning may also remove ecologically relevant
shared variation, and the effective meaning of each PC changes by scale.

## BIO18, elevation and bioclimatic extremes: stage 2 (not yet computed)

BIO18 (warmest-quarter precipitation) and elevation are **absent** from the
current nine-variable process environment. Their covariances and effects
cannot be inferred from the frozen RV/atlas tables.

Predeclare a *candidate-pool* diagnostic with continuous climate and terrain:
- original BIO1, BIO4, BIO12, BIO15, shortwave radiation, VPD, wind, GSP, NPP;
- potential BIO5/6 (warmest/coldest extremes); BIO13/14 (wettest/driest month);
  BIO17/18 (driest/warmest quarter precipitation);
- elevation (plus terrain complexity only with a distinct biological question).

For exact same mapped observations, first produce Pearson/Spearman correlation
matrices, VIF, effective rank, missingness and nonlinearity by **within-taxon
demeaned**, **taxon-median**, and raw levels. Example pre-existing warning:
the older 9-variable audit gives BIO1 ~ VPD Pearson r 0.872 and a joint-model
VIF of 29.38 for BIO1 / 21.00 for VPD, so do NOT put all correlated layers
uncritically in one regression.

The subsequent tests must compare *substantively distinct models*:
- annual rainfall (BIO12) vs warmest-quarter rainfall (BIO18), on the same taxa;
- seasonal/extreme water variables vs annual totals;
- with and without elevation, while comparing collinearity and effective sample;
- climate-associated trait structure before and after controlling geographic
  or phylogenetic structure, preserving separate uncertainty statements;
- at least a leave-region-out or balanced-biogeographic sensitivity.

Model selection using the traits must not be relabelled as preregistered.
Exploratory new predictors require a **new multiplicity family** and cannot
inherit the old q-values or spatial/tree gates.

## Revised main-paper figure logic (proposal)

1. Measurement definitions and trait assessability.
2. Global coverage and biases.
3. **Within/among covariance of all 9 constructs**, accompanied by equal-n/
   null-centering cautions; colour, form and armature fully visible.
4. **All-construct x environment atlas**: frozen marginal coefficients with
   explicit evidence-tier and sample-size labels.
5. **Robustness of the two retained ecological anchors**. If new climate
   variables yield robust independent results, they require a *separate*
   prospective or explicitly exploratory reanalysis before Main promotion.

This retains the core v3 result (scale-dependent integration) and uses the
existing colour and form measurements rather than treating their non-significant
individual climate slopes as failed data. It does not imply adaptation,
genetic/developmental modularity or physical trait measurement accuracy.


## Executed first-stage results — 2026-10-24 (new, exploratory)

**GitHub Actions succeeded**: run `40055211966`, job `120756591327`,
artifact `12544872267` (`ch1-multitrait-climate-conditional-40055211966`);
artifact ZIP SHA-256
`155bf1bc7734829b543309f50666d37ec4b030503d758d312d3de44a8648df39`.
Input archive hashes, independent unit tests and exact reconstitution of all
36 baseline raw RV rows passed.

The first three environmental PCs explain **62.04%** of the weighted within-taxon
nine-predictor variance and **82.82%** of the among-taxon nine-predictor
variance. These are distinct scales and distinct PC bases.

| Diagnostic | Before climate PC adjustment | After 3-PC adjustment |
|---|---:|---:|
| median within RV | 0.002238 | 0.002286 |
| median among RV | 0.043212 | 0.045955 |
| pairs stronger among (raw RV, descriptive only) | 33/36 | 34/36 |
| median stronger-among count from equal-n resampling | 23/36 (older estimator audit) | 24/36 (500 new replicates) |
| median among-minus-within matched-n RV difference | +0.019316 (older audit) | **+0.021654** |
| 95% equal-n resampling range of difference | +0.001884 to +0.030352 | **+0.002326 to +0.033197** |

In the new adjustment, **98.2%** of matched-n resamples have a positive median
among-minus-within difference. These are Monte Carlo resampling *distributions
conditional on the fixed taxa and among-taxon medians*, not inferential
confidence intervals for a new population of taxa or independent confirmations.

Three high-difference examples after environmental adjustment:

| Pair | within adjusted RV | among adjusted RV | matched-n difference median |
|---|---:|---:|---:|
| presentation angle × head elongation | 0.002414 | 0.365192 | +0.351387 |
| floral chroma × projection prominence | 0.000328 | 0.311041 | +0.296904 |
| floral lightness × floral hue | 0.181528 | 0.485623 | +0.291583 |

Exception:
projection prominence × projection pattern has adjusted within RV
**0.351340** and among RV **0.067699**, matched-n difference median
**-0.273725**.

**Safe interpretation:** The qualitative tendency for stronger among-taxon
visible-trait integration is not eliminated by this particular linear,
three-PC shared-climate adjustment, after matching the number of rows at the
two scales. The contrasting bract-projection relation persists. Neither
the 34/36 count nor the selected-pair magnitudes are confirmatory effect
tests. Residual RV can change in magnitude simply by changing marginal
variances; do not describe a larger post-adjustment RV as "climate makes
integration stronger."

**Unresolved:** no BIO18 or elevation in these inputs; no geographic/phylogenetic
joint correction; no image-viewing-angle/colour-error correction;
no inference that the remaining associations are adaptive strategies.
The three-PC dimensional choice itself needs a sensitivity comparison
(e.g. 2 and 4 PCs) before Main promotion.

Source: [Actions run 40055211966](https://github.com/zuizui0223/azami/actions/runs/40055211966)
and archived output artifact 12544872267.

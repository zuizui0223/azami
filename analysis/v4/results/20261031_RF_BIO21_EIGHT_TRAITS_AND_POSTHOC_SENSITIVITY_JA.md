# 12対21変数Random Forest：原S1全8構成形質と事後的サンプリング感度（2026-10-31）

**STATUS: source-real empirical, exploratory, original GEB science remains frozen.** [Successful protected run #41187277022](https://github.com/zuizui0223/azami/actions/runs/41187277022), aggregate artifact **13037606811**. This report DOES NOT replace the original positive Holm-tested regression results or prove a lack of ecological effects.

## Identified source and scope
- Original protected S1 313,350 photographs/observations; full 21-variable source complete cases **311,585**. 7 historically external CHELSA/WorldClim axes + 9 new exact-period CHELSA BIOCLIM layers were source-pixel sampled, each checked for source SHA, raster identity, period, units and missingness. New BIO9 complete at all 313,350 original points.
- Strict RF cohort: **one taxon × 5° cell** median from at least three eligible photo records per cell, taxa with ≥8 cells and ≥40 contributing photo records. **26 taxa** meet the final phenotypic/geographical gates per endpoint below (not the entire 178-taxon S1 angle cohort).
- Core comparisons: identical complete-case photo/cell cohort, train-mean constant baseline, inverse taxon-cell weighting; n_estimators=80, max_features=0.75, min_samples_leaf=5, 3 held-out permutation repeats, 5-fold unseen-taxon GroupKFold and 3 held-out longitude sectors.
- Original camera-relative head angle, JPEG colour, and projected outlines are NOT calibrated physical orientation/pigment/anatomy; taxon is not an RF predictor. H3c, S1/S2 Holm families and the existing manuscript are unchanged.

## Main exploratory strict 3-photo/cell result
Skill = 1 − held-out weighted model SSE / train-fold constant-mean SSE; negative values imply worse predictions than the training mean. **ΔSkill=Skill21−Skill12**. Not a new family of significance tests.

| Visible construct | Taxon-cells/taxa | Unseen-taxon 12 | Unseen-taxon 21 | Δtaxa | Unseen-sector 12 | Unseen-sector 21 | Δsector |
|---|---:|---:|---:|---:|---:|---:|---:|
| Image-frame angle | 647/26 | −0.104 | −0.120 | −0.017 | −0.058 | −0.039 | +0.019 |
| Lightness | 659/26 | −0.187 | −0.206 | −0.018 | −0.024 | −0.016 | +0.009 |
| Chroma | 659/26 | +0.038 | +0.040 | +0.002 | −0.268 | −0.205 | +0.063 |
| Joint circular hue (sin+cos) | 659/26 | −0.351 | −0.353 | −0.002 | −0.331 | −0.341 | −0.010 |
| Projected aspect | 653/26 | −0.258 | −0.277 | −0.019 | −0.115 | −0.060 | +0.055 |
| Circularity | 653/26 | −0.227 | −0.183 | +0.043 | −0.224 | −0.183 | +0.042 |
| Solidity | 653/26 | −0.195 | −0.208 | −0.012 | −0.126 | −0.165 | −0.039 |
| Width profile CV | 653/26 | −0.137 | −0.172 | −0.035 | −0.002 | −0.082 | −0.081 |

Strict cohort: **21-variable RF underperforms constant-mean predictions for 7/8 traits in unseen taxa and all 8/8 in unseen geographic sectors**; chroma unseen taxa is only slightly above zero (Skill +0.040). Predictive improvement from new 9 BIOCLIM variables is mixed: +2/8 unseen taxa, +5/8 unseen sectors, but absolute scores remain mostly negative. DO NOT call improved Δ a validated ecological effect when final Skill remains negative.

## Posthoc support and balanced-taxon sensitivity (separate estimand)

After observing the first 3-trait run, we explicitly changed minimum photo support per 5° cell from **3 to 1** (still ≥8 cells and ≥40 photographs per taxon), and replaced cell-balanced GroupKFold with **taxon-count-balanced unseen-taxon folds**. The original 3-photo/cell results above remain unchanged; this is *not* a preregistered confirmation.

| Construct | New taxon-cells/taxa | Taxa-balanced 12 | Taxa-balanced 21 | Δtaxa | Held-out-sector 12 | Held-out-sector 21 | Δsector |
|---|---:|---:|---:|---:|---:|---:|---:|
| Angle | 818/31 | −0.174 | −0.144 | +0.030 | −0.082 | −0.050 | +0.032 |
| Chroma | 849/34 | +0.028 | −0.011 | −0.038 | −0.053 | +0.059 | +0.112 |
| Aspect | 827/32 | −0.216 | −0.236 | −0.020 | −0.070 | −0.069 | +0.001 |

Strict GroupKFold had a fold containing only **1 held-out species** (139 angle, 141 chroma cells) despite species-equal weighting. Sensitivity now has 6–7 held-out species in each fold. After loosening photo support, taxon coverage increases to 31/34/32, still a selective subset of S1. **Chroma sector Skill +0.059 in the relaxed sensitivity** is not a reproducible positive result across these cohort designs (strict Skill −0.205). The 21-variable uplift for angle unseen taxa switches sign (−0.017 strict to +0.030 sensitivity), confirming selection/split dependence.

## Scientific decision
The RF does **not** presently support a general claim that 21 climatic variables predict these photographic traits across unfamiliar taxa or broad geographic sectors. The conditional geographical trait resemblance and specific original regression associations in the GEB manuscript ask distinct questions and are not invalidated by an RF with low out-of-domain skill. A more informative next test adds mechanistically distinct, independently measured climate fields or landscape variables—not further highly redundant BIOCLIM indices—under the same paired holdout logic.

## Frozen provenance, data privacy
- [Artifact](https://github.com/zuizui0223/azami/actions/runs/41187277022), ID **13037606811**; full ZIP SHA256 `8afc24640b652842ac08361c8a27fda43305cfa6e5e0f48d1ddc448b882b5907`.
- `rf_12_vs_21_full8_aggregate.json`: `13350d0ffc54ff2a09ba004e50eae5d498aaaf0bee115e7a56c2bbf8b7591072`.
- `rf_12_vs_21_posthoc_balanced_support_aggregate.json`: `3647c7e6452b0f1792d85c54bc2265c9b2027e836d578ef1c1ab7b2d96193f04`.
- `rf_12_vs_21_anchor3_aggregate.json`: `26081c5d30766ee3ff30cd74b4288f32b5606d0a23f120f378a622d53defa903`.
- `bio21_source_receipt.json`: `26b80b279fa2b372099ef945972d3e4fdf58127fe9f8fd92530ac08a931d9545`.
- Protected originals, photograph IDs, photographer IDs, locations, taxon identities and individual predictions never appear in this report or public Actions artifacts. These aggregate receipts cannot reconstruct confidential source rows.

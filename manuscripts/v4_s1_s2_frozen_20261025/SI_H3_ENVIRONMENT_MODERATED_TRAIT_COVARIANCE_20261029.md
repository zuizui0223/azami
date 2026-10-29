# Supporting Information: H3d — trait-specific climatic breadth and conditional covariance

## Question and design

H3a found no clear association of sampled **3D climatic volume** with **3D photographed capitulum morphospace volume**. H3d asks two separate biological questions: (i) whether that global volume metric obscures diversity differences in any **single visible trait**; (ii) whether **trait covariation itself** depends on precipitation or solar radiation.

The original S1 private source (313,350 observations, 157,822 matched photographs, 242 taxon keys) was decoded and analysed inside [Actions #40847696537](https://github.com/zuizui0223/azami/actions/runs/40847696537). The 1° analysis used **63** geographically adequately sampled taxa; the coarser 5° sensitivity used **31**. Anonymous [full aggregate](../../analysis/v4/results/20261029_H3D_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_SOURCE_VERIFIED.json) and [source code](../../analysis/v4/s1_h3_environment_moderated_trait_covariance.py). A [pre-analysis contract](../../analysis/v4/H3_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_DESIGN_20261029_JA.md) fixed the axes and Holm families.

**Axis-specific diversity:** same climate ellipsoid volume as H3a predicts 3 separate rarefied log standard deviations (image angle, colour chroma, outline aspect) with sampling/geographic controls. Three permutation comparisons, adjusted together by **Holm-3**.

**Environmental moderation of trait covariance:** remove each taxon's mean and its additive BIO12, BIO01 and RSDS main effects from all three cell-level shape measurements. Compare each pair's cross-product of the residuals under **higher vs lower than the same taxon's median BIO12** and independently **higher vs lower RSDS**. Three photographic trait pairs × two exposures = **six planned contrasts**, all corrected together by **Holm-6**. Each taxon needs ≥3 cells on each side; taxa are inferential units. Sign-flip 999 and taxon-level bootstrap 999 repetitions quantify observational variability, **not** experimental selection or genetic trait interaction.

## Empirical findings

### Single-trait diversity versus climatic breadth

| Morphological dimension | 1° slope (HC3 pointwise 95% CI) | Permutation P | Holm-3 P |
|---|---|---:|---:|
| Image angle | +0.0251 [−0.0506,+0.1009] | 0.4475 | 1.0000 |
| Corolla chroma | +0.0336 [−0.0675,+0.1348] | 0.4150 | 1.0000 |
| Outline aspect | −0.0124 [−0.1018,+0.0770] | 0.7625 | 1.0000 |

The 5° sensitivity (31 taxa) likewise had no Holm-3 supported individual trait response to climate breadth (P_adj 0.3825, 0.7900, 0.7900, respectively). **Non-detection does not establish equal trait breadths or morphological canalisation.**

### Trait covariance conditional on local climate

| Cross-domain pair | Climatic axis | 1° high–low residual covariance | Permutation P | Holm-6 P |
|---|---|---:|---:|---:|
| Angle × chroma | BIO12 annual rain | +0.0239 | 0.635 | 1.000 |
| Angle × outline aspect | BIO12 | +0.0828 | 0.389 | 1.000 |
| Chroma × outline aspect | BIO12 | +0.0710 | 0.172 | 0.860 |
| Angle × chroma | RSDS solar radiation | −0.0221 | 0.589 | 1.000 |
| Angle × outline aspect | RSDS | **+0.1858** | 0.016 | **0.096** |
| Chroma × outline aspect | RSDS | −0.0307 | 0.462 | 1.000 |

The **nominal** 1° RSDS × angle/outline covariance contrast does **not** survive the six-test Holm correction, and under 5° cell aggregation its point estimate **reverses sign** (−0.1897, uncorrected P=0.143). It must **not** be promoted as a positive ecological interaction. All six 5° Holm P values exceed 0.05.

**Conclusion:** conditional covariance of photographically derived head domains remains unconfirmed. Weak original pairwise correlations should not be rescued by describing them as climate-dependent functional integration. Compared with the source-confirmed **H3c climatic-analogue local morphological alignment**, the H3d interaction and niche-breadth response analyses are currently **negative/inconclusive exploratory tests**, not competing independent discoveries.

These values are not genetic covariances, joint fitness functions, independently evolving modules or a test of selection. Camera angle, lighting and sampling by observers can affect apparent shape and colour even after geographic-cell aggregation; climate and geography cannot be causally separated.

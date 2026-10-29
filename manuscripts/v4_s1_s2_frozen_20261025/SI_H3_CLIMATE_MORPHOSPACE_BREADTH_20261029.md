# Supporting Information — H3 climatic occupancy and multidimensional floral morphospace

## Motivation and source

The third ecological question in the V4 GEB Research Article is whether taxa occupying a wider **realised geographic climatic niche** also occupy a larger volume in a common **photographically observed floral morphospace**, and whether morphological similarity recurs between taxa from different geographic sectors when their climatic centroids are similar. These are **hypotheses about observed ecological organisation**, not genetic plasticity, fitness, phylogenetic convergence or adaptive radiation.

Original source: encrypted, SHA-verified V4 native source **313,350 observations** and CHELSA pixel-decoded climatology. In the original records, exactly **157,822** observations from **242** source taxon keys passed matched angle, colour-chroma and outline-aspect quality gates. Original observer identities and site coordinates were retained only in the Actions runner, not uploaded.

Reproducibility:
- Code: `analysis/v4/s1_climate_trait_hypervolume.py` and source-free unit tests `tests/test_v4_s1_climate_trait_hypervolume.py`.
- [Empirical source execution (success)](https://github.com/zuizui0223/azami/actions/runs/40841007115), [aggregate-only artifact 12883021095](https://github.com/zuizui0223/azami/actions/runs/40841007115/artifacts/12883021095).
- Permanent, identifier-free [machine-readable result](../../analysis/v4/results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json).
- Design was stated before viewing H3 results in [the independent H3 hypothesis and analysis contract](../../analysis/v4/HYPERVOLUME_CLIMATE_MORPHOSPACE_H3_PREDECLARED_20261029_JA.md).
- Original S1 80 and S2 128 Holm tests and the original colour/presentation coefficients were **not refitted or modified**.

## Operational definitions

Three climatic axes: CHELSA annual precipitation (BIO12), mean temperature (BIO01), mean shortwave radiation (RSDS). Three photographic axes chosen across biological domains: head image angle, CIELAB corolla chroma and two-dimensional outline aspect ratio. These are **putatively functionally distinct visible traits**, not directly measured functional performance; aspect and angle are camera-view dependent.

Within each taxon, each occupied 1° geographic cell was represented by the median of its quality-eligible photographs and the corresponding climatic values. Eligible taxa needed **at least 40 photographs and 8 cells**. A taxon-equal pooled standardisation fixed the axes; then exactly eight occupied cells were sampled without replacement **96 times per taxon** for both spaces. The volume metric is **1/2 log det(Covariance + 0.025 I)**, a regularised three-dimensional Gaussian-ellipsoid log-volume *proxy*. A convex hull, KDE Hutchinsonian hypervolume, climate physiological tolerance and morphological adaptability were **not** estimated. The matching 5°-cell analysis applied identical eligibility and estimation rules as a spatial-resolution sensitivity. All eight photographs/variable subsets were predeclared.

For H3a (within-taxon climate breadth versus phenotypic breadth), species-level regression adjusted for geographic spread, numbers of photographs and occupied cells, and fixed longitude sectors. The inference included an exploratory HC3 pointwise CI plus **399 sector-stratified Freedman–Lane residual permutations**, with the taxon as the unit.

For H3b (similar climate and similar morphology across regions), taxon climate-centroid distance was correlated with photographic morphological-centroid distance **only between different predeclared longitude sectors**. Source taxa, not individual pairwise distances, were permuted as whole three-dimensional trait vectors **within sector** 399 times, so the pairwise comparisons were not falsely counted as independent samples. Longitude sectors are **not independent radiations**; no source-matched large phylogeny was fitted.

## Original empirical findings

| Metric | Primary 1° cells | Secondary 5° cells |
|---|---:|---:|
| Taxa passing cell/photo gates | **63** | **31** |
| Climate log-volume → trait log-volume coefficient (HC3 95% CI) | **+0.0669** [−0.1046,+0.2384] | **+0.0700** [−0.4403,+0.5803] |
| Sector-stratified permutation P | **0.3875** | **0.7100** |
| Cross-sector taxon pairs | **1,098** | **295** |
| Climatic centroid distance × morphological centroid distance (Spearman ρ) | **−0.0081** | **−0.0954** |
| Whole-taxon within-sector permutation P | **0.9125** | **0.4250** |

On these photo-geographic samples and under the specified volume proxy, **no clear positive association** was detected between realised climatic breadth and photographic floral morphospace breadth. The cross-sector centroid distances likewise show no supported covariation in the direction expected for a simple shared-climate → similar-flower hypothesis. The 5° sensitivity had fewer qualifying taxa and broad uncertainty for H3a.

These findings cannot prove **equivalence, ecological decoupling, lack of convergent evolution or adaptation to niche breadth**. Insufficient within-species geographic coverage excluded many taxa. Realised niche breadth is limited by biased photographic occurrence sampling, and the photographic axes are incomplete proxies for performance. A source-matched nuclear phylogeny and replicated functional measurements are prerequisites for testing independently repeated morphological radiation or evolutionary convergence. The EAzami repository's verified nuclear 20-tip panel has substantial global sampling limitations and cannot be silently substituted for a genus-wide tree.

## Evidence hierarchy

The descriptive volume association is an **exploratory H3** analysis, with exactly two declared ecological comparisons (breadth relationship and cross-region centroid relationship). Its permutation probabilities are **not part of S1/S2's 208 original Holm comparisons** and must not be pooled with them as a larger confirmatory family. The original GEB 3.1 and 3.2 coefficients remain unchanged. See Blonder et al. (2014) for the conceptual distinction between a Hutchinsonian hypervolume and approximations to volume from multivariate samples: https://doi.org/10.1111/geb.12146.

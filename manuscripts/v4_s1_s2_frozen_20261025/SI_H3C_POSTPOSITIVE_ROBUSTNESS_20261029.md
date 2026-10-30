# Supporting Information — H3c post-positive sensitivity and random-sequence stability

**Evidence status: empirical, source-verified, exploratory.** [GitHub Actions #40859868003](https://github.com/zuizui0223/azami/actions/runs/40859868003) completed successfully. [Aggregate artifact #12891286463](https://github.com/zuizui0223/azami/actions/runs/40859868003/artifacts/12891286463) contains the full per-condition and six-seed result objects. The original decoded native S1 source is unchanged (SHA256 `342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed`, 313,350 total source observations, 157,822 photographs with the three matched visual traits); 63 taxa satisfied the primary 1° support gate and 31 the 5° gate. The tested model and original exploratory Holm-2 result are [preserved separately](SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md).

## A. Paired scenario comparison after the original positive H3c result

The original matching uses three climate dimensions (BIO12, BIO01, RSDS), three photographed head dimensions (image-referenced angle, CIELAB chroma, 2D outline aspect), climatic Euclidean z-distance ≤0.9 and fixed KDE width h=0.75. We examined alternate match thresholds (0.7 and 1.1), omission of each floral dimension, leave-one-longitude-sector-out, and KDE bandwidth (0.50 and 1.00). All sensitivity conditions used **24 independent eight-cell rarefaction draws and 99 taxon-preserving within-taxon phenotype-vector permutations**. **Where the taxon roster was the same, the random cell draws and permutation streams were identical across conditions**, so a KDE bandwidth-only change did not affect the three-trait distance estimates. This is a **post-positive sensitivity**, not a newly preregistered, independent family of confirmatory tests.

The following is the **observed distance improvement minus the corresponding permutation-null mean** in global standardised floral coordinates; positive means greater climate-conditioned morphological similarity relative to the null. P values are **unadjusted diagnostic permutation P** and are not additional independent evidence.

| Scenario | 1° effect | 1° nominal P | 5° effect | 5° nominal P |
|---|---:|---:|---:|---:|
| Same threshold 0.9; independent sampling stream | +0.00725 | 0.540 | +0.09828 | 0.010 |
| Narrower climate match 0.7 | +0.02734 | 0.010 | +0.11865 | 0.010 |
| Wider climate match 1.1 | +0.01112 | 0.270 | +0.07286 | 0.010 |
| Omit photo angle | +0.01009 | 0.330 | +0.07788 | 0.010 |
| Omit floral chroma | **−0.01158** | 0.330 | +0.07916 | 0.010 |
| Omit outline aspect | +0.01999 | 0.070 | +0.07476 | 0.010 |
| Omit western sector | +0.08182 | 0.010 | not examined | — |
| Omit central sector | **−0.00742** | 0.760 | not examined | — |
| Omit eastern sector | +0.03287 | 0.030 | not examined | — |

Changing the Gaussian-KDE bandwidth while keeping exactly the same climate/trait cell draws left the 1° **distance effect at +0.00725 (P=0.540)**. The **KDE density-overlap difference relative to its null** was +0.00544 (h=0.50; diagnostic P=0.020), +0.00471 (h=0.75; P=0.070), and +0.00403 (h=1.00; P=0.100). These density overlaps remain normalised L2 inner products, **not the physical intersection volume** of two bounded phenotypic hypervolumes.

**Interpretation of A:** The positive three-trait distance pattern is more stable in the geographically coarser, smaller 5° supported sample. At 1°, the finding is sensitive to the matching threshold, local trait composition (especially chroma) and which longitude sectors are compared. It should not be described as an invariant environmental convergence across all regional radiations.

## B. Original climate threshold and KDE width across six random-sequence replicates

We subsequently retained **exactly** the original z threshold 0.9, KDE h=0.75 and all three traits; six independent, predetermined random streams each ran 24 rarefaction draws and 99 within-taxon permutations. The six streams are different Monte Carlo realisations of **the same taxa and observations**, not independently replicated taxonomic studies.

| Monte Carlo stability | 1° (63 taxa) | 5° (31 taxa) |
|---|---:|---:|
| Streams with positive 3D distance effect relative to null | **6 of 6** | **6 of 6** |
| Median effect across streams | **+0.02123** | **+0.08388** |
| Range across streams | +0.01834 to +0.03275 | +0.06538 to +0.11148 |
| Streams with unadjusted nominal P below 0.05 | **2 of 6** | **6 of 6** |
| Range of nominal P | 0.010–0.160 | 0.010–0.010 |
| Streams with positive KDE overlap effect relative to null | 6 of 6 | 5 of 6 |
| Streams with nominal KDE P below 0.05 | 3 of 6 | 1 of 6 |

The original source-frozen H3c estimate remains **+0.04096, exploratory Holm-2 P=0.010 at 1°** and **+0.09433, exploratory Holm-2 P=0.010 at 5°**. The six-seed audit does **not** create six new independent positive findings. It demonstrates a **stable positive direction** of local climate-matched **distance** in both resolutions under these fixed seeds, but markedly **weaker stability of nominal statistical support at 1°** and weak seed stability of the density-overlap metric. This is not an independent out-of-sample validation.

## C. Boundaries and requirements for any stronger claim

Neither the 1° nor 5° analyses identify whether phenotypes have evolved convergently, whether selection favours similar forms, or whether camera orientation, light balance, habitat composition and shared lineage history account for some of the observed correspondence. Restricting to eight sufficiently represented cells per species also excludes taxa with limited photo-geographic sampling, particularly in the eastern longitude sector. The 5° sensitivity is a **different coarser aggregation with fewer qualifying taxa**, not a second independently sampled radiation.

The appropriate GEB conclusion is **an exploratory, spatial-resolution-dependent local photographic climate–phenotype association** with a reproducibly positive distance direction, **not** a globally uniform law or a functional trait interaction. The [parent H3c SI](SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md), [H3d trait covariance SI](SI_H3_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_20261029.md), [source code](../../analysis/v4/s1_h3c_seed_stability.py), and the archived full source-verified artifact provide the audit trail.

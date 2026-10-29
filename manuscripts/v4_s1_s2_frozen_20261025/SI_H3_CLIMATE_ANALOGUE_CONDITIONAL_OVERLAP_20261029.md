# Supporting Information: H3c — local climatic analogues and floral phenotypic alignment

## Question and empirical evidence

**Do taxa from geographically separated regions have more similar photographic capitula when their actual sampled climatic conditions are matched?** This is a test of **climate-conditioned phenotype association**, **not** a test of independently evolved similarity (phylogenetic convergence).

Code: `analysis/v4/s1_h3_climate_analogue_alignment.py`, exact original S1 decrypted only inside [Actions #40847696537](https://github.com/zuizui0223/azami/actions/runs/40847696537); archive `12885860903`; [anonymous machine-readable result](../../analysis/v4/results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json). The source-observation SHA256 is `342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed`, with 313,350 original native observations, 157,822 source-matched eligible photos/242 taxon names. A prior [predeclared H3c design](../../analysis/v4/H3_CLIMATE_ANALOG_PHENOTYPE_CONVERGENCE_DESIGN_20261029_JA.md) fixed the climate matching threshold and taxon-level null before looking at results.

Three environmental axes are BIO12, BIO01 and RSDS; three visible photographic floral dimensions are camera-referenced angle, CIELAB chroma and 2D outline aspect. We computed geographic-cell medians from *the same* source photo QC cohort as H3a/b. Eligible taxa required 40 photos and 8 cells; each was resampled to **8 equal-opportunity cells 24 times**. A climatic analogue between observations from *different longitude sectors and different taxa* was predeclared as multivariate climatic distance ≤**0.9** on the same global z scales, with at least 3 matched cross-cell pairs and at least two source cells on each side.

### Two independently interpretable metrics, both within a single exploratory Holm-2 family

1. **Three-trait photographic distance matching advantage**: for each qualifying cross-sector taxon pair, subtract the mean 3D floral distance of climate-analogue cell pairs from the distance of all 8×8 cell combinations. An advantage above the null suggests local climate-conditional alignment of photographed traits.
2. **3D floral Gaussian-KDE density similarity/overlap**: calculate the normalised L2 inner product between fixed-bandwidth (h = 0.75 global z units) Gaussian KDEs for each pair of taxa, using either **all 8 cells** or **the cells participating in climatic analogues**. This is a bounded 0–1 density similarity and **is not the intersection/Jaccard area of thresholded Hutchinsonian hypervolumes**. Compare the matched-minus-all difference against the same within-taxon permutation null. A negative absolute difference may still be **less negative than the null**; do not turn that into a claim that conditioning increases overlap over the unconditional set.

The null permutes full 3D image phenotype vectors **among cells within each taxon** (199 permutations, identical taxon permutation across every pair in a draw). The source climate, taxon mean and overall phenotype distribution, and regional taxon assignments are retained. The 24 resamplings are averaged *before* the p calculation; they are not 24 independent experiments. Holm adjusts the two predeclared matching metrics, not S1's 80 or S2's 128 original tests.

## Source-real results

| Metric | 1° primary | 5° sensitivity |
|---|---:|---:|
| Geographically eligible native taxa | **63** | **31** |
| Dominant-sector composition (west / middle / east) | 36 / 21 / 6 | 15 / 11 / 5 |
| Qualified cross-sector taxon pairs per rarefaction, mean | **229.2** | **103.3** |
| Morphological-distance advantage (unconditional minus climate analogue) | **+0.04096** | **+0.09433** |
| Mean null distance advantage | −0.00127 | +0.00018 |
| Distance: permutation P / Holm-2 P | **0.005 / 0.010** | **0.005 / 0.010** |
| KDE L2 overlap: full 8-cell / climate-analogue subset | 0.61509 / 0.50742 | 0.53554 / 0.46046 |
| Conditional − unconditional KDE overlap | −0.10767 | −0.07508 |
| Mean null conditional − unconditional overlap | −0.11423 | −0.08349 |
| KDE relative-to-null: permutation P / Holm-2 P | **0.010 / 0.010** | **0.015 / 0.015** |

**Positive finding:** cross-region photographs at matched local climate show closer **three-dimensional visible head traits than expected under permuting the climate/photo-trait relation within taxa**, even though the previous species-climate-centroid correlation (H3b) was not supported. The KDE finding is **relative to its permutation null**, not a greater overlap than all cells.

**Current limits:** observational source photos are not gravity-calibrated or illumination-standardised; 8-cell distributions yield bandwidth-sensitive KDE estimates; climatic analogues are based on sampled locations rather than available backgrounds. East-sector taxa are only **6** in the primary panel and there is no globally comprehensive matched nuclear phylogeny. The result may reflect local within-taxon change, geography or photographic conditions. It **does not establish adaptive filtering, plasticity, fitness, evolutionary convergence or independent origins**. Additional region-pair and bandwidth sensitivity analyses will have to be considered exploratory rather than retroactively included in the original Holm-2 family.

For KDE overlap versus set-intersection hypervolumes see Blonder et al. (2014), *Global Ecology and Biogeography*, https://doi.org/10.1111/geb.12146.

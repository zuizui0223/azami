# GEB v2 within- versus among-taxon environmental structure — 2026-08-27

## Provenance and frozen design

This analysis extends the completed GEB v2 artifact without remeasuring images.

- source GEB v2 artifact ID: `9612943217`
- source GEB v2 run: `32975451732`
- source artifact digest: `sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`
- analysis branch: `analysis/geb-v2-within-among-alignment`
- pull request: `#72`
- within/among contract commit before outcome inspection: `16c329701fdf0c1470e112d067e7927461f861ad`
- expanded niche-sorting contract commit before outcome inspection: `250c35d42ffde6a470faf4babb8effa7916a1fe0`

The cross-scale analysis uses the same four CHELSA predictors as the GEB v2 within-taxon screen. Linear endpoints are summarized by taxon medians and tested with standardized taxon-level slopes; circular hue remains one joint sine/cosine inferential unit. The primary among-taxon scope requires at least five usable observations for that endpoint in a taxon and at least 20 taxa per model. A minimum-two-observation scope is a frozen sensitivity. P values use 10,000 taxon-label permutations and BH correction is separate for primary and candidate tiers.

Taxon-label permutation is not phylogenetic correction. Cross-scale agreement or disagreement is a spatial-pattern result and does not demonstrate plasticity, adaptation, selection or mechanism.

## Predictor-aligned within/among comparison

The measured GEB v2 universe contains 17 inferential units after the two hue components are treated jointly: eight primary units and nine candidate units. Across four CHELSA predictors this gives 68 among-taxon endpoint-predictor models in the primary replication scope.

Four of 68 among-taxon rows passed the frozen FDR rule. All four belonged to the primary tier; no candidate endpoint passed among-taxon FDR.

| Inferential unit | Predictor | Among-taxon effect | Among q | Within-taxon status | Cross-scale interpretation |
| --- | --- | ---: | ---: | --- | --- |
| `orientation_image_vertical_angle` | BIO12 annual precipitation | beta = +0.304359 | 0.006399 | beta = +0.005326, q = 0.874244 | among-only |
| `corolla_hue` | BIO1 annual mean temperature | joint magnitude = 0.303788 | 0.010399 | joint q = 0.015993 | both scales; colour calibration still open |
| `corolla_hue` | BIO4 temperature seasonality | joint magnitude = 0.334766 | 0.006399 | joint q = 0.010667 | both scales in the >=5 scope, threshold-sensitive among taxa |
| `corolla_hue` | BIO15 precipitation seasonality | joint magnitude = 0.410191 | 0.003200 | joint q = 0.010667 | both scales; colour calibration still open |

Using the primary >=5 among-taxon scope, the full endpoint-predictor alignment classified 3 rows as `both_scales`, 8 as `within_only`, 1 as `among_only` and 56 as `neither`.

The >=2-observation sensitivity again produced four among-taxon FDR rows. Three were shared with the primary scope: orientation-BIO12, hue-BIO1 and hue-BIO15. Hue-BIO4 lost among-taxon FDR support, while hue-BIO12 gained it. Thus the stable among-taxon pattern across the two frozen replication thresholds is concentrated in orientation-BIO12 plus joint hue associations with BIO1 and BIO15.

## The expanded involucral signals are within-taxon, not detected among taxa

The three candidate rows that entered the final GEB v2 evidence atlas did not show FDR-supported among-taxon sorting in the predictor-aligned analysis.

| Candidate endpoint | Predictor | Within-taxon GEB v2 result | Among beta, >=5 | Among q, >=5 | Among q, >=2 | Cross-scale result |
| --- | --- | --- | ---: | ---: | ---: | --- |
| `involucre_length_width_ratio` | BIO15 | +0.054923, q = 0.000282; quality-robust C grade | -0.010161 | 0.985601 | 0.733967 | within-only |
| `involucre_apical_taper_ratio` | BIO12 | -0.049721, q = 0.000821; quality-robust C grade | +0.004839 | 0.985601 | 0.467053 | within-only |
| `bract_projection_p95` | BIO15 | -0.033673, q = 0.014377; image-sensitive C grade | +0.090504 | 0.897562 | 0.733967 | within-only |

Because the among-taxon estimates are unsupported, their signs should not be biologically interpreted. The important result is the absence of detectable taxon-level sorting for these candidate associations under either frozen replication threshold. At present, the finer involucral signals are therefore better described as within-taxon geographic phenotype-environment structure rather than as differences among taxa occupying different climates.

## Expanded environmental-PCA sorting across the complete measured trait universe

A second frozen analysis extends the legacy low/high-trait environmental-PCA permutation test to the full measured endpoint set while keeping all linear endpoints on the same taxon set. Circular hue is tested jointly and is never converted to ordered colour quartiles.

The primary common-taxon scope required at least five usable measurements for every one of the 18 measured inferential endpoints. This retained 44 taxa. Across 16 linear endpoints, no endpoint passed BH correction for centroid distance or Bhattacharyya overlap, either jointly or singly. Joint circular hue retained an environmental association (`R2 = 0.175622`, permutation `p = 0.015198`).

The >=2-observation common-taxon sensitivity retained 78 taxa. Orientation was the only linear endpoint passing both environmental-sorting metrics: centroid distance = 1.294203 (`q = 0.046195`) and Bhattacharyya overlap = 0.627559 (`q = 0.031497`). Joint hue was also supported (`R2 = 0.145357`, permutation `p = 0.000800`). No candidate linear endpoint passed the expanded niche-sorting FDR rule.

## Cross-scale biological synthesis

The new symmetric analysis strengthens the GEB v2 conclusion that environmental structure is trait- and scale-specific rather than one integrated climate-associated capitulum syndrome.

First, the strongest expanded involucral associations are currently a within-taxon phenomenon. They remain useful exploratory spatial signals, but there is no evidence in this analysis that taxa with different median involucre length/width, apical taper or p95 projection are systematically sorted among the same CHELSA gradients.

Second, orientation has different environmental signatures at different biological scales. The established within-taxon result is orientation-BIO1, whereas the strongest replicated among-taxon coefficient is orientation-BIO12. The latter is consistent with the earlier taxon-level historical-sensitivity result for orientation and annual precipitation, but the present taxon-label test is not a phylogenetic correction.

Third, visible hue contains the clearest cross-scale environmental structure, but it remains constrained by camera-colour calibration and must be interpreted jointly. BIO1 and BIO15 are the most replication-threshold-stable among-taxon hue associations.

The resulting conceptual model is therefore:

> continuous capitulum phenotype -> distinct within-taxon geography + distinct among-taxon sorting -> partial cross-scale overlap, not one universal environmental syndrome

This is stronger than treating either the within-taxon or among-taxon layer alone as the biological conclusion.

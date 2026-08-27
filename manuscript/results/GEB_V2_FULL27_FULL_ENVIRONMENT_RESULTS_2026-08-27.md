# GEB v2 full-27 / full-environment Chapter 1 result ledger

## Decision

This is the canonical result ledger for the new Chapter 1 spatial-scale lane. It starts from all 27 endpoints registered in the GEB v2 continuous-trait contract and all nine frozen environmental predictors. The older nine-endpoint result set was not used to select traits, and analysis-tier labels were metadata only. No A/B/C layer removed rows or defined a multiplicity family.

The analysis is retrospective exploratory because related GEB v2 outcomes had already been inspected before the full 27 x 9 atlas was declared. It therefore establishes a reproducible present-geography result and candidate hierarchy, not confirmatory proof of adaptation.

## Canonical flow

> photograph -> 27 registered continuous endpoints -> variation hidden by taxon means -> endpoint-specific geographic support -> six environmental-hypothesis blocks -> within/among scale comparison -> sampling-composition sensitivity + broad-space sensitivity -> historical-placement sensitivity -> adaptive-pattern candidates under current controls

## Cohort and measurement status

- Frozen spatial cohort: 46,276 coordinate-bearing observations from 259 source-assigned taxa.
- Registered endpoints at entry: 27.
- Endpoints with finite measurements: 22.
- Joint inferential units: 26, because hue sine and cosine were tested together as one circular unit.
- Five endpoints remained explicit `unexecuted_no_measurement` rows: `visible_floret_fraction`, `corolla_white_pixel_fraction`, `corolla_redmagenta_pixel_fraction`, `corolla_purple_pixel_fraction` and `corolla_yellow_pixel_fraction`.
- Measurement coverage ranged from 909 observations and 85 taxa for validation-only surface proxies to 40,510 observations and 253 taxa for Lab lightness and chroma.
- Measured endpoints occupied 190-431 five-degree cells. The large primary colour/orientation/outline layer extended from approximately 55.1 degrees S to 71.0 degrees N and crossed the dateline. The expanded involucral layer occupied 259 cells, and the surface-proxy layer occupied 190 cells.

Missing measurement was never encoded as biological absence or zero. The 27-endpoint starting universe therefore became 22 measured endpoints plus five visible unexecuted rows, not a selected nine-endpoint universe.

## What taxon means lose

Across the 22 measured endpoints, the fraction of observation-level sum of squares below source-assigned taxon means ranged from 0.8137 to 0.9863. The smallest value was for Lab lightness and the largest for capitulum outline aspect ratio. Thus, in the observed, strongly unbalanced public-image sample, a one-row-per-taxon mean compressed 81.4-98.6% of visible variation.

Because abundant taxa can dominate that raw partition, an equal-replication sensitivity sampled two observations per eligible taxon. Its endpoint-specific median below-taxon fraction ranged from 0.2869 for hue sine to 0.5846 for basal taper ratio. Even under equal replication, taxon means therefore compressed a substantial 28.7-58.5% of visible variation.

These are visible image-phenotype partitions, not heritability or genetic-variance estimates. Their methodological implication is narrower and direct: repeated, coordinate-bearing photographs retain within-taxon phenotype information that categorical and taxon-mean databases cannot represent.

## Six frozen environmental-hypothesis blocks

| Block | Predictors | Interpretive boundary |
|---|---|---|
| Thermal | BIO1, BIO4 | annual mean temperature and temperature seasonality |
| Hydric | BIO12, BIO15 | annual precipitation and precipitation seasonality |
| Radiative / atmospheric | rsds, VPD | shortwave radiation and atmospheric drying power |
| Mechanical | sfcWind | mean near-surface wind exposure |
| Growing-season water input | GSP | climatic growing-season precipitation, not taxon-specific flowering-season rain |
| Resource / productivity | NPP | modelled potential productivity, not local measured resource supply |

Coverage was 100% for eight predictors and 99.989% for wind; all exceeded the frozen 98% gate. All 1,874 rows overlapping the earlier complete-18 process artifact matched exactly for taxon identity and all nine environmental values (maximum absolute difference 0).

## Endpoint-level atlas before spatial control

BH correction was applied once across every successful inferential unit x nine-predictor row, separately for the within-taxon, among-taxon min5 and among-taxon min2 scopes. Unsupported rows were retained.

The within-taxon atlas contained 189 executed tests and 45 unexecuted rows; 26 executed rows had global q < 0.05. The signals were selective rather than a whole-capitulum syndrome:

- Thermal: orientation increased with BIO1; outline aspect ratio increased with BIO4; joint hue varied with BIO1 and BIO4; surface edge density and high-frequency energy increased with BIO1.
- Hydric: joint hue varied with BIO12 and BIO15; apical taper decreased with BIO12, projection p95 decreased with BIO15 and involucre length/width increased with BIO15.
- Radiative / atmospheric: joint hue and projection roughness varied with rsds; chroma and projection-peak density decreased with VPD, while surface high-frequency energy increased with VPD.
- Mechanical: joint hue varied with wind.
- Growing-season water input: projection maximum, projection asymmetry, involucre length/width, surface edge density and specular fraction decreased with GSP; surface high-frequency energy and LBP entropy increased with GSP.
- Resource / productivity: joint hue varied with NPP and basal taper ratio decreased with NPP.

The primary among-taxon min5 atlas contained 180 successful tests, 45 unexecuted rows and nine no-finite-variation rows for the surface specular proxy. Ten rows had global q < 0.05:

| Endpoint / joint unit | Predictor | Block | Standardized effect | Global q |
|---|---|---|---:|---:|
| orientation | BIO12 | hydric | +0.3044 | 0.0210 |
| orientation | GSP | growing-season water input | +0.2581 | 0.0338 |
| chroma | rsds | radiative / atmospheric | -0.3454 | 0.0060 |
| bract projection-peak density | VPD | radiative / atmospheric | +0.4437 | 0.0468 |
| joint hue | BIO1 | thermal | magnitude 0.3038 | 0.0334 |
| joint hue | BIO4 | thermal | magnitude 0.3348 | 0.0210 |
| joint hue | BIO15 | hydric | magnitude 0.4102 | 0.0060 |
| joint hue | rsds | radiative / atmospheric | magnitude 0.5004 | 0.0060 |
| joint hue | wind | mechanical | magnitude 0.2974 | 0.0340 |
| joint hue | NPP | resource / productivity | magnitude 0.3155 | 0.0210 |

The min2 replication scope supported seven rows. Orientation-BIO12, orientation-GSP, chroma-rsds and joint hue with BIO1, BIO15, rsds and NPP repeated there; joint hue-BIO4, joint hue-wind and projection-peak-density-VPD did not pass its separate global family.

## Within- versus among-taxon scale classes

Using min5 as the primary among-taxon scope, the complete 234-row grid contained:

- 7 supported at both scales;
- 18 within-taxon only;
- 3 among-taxon only;
- 152 supported at neither scale;
- 54 not comparable because measurement was unexecuted or the min5 response had no finite variation.

The seven both-scale rows were joint hue with BIO1, BIO4, BIO15, rsds, wind and NPP, plus projection-peak-density with VPD. The three among-only rows were orientation-BIO12, orientation-GSP and chroma-rsds. Most supported involucral and surface-proxy rows were within-only.

This rejects scale invariance. A within-taxon gradient cannot be read automatically as the process that produced taxon turnover, and an among-taxon gradient cannot be read backward as plasticity.

## Sampling-composition sensitivity

The 46,276-observation cohort was geographically and taxonomically extensive but strongly unbalanced. Europe contributed 23,337 observations (50.43%) and North America 19,357 (41.83%); together they contributed 92.26%. The Northern Hemisphere contributed 43,920 observations (94.91%). *Cirsium vulgare* contributed 13,207 observations (28.54%), *C. arvense* 11,795 (25.49%), and the top ten taxa 36,990 (79.93%). The taxon-level median was six observations; 51 of 259 source-assigned taxa had one observation and 150 had fewer than ten.

An independent contract declared a retrospective selected-row sensitivity without changing the main atlas contract or its hash. Every global-q-supported row entered: 26 within taxa and ten among taxa. Within rows were checked under equal total taxon weight; both scales omitted each top-ten taxon separately and the top two jointly, omitted each of six Natural Earth broad regions separately, and retained explicitly native observations only. This produced 674 scenario rows, all evaluable. No post-selection P value or q value was calculated.

All ten among-taxon pairs retained direction in every applicable scenario. Sixteen of 26 within-taxon pairs did so; ten reversed in one or more scenarios. Five reversals appeared after equal taxon weighting, three after dominant-taxon omission, three after removal of North America and two in native-only data, with overlap among pairs. The direction-unstable set was:

- projection p95-BIO15;
- projection roughness-rsds;
- joint hue with BIO12, BIO15 or rsds;
- apical taper-BIO12;
- basal taper-NPP;
- surface edge density-GSP;
- surface LBP entropy-GSP;
- surface specular fraction-GSP.

These rows remain in the complete atlas with a sampling-composition-sensitive annotation. Directional stability is not significance preservation, a correction to probability sampling or evidence of causation.

## Broad-space sensitivity

Only global-q-supported rows entered this gate: 26 within-taxon and ten primary among-taxon rows. Each was refitted with a second-order spherical-coordinate basis. Passage required a spatial permutation P < 0.05, preserved sign for a linear endpoint and residual eight-nearest-neighbour Moran P >= 0.05.

Five within-taxon rows passed:

| Endpoint | Predictor | Spatial beta | Spatial P | Residual Moran P |
|---|---|---:|---:|---:|
| orientation | BIO1 | +0.0265 | 0.002 | 0.324 |
| chroma | VPD | -0.0421 | 0.001 | 0.381 |
| involucre length/width | BIO15 | +0.0602 | 0.022 | 0.411 |
| projection roughness | rsds | -0.1392 | 0.005 | 0.626 |
| surface high-frequency energy | BIO1 | +0.1235 | 0.026 | 0.288 |

Two among-taxon rows passed:

| Endpoint | Predictor | Spatial beta | Spatial P | Residual Moran P |
|---|---|---:|---:|---:|
| chroma | rsds | -0.7124 | 0.001 | 0.217 |
| orientation | BIO12 | +0.2861 | 0.017 | 0.082 |

All joint-hue rows failed the full gate because residual spatial autocorrelation remained detectable, even where the spatial-basis effect itself retained P < 0.05. They remain geographically structured colour patterns, not current adaptive-pattern candidates. The broad spherical screen is not equivalent to the older SPDE analysis and does not remove every possible spatial process.

Four of the five within-taxon broad-space passes retained direction in every sampling-composition scenario. Projection roughness-rsds reversed after equal taxon weighting and after jointly omitting *C. vulgare* and *C. arvense*. It therefore remains a broad-space pass but not a sampling-composition-stable result.

## Historical-placement sensitivity and final candidates

The two spatially retained among-taxon pairs were evaluated on 52 audited dated-backbone placement trees: deterministic scenarios S1 and S3 plus 50 randomized S2 trees. Both pairs were successful, retained direction and had P < 0.05 on all 52 trees.

| Candidate | Among beta | Broad-space beta | PGLS beta range | PGLS P range | Pagel lambda range | Final status |
|---|---:|---:|---:|---:|---:|---|
| lower chroma with higher rsds | -0.3454 | -0.7124 | -0.3454 | 0.000024 | 0 | adaptive-pattern candidate under current controls |
| larger signed head-axis angle relative to EXIF-image vertical with higher BIO12 | +0.3044 | +0.2861 | +0.3044 to +0.3045 | 0.000201-0.000231 | 0-0.0537 | adaptive-pattern candidate under current controls |

These are sequential spatial and historical-placement sensitivities, not one resolved joint spatial-phylogenetic causal model. Only 54 of the 216 historical-atlas taxa were direct backbone tips; most were genus-grafted. Passing therefore means stability under the tested broad-space and placement uncertainties. It does not demonstrate heritability, selection, fitness advantage, functional mechanism, independent origin or convergence.

Both candidates also retained direction in all 18 applicable sampling-composition scenarios. The smallest absolute effect-magnitude ratio was 0.615 for chroma-rsds and 0.781 for orientation-BIO12, in each case under broad-region omission. This is a retrospective robustness annotation and does not change the candidate claim ceiling.

Five within-taxon rows passed the broad-space screen, and four of those were also directionally stable throughout the sampling-composition audit. They do not enter the phylogenetic candidate gate and are not evidence of plasticity or local adaptation.

## Chapter 1 conclusion

Coordinate-bearing photographs recover continuous, repeated capitulum phenotypes that categorical and taxon-mean databases compress. The recovered variation is not organized by one head syndrome or one environmental axis. It is endpoint specific and scale dependent: many involucral relationships occur only within taxa, ten selected within-taxon rows are sampling-composition sensitive, several colour relationships appear at both raw scales but remain spatially structured, and only chroma-radiation and orientation-annual-precipitation survive the current sequential among-taxon spatial and historical-placement gates while also retaining direction under the new sampling perturbations.

The methodological contribution is therefore inseparable from the biological result: retaining continuous observation-level values and coordinates reveals both the variation hidden by taxon means and the scale at which each environmental alignment occurs.

## Provenance and validation

- GEB v2 continuous-trait universe: GitHub Actions run `32975451732`, artifact `9612943217`, artifact digest `sha256:101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e`.
- Within/among source artifact: run `33035785120`, artifact `9632715852`, digest `sha256:51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6`.
- Historical placement resource: run `29091625109`, artifact `8227254443`; tree audit recorded 54 direct tips, 162 within-genus grafts and 50 randomized S2 trees.
- Trait-long SHA-256: `d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344`.
- Full nine-variable environment SHA-256: `e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7`.
- Main atlas validation: 24/24 checks passed.
- Spatial/historical validation: 11/11 checks passed; 104/104 historical models completed.
- Sampling-composition validation: 7/7 checks passed; 674/674 declared scenario rows evaluated.
- Frozen result tables: `analysis_outputs/v2_full27_environment_atlas_2026-08-27/`.

## Claim ceiling

This ledger supports an exploratory, endpoint-specific present trait-environment atlas; a quantified taxon-mean information-loss result; explicit scale classes; five within-taxon broad-space passes of which four are also directionally sampling-composition stable; and two adaptive-pattern candidates under current sequential controls that additionally retain direction in the retrospective sampling audit. It does not support a globally representative probability sample, demonstrated plasticity, local adaptation, causal environmental mechanisms, botanical identity of validation-only image proxies, resolved phylogenetic correction, independent origins, functional convergence or adaptive convergence.

# V4 S1–S2: expanded BIOCLIM, elevation and capitulum LEGO — real source result (2026-10-24)

## Verified execution / original units

- **Status:** GitHub Actions run [40057456480](https://github.com/zuizui0223/azami/actions/runs/40057456480) completed successfully; job 120763479106, artifact **12545642886** (`v4-s1-s2-extreme-lego-40057456480`), ZIP SHA256 `eee2c03a8d090db1d228dffb2ca23b70d99eeabe40db1b867870a2e93d3dcc5f`.
- Frozen original S1: **313,350 native observation IDs**. S2: **297,622 source-authenticated wild native or explicitly introduced observations** of *C. arvense* and *C. vulgare*.
- Original source positions combined without double-counting are **469,958 distinct observations**.
- Six actual CHELSA V2.1 BIO05/BIO06/BIO13/BIO14/BIO17/BIO18 rasters sampled at original pixels: **469,958/469,958** finite each. Separate WorldClim V2.1 elevation (m) **466,546/469,958** finite (99.27%). Sources, original raster SHA and decoding in `environment_source_receipt.json`, SHA `7f16950d8ca2b355c901ead965da152d22296e7b1ad15515cb99496649917d05`.
- The complete aggregate result JSON SHA is `f1fc39fa829318ce7289b27857255992f00b941094faf16d3b43812a77c04209` (`S1_S2_extreme_LEGO_results.json`).
- Both source restoration and numerical tests passed. No original observation IDs, uploader IDs or coordinates were added to this report / public artifacts.

## 1. The requested multicollinearity diagnostic comes first

All **12 predictors** = original BIO01/04/12/15/RSDS + BIO05/06/13/14/17/18 + elevation.

| Scope (predictors only, before trait filtering) | n | high-correlation example | Pearson r | max VIF (12-variable hypothetical joint model) |
|---|---:|---|---:|---:|
| raw original observations | 311,585 | BIO14—BIO17 | +0.991 | 104.98 (BIO06) |
| equal-taxon climate means | 303 taxa | BIO14—BIO17 | +0.997 | 218.56 (BIO17) |
| within-taxon centred | 311,585 | BIO14—BIO17 | +0.984 | 68.58 (BIO06) |

Additional contrasts: BIO12—BIO13 r +0.923 raw / +0.935 among taxa; BIO01—BIO06 r +0.903 raw; BIO12—BIO18 r **+0.798 raw, +0.896 among, +0.672 within**. Therefore BIO12 and BIO18 are conceptually different, but **not independent** on current species-level climate coverage. Avoid simultaneous arbitrary insertion of all 12.

Outcome-blind VIF + |r| screening produced nine alternative supported scenario designs. In the annual-rain and warmest-quarter scenarios, the optional elevation term was **dropped**, so the with-/without-elevation versions collapse to the same selected model. This is not evidence of zero elevation biology. The cold-BIO06/warm-rain scenario retained elevation and dropped BIO04 instead. No environmental model was selected by looking at a trait p-value.

## 2. V4 S1 phenotype LEGO: eight features, not V3 nine constructs

Complete linked S1 **eight-observable-feature × twelve-exposure** cohort: **156,683 original observations / 174 taxa** (before 12-exposure restriction: 157,654 / 174). V4 core = photo-relative presentation angle; lightness, chroma and joint circular hue; aspect ratio, circularity, solidity and width-profile CV. Exactly **28 unordered pairs** on identical observations/taxa, grouped into presentation / colour / outline modules.

- **Raw RV:** **21/28** pairs have stronger among-taxon than within-taxon RV. This is descriptive raw RV, **not V3 33/36**.
- Colour module has strong internal coupling (e.g. lightness×hue within 0.1042 vs among 0.4273; chroma×hue 0.0601 vs 0.3482).
- Cross-module relationship angle×aspect has within 0.0048 vs among 0.0722; not a genetically distinct orientation module.
- Circularity×solidity is an exception: within 0.734 vs among 0.690.
- **Nine same-module pairs:** median RV within 0.0850, among 0.2238; **19 cross-module pairs:** median within 0.00140, among 0.01062. This is consistent with phenotypically modular *appearance*, NOT verified independently inherited components.

| Identical S1 cohort, exploratory conditioning | among-stronger raw RV pairs | matched-n median among-minus-within RV | 95% matched-n resampling range | matched-n median among-stronger pairs |
|---|---:|---:|---|---:|
| before new climate adjustment | 21/28 | +0.01404 | [+0.00383, +0.01900] | 19/28 |
| adjusted BIO12 annual climate model | 23/28 | +0.01023 | [+0.00030, +0.01538] | 19/28 |
| adjusted BIO18 warmest-quarter model | 23/28 | +0.01246 | [+0.00251, +0.01769] | 19/28 |
| adjusted BIO17 driest-quarter candidate | 24/28 | +0.01173 | [+0.00179, +0.01685] | 20/28 |

These 500 matched-n resamples retain *fixed observed taxa and among-taxon means* and are **Monte Carlo sensitivity distributions, not independent biological replications or population confidence intervals**. Do not use raw pair counts as significance tests or claim that these data prove an evolutionary increase in integration.

## 3. Beyond pairs: aza3-inspired three-part phenotypic configuration

The first three-way test is photo-relative angle × image chroma × image aspect. With taxon-specific median splits, **124/135** taxa with sufficient distinct observations have all eight recorded three-variable combinations (135 eligible taxa, 156,384 observations). A taxon-balanced three-way residual parity contrast relative to the pairwise-only maximum-entropy reference was +0.03849, taxon-bootstrap range [+0.01424,+0.06311].

**Crucial second coding sensitivity:** under source taxon-centred continuous feature coding and subsequent binary sign conversion (another legitimate unselected threshold), the triple co-positive proportion was **0.08696 observed vs 0.08995 pairwise-only max-entropy**, difference **−0.00299**, taxon-bootstrap 95% range **[−0.00663,+0.00078]**. Among-taxon third-order difference was essentially zero. Thus **nonzero higher-order excess is NOT robust to the phenotypic state threshold**. Do not promote a unique high-order LEGO interaction or all-eight-configurations prevalence as genetic recombination.

A separate *candidate* protective-envelope proxy was actually measurable in **942 source-observations / 37 taxa** after 12-exposure QC (1,075 before all restrictions); colour-chroma × outward bract-projection geometry raw among RV 0.1915 vs within RV 0.00119, but the projection is **neither true spine length nor stickiness** and this cohort is not the main 156,683 cohort.

## 4. V4 S2 native-to-introduced original climate and observable shape/colour

WCVP original wild source, same geographic-source-denominator per taxon:

| Species | native source n | explicitly introduced source n | introduced minus native BIO12 median (mm/year) | introduced minus native BIO18 median (warmest-quarter mean monthly mm) |
|---|---:|---:|---:|---:|
| *C. arvense* | 78,646 | 70,964 | +263.0 | +48.2 |
| *C. vulgare* | 62,368 | 85,644 | +298.0 | +49.7 |

The original S2 archive contains **actual linked observation-level angle, Lab colour and 2D head-shape variables** in both ranges; all eight V4 constructs can be audited rather than imputed. S2 models were fit on *separate original photo-matched* subsets within each species×range (date/photographer/5° cell balancing and cell-robust errors), comparing annual BIO12 vs BIO18 on identical model observations. S2 does **not** establish independent invasion replicate lineages, genetic/temporal change, causal rainfall effects or adaptation.

Descriptive candidate examples, **uncorrected post-hoc tests only**: native *C. vulgare* chroma vs BIO18 +1.743 image-chroma units / within-range climate SD (cell-cluster uncorrected P=.005), but the native/introduced gradients differ and photographer/pose effects remain unresolved. Its native head-aspect vs BIO18 slope +0.007 image ratio units / climate SD (P=.018 uncorrected), NOT a confirmatory head-shape ecological discovery. Some S2 with-elevation fits were stopped due to VIF rather than silently inferred.

## Safe positive conclusion and stop conditions

**Supported observational pattern:** The V4 cohort contains nontrivial and scale-dependent **image-visible floral colour / capitulum presentation / 2D shape configurations**, strongest within-modules, with several cross-module pair contrasts. The qualitative tendency for stronger among-taxon RV survives alternative original-pixel hydric adjustment under matched-n Monte Carlo checks. Main-vs-SI eligibility, measurement reliability, geographical/phylogenetic residuals and multiple-testing families remain to be reviewed before elevating new post-hoc climate regressions.

**Not established:** all phenotypes from the exact same physical head (some original photo/observation values are head medians), inherited independent genetic LEGO modules G, higher-order functional resource retention G123, trait reassembly through recombination, independent evolutionary origin, adaptive rain protection, pollinator or enemy fitness advantage, or lineage history.

To prove aza3 *genetic/functionally reconfigurable* LEGO specifically: same individual/head linked floret colour chemistry and regulatory/catalytic G; actual gravity-referenced angle & development; physical phyllary spines, secretion and their genetic loci; ancestry/cytotype; functional insect access & viable achenes. Independent genealogy/segregation is indispensable.

**Provenance:** V4 PR #146 (base V4 PR #144), source-run 40057456480, numerical aggregate ZIP artifact ID 12545642886. The V3 33/36 figure is not reused and the V4 source parent is not merged.

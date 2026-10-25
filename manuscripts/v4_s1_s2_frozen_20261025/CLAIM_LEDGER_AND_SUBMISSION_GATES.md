# V4 S1–S2 manuscript claim and source ledger

**Frozen science reference:** [FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md](../../analysis/v4/FINAL_V4_BIO18_HOLM_LEGO_FREEZE_20261024_JA.md).  
**Frozen commit:** 7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159 (2026-10-24).  
**Source computations:** d2c45b00018a0c04e7c69a29857207dc1ac06d0c and Actions 40062126006, artifact 12547876704, ZIP SHA-256 6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205.  
**Aggregate JSON SHA-256:** 3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70.  
**LEGO JSON SHA-256:** bb01d3b5da03148d84b2521c9c92cf14f7f6177435e145118f6ff6df8846733f.  
**Independent frozen audit:** Actions 40063813578 (success), artifact 12548870574.  
**Immutable base:** snapshot/v4-final-bio18-holm-20261024. This manuscript branch adds text/data-to-figure scaffolding only.

| Claim ID | Manuscript section | Evidence fixed in the freeze | Allowed claim | Stop / qualifier |
|---|---|---|---|---|
| C01 | Abstract / 2.1 | S1 313,350 source WCVP-native records; S2 297,622 wild native/introduced records in two species | Source study populations | These are not all measured head sample sizes; native sources overlap |
| C02 | 2.2 / 3.1 | S1 168,907 angle-ready records / 178 taxa, 156,683 all-8 complete photographs / 174 taxa | Trait-specific denominators | Do not substitute original 170,066 / 244 from a distinct S1 cohort |
| C03 | 2.3 / 3.1 | BIO18/BIO05/RSDS/elevation model; max source VIF 4.9559 | Fixed four-exposure observational model | BIO18 is mean monthly precipitation during warmest quarter; not measured flowering rain |
| C04 | 2.4 / 3.1 | Separate 80 S1 and 128 S2 joint-hue Holm families; four pass | Four retained associations | CIs pointwise, P adjusted only within selected model/family; model chosen after exploration |
| C05 | 3.1 Table 1 | S1 between-taxon angle × BIO18 β+0.045447, CI +0.022349,+0.068545, Holm 0.00920635 | Larger image angle at wetter warm-quarter climates across taxa | Not gravity-validated nodding or rainfall protection |
| C06 | 3.1 Table 1 | Introduced arvense lightness × elev β−0.173732, Holm 0.00131233 | Within introduced stratum image lightness declines with elevation | Not pigment assay or introduction-induced evolutionary change |
| C07 | 3.1 Table 1 | Introduced vulgare chroma × RSDS β+0.134504, Holm 0.01456813 | Within introduced stratum chroma is higher with RSDS | Same photo colour extraction; not independent mechanism |
| C08 | 3.1 Table 1 | Introduced vulgare lightness × RSDS β−0.142801, Holm 0.01925414 | Within introduced stratum lightness is lower with RSDS | Same RSDS/sample/photos as C07 |
| C09 | 3.2 Table 2 | S1 28 pairs, raw between>within 21; BIO18-adjusted 24; within module RV median 0.085245 vs cross 0.001419; among 0.197455 vs 0.015453 | Partially separable observable trait association groups | RV baseline/sample size, common pixels/derived metrics |
| C10 | 3.2 | Size-preserving 280 trait-label allocations; rank 1/280 within, 3/280 among | Registered domains distinctive in frozen measured covariance | NOT genetic modularity test nor independent replication |
| C11 | 3.3 | 123/135 taxa possess all eight three-trait high/low combinations after adjustment | Broad photographed combinations exist | Sampling intensity/threshold dependence; not viability/genetic recombination |
| C12 | 3.3 | Parity +0.03640 [0.01456,0.05762]; all-three-high −0.003595 [−0.007336,−0.000049]; between contrast +0.011494 [−0.004680,+0.027401] | Higher-order statistic depends on coding | No threshold-invariant synergistic G123 |
| C13 | 3.4 | S2 adjusted RV medians for all 4 species×range strata | Descriptive colour/outline/photo covariances within strata | No replicated evolutionary module transition |
| C14 | 2.6 / 3.4 | Bract 1,168 each/39 taxa; chroma-bract 1,007/38; adjusted between RV0.153365 /within0.001796; subset max VIF6.14 | Limited supplementary image proxy | Not bract spine/stickiness/defense; not four-exposure primary claim |
| C15 | Discussion | No independent gravitational marker, pigment measurement, genetic module, field fitness, systematic phylogenetic correction | Observational claim boundary | Never call adaptation, LEGO genetic reuse or independent gains established |
| C16 | Abstract / Results 3.4 | Original Actions 12547876704 `S1_S2_extreme_LEGO_results.json`: *arvense* BIO18 222.2 native n=78,646 → 270.4 introduced n=70,964; *vulgare* 218.0 native n=62,368 → 267.7 introduced n=85,644 (mm per month) | Higher source-site median BIO18 in introduced ranges, +48.2/+49.7 mm month⁻¹ | Occurrence-site climate distribution, not matched photo-trait regression cohort or climate-niche evolution; source at same verified ZIP SHA |

## Frozen statistics used to regenerate summary figures

The machine-readable manuscript-local file [table1_holm_associations.csv](table1_holm_associations.csv) contains only the four frozen Holm-passing comparisons. It is a **display subset**, not the complete 208-test result table. Any final Figure 2 caption must name the selection rule and direct readers to the complete archived model output. [s2_bio18_range_occupancy_summary.csv](s2_bio18_range_occupancy_summary.csv) contains source-location BIO18 n, 5th/95th quantiles and medians for all four range strata. [module_summary.csv](module_summary.csv) contains descriptive module RV medians and [s2_module_summary.csv](s2_module_summary.csv) the four S2 stratum medians. These values must never replace the complete 28-pair artifact when making a full integration heatmap.

## Editorial completion gates

1. Recover the full 208-test model tables, full 28-pair raw/conditional RV outputs and original species/photo/site eligibility appendix from checksum-pinned statistical artifacts. Verify Table 1 numbers byte-for-byte.
2. Add provenance-safe figures: geographic sampling coverage, head-crop image example with original licences, complete 28-edge heatmap, and the three-way threshold sensitivity; never infer pixels/locations.
3. Transcribe image/EXIF preprocessing, classifier weight identity, observer balancing, S1 between/within design formula, software versions, spatial residual limitations and source filtering into a reproducible Supplement.
4. Audit every external scientific reference, journal formatting, acknowledgements, complete author metadata, permission/photograph attribution and public accessibility of reproducibility artifacts before claiming submission-readiness.
5. Obtain an independent critical review of whether selection-on-significance and trait-code shared geometry still permit the chosen conceptual framing. If not, revise claims while keeping the locked numerical science unchanged.

## Explicitly excluded interpretations

- adaptive rain shielding or gravity-axis nodding;
- causal effects of BIO18, RSDS, BIO05 or elevation;
- genetically independent reusable LEGO units, three-way functional G123 or independent evolutionary transitions;
- introduced-range evolution inferred from two species;
- phylogenetically corrected T2 or cross-tribe T3 replication (separate projects);
- claims that Holm adjustment also corrects the entire preceding model-selection process.

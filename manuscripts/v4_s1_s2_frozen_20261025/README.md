# GEB V4: global photographic phenomics Research Article (science frozen 2026-10-24)

**Science locked:** 24 October 2026, original V4 S1–S2 commit `7e25e409d39a7e8bb8a1ce1043d6bbfd41b7d159`. The manuscript has been revised for GEB after this freeze; numerical estimates were not modified.

**Current article title:** *Global photographic phenomics reveals precipitation-associated capitulum presentation and floral colour gradients in Cirsium*.

**Spatial-model sensitivity (SI only; 2026-10-29):** [source-verified BIO18 summary](SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md) and [original runner](https://github.com/zuizui0223/azami/actions/runs/40837107465) show that spatial terms attenuate the original BIO18 coefficient. Because climate and geography overlap, this does not causally separate or negate environmental associations. The original S1 Holm result and Main Figure 3 are retained unchanged.

**Current central argument (2026-10-29 GEB published-article editorial pass):** Global photograph-derived continuous traits allow climate–trait comparisons across a genus diversified in multiple regions. The paper foregrounds **the corrected native among-taxon BIO18 × image-angle relationship** and **three introduced-stratum colour–radiation/elevation relationships in two source-eligible species**, then asks what a historical native/introduced comparison contributes beyond among-taxon turnover. The directional rain-shielding interpretation is a hypothesis with independent experimental precedent, not measured fitness or an adaptation claim. Within-/across-domain image metric covariance and eight-state threshold configurations are descriptive context and their complete results remain in SI.

## Three independent ecological questions (two original source cohorts)

1. **3.1: Global climate association.** Native-taxon continuous photographed angle versus BIO18 among 178 taxa. Original source-verified β +0.045447, Holm P0.00920635; [spatial-model dependence](SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md) is reported in **SI only**, not recast as a new main finding.
2. **3.2: Conspecific introduction comparison.** Two taxa and four native/introduced strata, with three source-verified, Holm-supported image colour×radiation/elevation associations. This is not a test that evolution occurred after introduction.
3. **3.3: Local climate-analogue matching of multivariate floral appearance.** The original 3D climate/niche breadth test and species-centroid comparison were not statistically supported, but cross-region *same-climate-site* head phenotypes aligned more strongly than taxon-preserving null expectations (1° +0.04096, 5° +0.09433; exploratory Holm-2 P=0.010 each). The matching KDE L2 overlap is significant only **relative to its within-taxon permutation null**. Environment-moderated 3-trait pairwise covariance interactions (six comparisons) and axis-level breadth tests (three comparisons) were not Holm-supported.
Source **S1/S2** labels are unchanged, and there is **no third independent source S3**. V3's measurement and ecological motivation is retained without treating its historical regression coefficients as V4 replication. “LEGO” was an informal project nickname, not the formal hypothesis.

[Japanese outline](GEB_THREE_PART_STRUCTURE_JA.md) | [V3/V4 evidence comparison](GEB_V3_V4_MINIMAL_DELTA_STORY_AND_FIGURES_JA.md) | [Contribution pitch kept in the cover-letter draft](GEB_COVER_LETTER_CONTRIBUTION_PITCH.md).

- [H3 source-verified climate/trait morphospace aggregate](../../analysis/v4/results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json) and [H3 methodology/results SI](SI_H3_CLIMATE_MORPHOSPACE_BREADTH_20261029.md) — new 3.3, not part of the original S1/S2 Holm families.

## New source-confirmed H3c/H3d endpoints (29 October 2026)

**H3c — positive ecological pattern:** 157,822 common valid photographs/242 source taxa; **63 geographically supported taxa at 1°** and 31 at 5°. For different-sector taxa, matching the three-dimensional climate conditions of occupied cells (BIO12/BIO01/RSDS) yields closer three-dimensional photographed capitulum phenotypes (angle/chroma/aspect) than expected under within-taxon climate–trait permutation: **distance advantage +0.04096, Holm-2 P=0.010** (1°; mean 229 supported cross-sector pairs) and **+0.09433, Holm-2 P=0.010** (5°; mean 103 pairs). A fixed-bandwidth conditional 3D KDE L2 overlap comparison also deviated from its taxon-preserving null (Holm-2 P 0.010 / 0.015). **The conditional subset's absolute overlap was LOWER than the full eight-cell overlap**; this is a relative-to-permutation result, not an absolute increase in niche intersection. The H3a species-level breadth slope +0.0669 (CI −0.1046,+0.2384) and H3b species centroid ρ−0.0081 are not supported.

**H3d — not a positive claim:** Three trait-specific climatic-breadth responses and six rainfall/radiation-modulated cross-trait covariance contrasts were prespecified and tested on the same matched source; none survived their respective exploratory Holm families. In particular, the nominal 1° solar-radiation × image-angle–outline covariance P=0.016 had Holm-6 P=0.096 and a non-replicating/opposite-signed 5° estimate. No genetically encoded module, correlational selection or evolutionary convergence is established.

[H3c empirical source](../../analysis/v4/results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json) · [H3d empirical source](../../analysis/v4/results/20261029_H3D_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_SOURCE_VERIFIED.json) · [Source CI #40847696537](https://github.com/zuizui0223/azami/actions/runs/40847696537) · [H3c SI](SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md) · [H3d SI](SI_H3_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_20261029.md).

**GEB placement:** H3c is now Results **3.3** and Figure 5; H3a/H3b as context; H3d weak covariance moderation, the original 28-pair matrix and spatial BIO18 sensitivity are **SI**, while original 3.1/3.2 results and original S1 80/S2 128 Holm tests remain fixed.

## Manuscript and scientific evidence
- [GEB_RESEARCH_ARTICLE.md](GEB_RESEARCH_ARTICLE.md) — scientific working text with GEB structured Abstract, Methods, Results, Discussion and cited natural-contrast interpretation.
- [GEB_RESEARCH_ARTICLE_BLINDED.md](GEB_RESEARCH_ARTICLE_BLINDED.md) — a separately prepared identity-concealed text draft. Actual anonymous reviewer access, affiliations, title page and final journal formatting are **pending**.
- [CLAIM_LEDGER_AND_SUBMISSION_GATES.md](CLAIM_LEDGER_AND_SUBMISSION_GATES.md) — each frozen estimate to its source and non-claims.
- [Full 28-pair and configuration SI context](SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md) — all same-domain metrics, three original cross-domain signed-r pairs, descriptive RV medians, 280 group-label allocations and threshold-sensitive eight-state counts.
- [GEB published-paper writing audit (Japanese)](GEB_PUBLISHED_RESEARCH_STYLE_AUDIT_20261029_JA.md) — five source-verified GEB articles, their rhetorical structure and the exact V4 editorial decisions.
- [GEB_SUBMISSION_CHECKLIST.md](GEB_SUBMISSION_CHECKLIST.md) — structure, required figures, full 208-test table and unresolved submission gates.
- [S2 seven-species eligibility table and stronger two-range gate](SI_S2_FOCAL_SPECIES_SELECTION_AUDIT.md) — source-real WCVP 350-name screen, all seven two-range candidates, and exactly two species meeting coverage in **each** native/introduced stratum.
- [S2 study-design synthesis: historical rationale × empirical eligibility](S2_SELECTION_RATIONALE_20261025.md) — two defensible reasons for selecting the same two focal species without claiming random selection or new niche/trait fits.
- [S2 range-comparison evolutionary rationale and source-verified literature (Japanese)](S2_NATIVE_INTRO_EVOLUTIONARY_RATIONALE_REFERENCES_JA.md) — why older native vs secondary introduced residence histories are compared, direct C. arvense experimental precedent, and why this is NOT rapid species radiation or a rainfall-niche analysis.
- [Interaction and landing-filter literature audit (Japanese)](GEB_INTERACTION_FILTER_LITERATURE_AUDIT_20261027_JA.md) — primary studies on visibility, landing posture, visitor selectivity, nectar access and capitulum seed predation; explicit V4 cross-domain signs versus hypotheses.
- [Conditional trait interactions and competing causal hypotheses (Japanese)](GEB_TRAIT_COVARIATION_BIOTIC_INTERACTIONS_20261027_JA.md) — distinguishes photo covariance, combinatorial occupancy, functional interactions, selection and genetic autonomy; proposes field tests for mutualist–antagonist–rainfall trade-offs.
- [S1 BIO18 image-angle full raw/BH/Holm context](SI_S1_BIO18_ANGLE_RAW_BH_HOLM_CONTEXT.md) — the two original scales and all 80 tests in SI; Main Figure 3 displays only the Holm-supported among-taxon result.
- [S2 native/introduced colour-coefficient context](SI_S2_COLOUR_NATIVE_INTRO_FIXED_HOLM_CONTEXT.md) — three original introduced Holm-positive colour associations shown next to the corresponding native coefficients; no new model or rainfall-niche comparison.
- [S1–S2 scope lock](SCOPE_LOCK_20261025_JA.md) — retain 24 October models, S1 macroecological results, three S2 introduced-stratum colour associations and phenotype covariance, without precipitation-niche comparisons.
- [Deferred source-real S2 climate-niche exploration](../../legacy/v4_s2_precipitation_niche_20261025/README.md) — original exploratory results preserved outside the active article.

## Exact source-level statistical evidence
- Actions original success **40062126006**, artifact **12547876704**. Original archive ZIP SHA-256:
  `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`
- `table1_holm_associations.csv`: four accepted associations after separate S1 (80 tests) and S2 (128 tests) Holm corrections.
- `module_summary.csv`, `s2_module_summary.csv`: display-only 24 October fixed phenotypic covariance summaries, not the complete hypothesis family.
- `build_full_supplement.py`: SHA-checked source-regeneration of **all 208** primary BIO18-family tests and 28 raw/adjusted RV construct pairs.
- `plot_frozen_summary.py`: the historical four-Holm forest and domain-median charts are **supplementary** descriptive displays.
- `plot_bio18_angle_part2.py`: the unchanged original Holm-supported BIO18×image-angle coefficient and the supplementary within/among scales; SI includes both scales and 80 raw/BH/Holm rows.
- `plot_s2_native_introduced_colour.py`: Main Figure 4 **source-verified six original coefficients**, generated from the exact original S2 source CSV by the existing SHA-locked validator. This is not a cross-range hypothesis test or rainfall-niche comparison.
- `build_full_supplement.py`: preserves and regenerates the **complete original 28-pair matrix in SI** (with no source-data changes). `plot_cross_functional_pairs.py` verifies and plots the **three original cross-domain signed-r pairs as a supplementary exploratory asset** (historical filename still begins `figure2_`); no new statistical tests. Main Figure 2 is reserved for a **not-yet-source-generated** geographic coverage map; no invented locations/images.

## Explicit scientific separation

**Part 1** reuses native common photographs from original source S1 for visible trait integration; **Part 2** models the native S1 cohort and emphasizes the Holm-supported image-angle×BIO18 association; **Part 3** uses the two-species×two-range source S2 for **within-stratum** photographed colour gradients. The parts are not independent datasets or new multiple-testing families. **No native-versus-introduced precipitation/climate niche comparison belongs to this article**, and all 24 October source results and Holm coefficients are unchanged.

Historical V3 numerical results belong to `legacy/v3/`. Separate 25 October S2 climate-niche analyses are preserved under `legacy/v4_s2_precipitation_niche_20261025/`; they are **not** current GEB findings, submission figures or a prerequisite for the frozen numerical replay. T2/T3 are not validation layers in this S1–S2 manuscript. Private identities, coordinates and unlicensed images are not redistributable.

- `plot_bio18_geographic_sensitivity.py`: **SI-only spatial-model sensitivity forest** (source-locked original baseline, fixed sector/sphere controls and eastern holdout under photo-weighted versus cell-balanced taxa). The source-audited aggregate JSON is committed separately; original individual metadata remain private.

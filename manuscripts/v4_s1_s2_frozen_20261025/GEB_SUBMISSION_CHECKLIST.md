# GEB Research Article preparation status — updated 30 October 2026

Journal source: [Wiley Global Ecology and Biogeography author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html) (rechecked 2026-10-30).

**Current stage:** source-verified figures and editable main/SI review candidates assembled. No journal submission, permanent data deposit or DOI publication has occurred.

**Repository merge / Zenodo boundary (2026-10-30):** PR #159 is merged into `main` at `e2229a2`; exploratory PR #160 remains open, unmerged and separately reviewed. [Zenodo V4 data-only preparation](../../docs/ZENODO_V4_RELEASE_HANDOFF.md) specifies frozen aggregate candidates, rights, creator metadata and byte verification. A V4 DOI and anonymous reviewer access remain absent. CC BY 4.0 is now approved for eligible author-owned aggregate results, and the 65-file rights-scoped candidate passes 9/9 offline replay commands with all five CSV tables byte-identical. Third-party rights clearance and creators still require confirmation. The remaining author and access decisions are listed below; scientific limitations remain explicit in the manuscript and SI.

## Latest explanation-variable decision (2026-10-31)

**For next submission revision (all elevation-free):** **3.1 / 3.2 = BIO18 + BIO05 + RSDS**; **3.3 = original BIO12 + BIO01 + RSDS**, keeping the original H3c endpoint and its exploratory Holm-2 separate. The earlier choice of BIO18 + BIO01 + RSDS for *all* Results has been superseded. See the [canonical updated decision and 3.2 five Holm-supported associations](GEB_CLIMATE_ONLY_MODEL_DECISION_20261031_JA.md). In the new 3.1, BIO18 between-taxon β=+0.04248349, new Holm P=0.01570433 (169,952 photo records / 178 native taxa). In new 3.2, five of 96 source-verified tests pass Holm across image lightness, hue, angle and outline; no native-versus-introduced interaction test was performed. The **existing** main article/figures/Table 2 and the original numerical SI still reflect the old four-variable freeze; they must be updated together before submission. The 2026-10-24 historical frozen inference and V4 Zenodo candidate remain separate and unchanged.

## Format and scope

- **Research Article**, manuscript body approximately 5,000 words, about 50 references and **6–8 display pieces** (figures and tables combined; flexibility allowed).
- **Structured Abstract ≤300 words** with the required headings: Aim, Location, Time period, Major taxa studied, Methods, Results, Main conclusions. The GEB Abstract must be revalidated against this current draft (≤300 words, including required subheadings).
- Running title <40 characters; 6–10 keywords in alphabetical order.
- Double-anonymous review: blinded main file separated from author-identifying title page (author names, affiliations, corresponding author email, acknowledgements/CRediT/funding as appropriate).
- Manuscript editable, high-resolution images and complete figure/table legends, continuous line numbers in the final editable document. Main figures embedded in the manuscript; Supporting Information in separate files.
- Main file order: Title, running title, abstract/keywords, body, references, data/code availability, tables, figure legends and embedded figures. Uniform author-date references.
- Cover letter uploaded separately as a PDF, including a <250-word paragraph explaining interest to GEB readership.
- Data and code availability declaration is mandatory. Pre-print/anonymous review availability are not substitutes for a durable, rights-compliant data/code release.

## Manuscript narrative: three ecological questions

**3.1 Climate and head presentation:** Original native-taxon BIO18 × image angle, n=168,907/178 taxa, β=+0.045447, Holm P0.00920635. Geographic sensitivity is reported in [SI](SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md) as model-dependence, not a causal separation of climate from space.

**3.2 Flower colour after colonisation of new regions:** Two species × four range strata, 297,622 original wild source observations; three introduced-stratum colour×radiation/elevation associations survived the original 128-test Holm family. No formal range×environment interaction was estimated.

**3.3 Local climatic analogues and floral morphospace:** H3c detected greater joint angle–chroma–aspect resemblance at climate-matched cells than expected from within-taxon vector reassignment (distance advantage 1° +0.04096; 5° +0.09433, exploratory Holm-2 P=0.010 at each resolution). Main Figure 5 presents this distance result. The primary analysis retains 63 taxa / 7,117 taxon-cells; the separately verified map shows 153,843 contributing source photographs before cell aggregation and rarefaction. Taxa are assigned to different dominant longitude sectors; their matched cells need not occur in different sectors. Climatic breadth, centroids, KDE and covariance endpoints remain in SI. [Detailed source audit](SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md).

**H3c robustness gate:** Six further fixed random-stream checks retained the positive *direction* at both resolutions (6/6), but nominal P<0.05 in only 2/6 of the 1° series versus 6/6 of 5° series. At 1°, excluding chroma or the central longitude sector reversed the distance contrast. KDE overlap was not consistently supported. See [full post-positive sensitivity SI](SI_H3C_POSTPOSITIVE_ROBUSTNESS_20261029.md). This is not another independent Holm-2 family.
The weak exploratory 28-pair shape/colour/angle covariance remains in SI; “LEGO” is not a formal study hypothesis. Original S1/S2 frozen results unchanged.

## What the frozen results can support

**Macroecological axis:** among native taxa, BIO18 relates positively to **camera-frame** capitulum angle in one of 80 S1 Holm-controlled tests. S2 estimates **phenotype–environment regressions within** two species×native/introduced strata: three introduced-stratum colour slopes pass Holm; no original S2 head-angle term passes Holm. **No wet/dry climate-niche occupancy comparison or introduction-driven phenotypic evolution is claimed.** S1 native photos can overlap S2 native photos.

**S2 colour asymmetry:** only introduced S2 strata contain Holm-surviving colour coefficients (3/128), but the native *C. vulgare* chroma×RSDS coefficient is similar to the introduced coefficient and passes BH, not Holm. *C. arvense* lightness×elevation has opposite-signed native/introduced point estimates. These are not tests of direct slope differences; [all six source-locked original rows](SI_S2_COLOUR_NATIVE_INTRO_FIXED_HOLM_CONTEXT.md) are in SI. The manuscript must never claim acquired colour responsiveness, identical native/introduced response slopes, or introduction-caused phenotypic evolution. **No precipitation niche is compared.**

**Phenotypic organisation:** eight visible constructs/28 pairs on the common S1 156,683-photo/174-taxon cohort show heterogeneous association structure; 24/28 adjusted RV are between>within; label-shuffle reference is conditional on measured coordinate definitions. **Do not claim** genetic autonomy of measured modules, functional synergy or independent evolutionary origins from common-source image traits.

**Now empirically tested geographic dependence:** [2026-10-29 original-cohort BIO18 sensitivity](../../analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED_JA.md) reprocessed the exact 168,907 images/178 taxa; sector adjustment β+0.020035 (HC3 CI includes zero), spherical geography β−0.003295 (CI includes zero), and eastern-region exclusion β−0.010372. The original Holm-supported β+0.045447 remains a global observational association; **it is not a geographically universal rain response**. The earlier BIO12 geographical result is a different model/cohort and not independent confirmation.

**Alternative explanation:** camera angle and lighting, taxon/geographic confounding now tested on original BIO18, shared-pixel mathematical covariance, imperfect taxonomic and temporal resolution. Selected-model Holm correction is **not** adjusted for the earlier explorations leading to the fixed model. All significant results and all 208 tests are visible in the numerical supplement; don't conceal the full family.

## Seven main display pieces for GEB (five figures + two tables)

**Methodology plus ecology:** show how continuous photographed phenotypes are measured and where each analytical cohort was sampled, then the precipitation, colour and local-resemblance results. The [cover-letter draft](GEB_COVER_LETTER_CONTRIBUTION_PITCH.md) has a single 206-word interest paragraph.

1. **Figure 1:** Actual archived source photographs illustrate detection, image-referenced angle, floral pixels and projected outline. Current cohort membership is not asserted for the illustrative images. All three photo licences were verified from their source pages and API records: one CC BY 4.0 and two CC0 1.0. Credits and modifications are in the caption and [rights manifest](FIGURE1_PHOTO_CREDITS_20261030.json). PNG 600 dpi; PDF/SVG masters; 168 × 202 mm.
2. **Figure 2:** Six actual source-verified sampling maps, pooled into 5° display cells and omitting cells with fewer than five records. Original S1 angle and six S2 colour cohorts were reproduced before export, and H3 geographic eligibility was verified. Mollweide equal-area projection, shared logarithmic scale, graticules and 10,000 km equatorial scale; 168 × 196 mm. [Coverage source](../../analysis/v4/results/20261030_SAMPLING_COVERAGE_SOURCE_VERIFIED.json) / [render manifest](figure2_render_manifest.json). The historical `figure2_cross_functional_trait_pairs` remains an SI figure.
3. **Figure 3 (native climatic association):** Main figure displays **only the frozen, Holm-supported BIO18 × image-referenced angle coefficient**. The extra sector/spherical-geography and regional-holdout models belong to [spatial-model sensitivity in SI](SI_S1_BIO18_SPATIAL_SENSITIVITY_20261029.md), not the main figure. The [full S1 SI](SI_S1_BIO18_ANGLE_RAW_BH_HOLM_CONTEXT.md) retains all original 80 unadjusted/BH/Holm rows and both observational scales.
4. **Figure 4 (native/introduced):** Three selected photographed colour×environment pairs in native/introduced strata, six exact original coefficients with CI; no range×environment interaction inference.
5. **Figure 5 (H3):** Source-verified climate-matched joint-trait distance advantage versus the within-taxon permutation-null mean at 1° and 5°. Resolution-specific standardisations prevent direct effect-size comparison between grains. Gaussian-KDE similarity, climatic breadth and the original eight high/low configurations remain in SI. No fictitious photo scatter or evolutionary convergence inference.
6. **Table 1:** Audited source totals and image-angle eligible observations by native/introduced stratum for *C. arvense* and *C. vulgare* (78,646/70,964 and 62,368/85,644; native/introduced totals 141,014/156,608). Weak r pairs and domain RV medians moved to SI.
7. **Table 2:** Four original Holm-supported environment associations, with Holm rather than raw P displayed in the main table.

**SI figures and tables:** three original cross-domain correlation pairs and  full 208-test family with original/BH/Holm P, the two-scale BIO18 angle forest, all 28 RV pairs (original full matrix and CSV), 280 fixed-size label ranks, eight-state configurations with threshold sensitivity, original domain RV medians, previous four-positive forest and S2 range-stratum RV medians. See SI_V4_FULL_COVARIATION_AND_CONFIGURATION_CONTEXT.md. Earlier V3 numerical plots are not re-labelled as V4.

## Completed preparation and remaining submission decisions

**Assembled:** five main figures in PNG/PDF/SVG; two editable main tables; complete captions with photographic attribution; editable main/SI Word candidates with continuous main-text line numbering; 11 SI appendices; five original numerical CSV exports; aggregate verification/rendering code and source receipts. Main/SI document assembly records source-text and numerical preservation, figure identity and reviewer navigation. Generated Word/ZIP packages are excluded from git and delivered separately.

**Source coverage verified:** [Actions #41020812565](https://github.com/zuizui0223/azami/actions/runs/41020812565), artifact 12961447511. Coarse-cell counts retain every cohort denominator and explicitly account for suppressed records. No exact coordinates, observation IDs, uploader IDs or environmental pixel records are included in the map export. Three photographic illustrations have separately documented attribution rights.

**Aggregate replay verified:** the candidate reproduces the 208-test table, 28-pair table, six original colour coefficients, the 80-test native family and two-scale BIO18 table. Nine verification/render commands pass. It starts from archived aggregates and **does not rerun observation-level fits, raw-image measurement, bootstrap or permutation analyses**. The review copy neutralises identifying repository references but retains the legally required MIT copyright notice; full anonymity is not certified. See [package scope](reproducibility_v4/AGGREGATE_RELEASE_CANDIDATE_README.md).

**Still required before actual submission:**

1. Confirm the author list/order, affiliations, corresponding author, acknowledgements, funding, conflicts and other journal declarations; create the separate author title page and final cover-letter PDF from those confirmed details.
2. The **CC BY 4.0 policy has been approved for original author-owned, redistributable numerical aggregates**. Complete file-level rights review and creator information, then establish durable deposition and suitable reviewer access. No public DOI or anonymous hosted link is yet claimed; MIT code and third-party photo/taxonomy/climate rights remain separate.
3. Perform the final author sign-off on the main/SI package and submission declarations. Scientific interpretation remains observational; independent image-measurement validation and experimental tests of protection/fitness are limitations or future work, not evidence manufactured by formatting the package.

**State:** the original numerical science is source-verified, and the manuscript is an editable review candidate. It has not been submitted or accepted. Author/access decisions above remain open.

## Frozen GEB inferential scope

**Included:** 24 October 2026 V4 S1 (native-taxon angle–climate and 28 photographic trait-pair analyses), original S2 128 within-stratum hypothesis tests, three Holm-supported introduced-stratum colour results, and original S2 phenotype-domain RV medians. **Excluded:** native-versus-introduced precipitation-niche, climatic centroid/novelty/expansion, background-matched niche analyses, and post-freeze within-species range × BIO18 interactions. Earlier source experiments remain archived at [`legacy/v4_s2_precipitation_niche_20261025/`](../../legacy/v4_s2_precipitation_niche_20261025/README.md), but should not be referenced as manuscript findings.

**Editorial comparison:** [GEB published-article writing audit](GEB_PUBLISHED_RESEARCH_STYLE_AUDIT_20261029_JA.md). GEB Research Article structured Abstract and broad-scale generality follow the [official author guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html).

**SI-only source-grounded spatial-model sensitivity:** [aggregate JSON](../../analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json), [secure-source CI](https://github.com/zuizui0223/azami/actions/runs/40837107465). `plot_bio18_geographic_sensitivity.py` generates an **SI sensitivity figure** in the standard CI artifact (separate from main Figure 3); no individual coordinates/observer data shared.

**H3 original-source CI:** [Actions #40841007115](https://github.com/zuizui0223/azami/actions/runs/40841007115), artifact 12883021095; [machine-readable anonymous source-verified H3](../../analysis/v4/results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json). All H3 hypotheses are exploratory and separate from the 24 October S1/S2 frozen Holm tests.

## New source-confirmed H3c/H3d endpoints (29 October 2026)

**H3c — positive ecological pattern:** 157,822 common valid photographs/242 source taxa; **63 geographically supported taxa at 1°** and 31 at 5°. For different-sector taxa, matching the three-dimensional climate conditions of occupied cells (BIO12/BIO01/RSDS) yields closer three-dimensional photographed capitulum phenotypes (angle/chroma/aspect) than expected under within-taxon climate–trait permutation: **distance advantage +0.04096, Holm-2 P=0.010** (1°; mean 229 supported cross-sector pairs) and **+0.09433, Holm-2 P=0.010** (5°; mean 103 pairs). A fixed-bandwidth conditional 3D KDE L2 overlap comparison also deviated from its taxon-preserving null (Holm-2 P 0.010 / 0.015). **The conditional subset's absolute overlap was LOWER than the full eight-cell overlap**; this is a relative-to-permutation result, not an absolute increase in niche intersection. The H3a species-level breadth slope +0.0669 (CI −0.1046,+0.2384) and H3b species centroid ρ−0.0081 are not supported.

**H3d — not a positive claim:** Three trait-specific climatic-breadth responses and six rainfall/radiation-modulated cross-trait covariance contrasts were prespecified and tested on the same matched source; none survived their respective exploratory Holm families. In particular, the nominal 1° solar-radiation × image-angle–outline covariance P=0.016 had Holm-6 P=0.096 and a non-replicating/opposite-signed 5° estimate. No genetically encoded module, correlational selection or evolutionary convergence is established.

[H3c empirical source](../../analysis/v4/results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json) · [H3d empirical source](../../analysis/v4/results/20261029_H3D_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_SOURCE_VERIFIED.json) · [Source CI #40847696537](https://github.com/zuizui0223/azami/actions/runs/40847696537) · [H3c SI](SI_H3_CLIMATE_ANALOGUE_CONDITIONAL_OVERLAP_20261029.md) · [H3d SI](SI_H3_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_20261029.md).

**GEB placement:** H3c is now Results **3.3** and Figure 5; H3a/H3b as context; H3d weak covariance moderation, the original 28-pair matrix and spatial BIO18 sensitivity are **SI**, while original 3.1/3.2 results and original S1 80/S2 128 Holm tests remain fixed.


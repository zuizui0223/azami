# azami — global continuous thistle image phenomics

This repository is organized around the Chapter 1 submission manuscript:

**Global citizen-science images reveal continuous within-taxon variation and environmental structure in thistle capitulum traits**

The canonical submission-facing manuscript entry point is [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md). The current paper no longer uses the former lability-axis story as a biological headline.

## Current paper in one sentence

Public biodiversity photographs are converted from heterogeneous species records into repeated **continuous head-level phenotype observations**, revealing that thistle capitulum orientation, visible colour, gross outline and fine involucral architecture show different environmental signatures within and among source-assigned taxa rather than forming one common climate-associated morphology.

## Scientific scope

YOLO11n is used only to localize visible capitula. Deterministic image-processing functions then return continuous numerical measurements for each assessable detected head.

Primary modules:

- **Orientation** — signed head-axis angle relative to image vertical;
- **Visible corolla colour** — CIELAB lightness, chroma and circular hue components;
- **Capitulum outline** — aspect ratio, circularity, solidity and width-profile CV.

A high-resolution layer additionally measures three exploratory two-dimensional involucral contour proxies:

- projection roughness;
- outward spread fraction;
- maximum spine-like projection relative to head radius.

The study therefore follows the explicit route:

`public photographs -> detected capitula -> continuous head-level measurements -> repeated within-taxon trait distributions -> multiscale environmental analysis`

## Current headline results

### 1. Species/taxon means conceal substantial visible variation

Across the nine primary endpoints, **0.589–0.931** of visible image variance occurs below source-assigned taxon means.

- among photographs within taxa: **0.440–0.691**;
- among heads within photographs: **0.143–0.379**;
- one-head-per-photo sensitivity: **0.582–0.899** below taxon means;
- equal-10-photo-per-taxon sensitivity: median **0.528–0.879**.

These are **visible image-phenotype variance components**, not genetic variance estimates.

### 2. Within-taxon environmental structure is trait specific

In the 46,276-observation spatially thinned primary cohort, eight of 36 endpoint-component × CHELSA rows pass BH correction.

Key non-circular results include:

- orientation angle increases with annual mean temperature (BIO1);
- corolla chroma increases with annual precipitation (BIO12);
- outline aspect ratio increases with temperature seasonality (BIO4) and decreases with annual precipitation (BIO12).

Circular hue components also vary with temperature seasonality and precipitation structure and are interpreted jointly rather than as independent named colours.

### 3. Spatial modelling concentrates the strongest support in orientation and colour

Grouped SPDE-INLA models retain the clearest stable environmental structure for:

- orientation angle;
- corolla lightness and chroma;
- circular hue components.

Temperature and soil pH contribute to visible-colour structure, whereas most gross-outline traits receive little global BH support after spatial modelling.

### 4. Fine involucral architecture tracks temperature seasonality

In the <=10 km high-resolution cohort (904 observations, 165 taxa), all three final auxiliary contour proxies increase with BIO4 temperature seasonality and pass BH correction across the 12-test auxiliary family.

These variables describe **outward contour architecture** only. They are not direct measurements of botanical spine length, bract recurvature, stiffness or defence.

### 5. Among-taxon environmental sorting is concentrated mainly in orientation and visible colour

Among 148 taxa complete for all primary traits and environmental variables, 10,000 trait-label permutations show non-random environmental sorting mainly for:

- orientation;
- chroma;
- circular hue components.

Most gross-outline traits are weaker. Direct-backbone sensitivity narrows the supported set further.

### 6. Historical analyses remain sensitivity tests, not a resolved species-tree result

Six taxon-level climate associations are BH-supported in all 50 randomized Pagel-lambda PGLS trees. However, only **54 of 216** atlas taxa are direct dated-backbone tips; the remainder require within-genus grafting.

Accordingly, the tree layer is reported as **historical-placement sensitivity**, not definitive removal of phylogenetic non-independence.

### 7. Repository-computable robustness checks do not overturn the main results

Completed audits include:

- residual Moran's I after diagnostic SPDE refits: approximately **-0.0096 to 0.0030**, with no permutation P < 0.05;
- leave-one-broad-region-out trait-rank stability: minimum endpoint-specific Spearman rho **0.856–0.972**;
- simultaneous collapse of eight WCVP synonym conflicts: no sign changes and no loss of the eight primary BH-supported component rows;
- 10,000-permutation null tests for among-taxon environmental sorting.

## Biological interpretation

The current manuscript does **not** support one universal warm-, dry- or pollinator-associated thistle syndrome.

Instead, the capitulum behaves as a **composite reproductive phenotype**:

- orientation is most clearly associated with reproductive microclimate/exposure axes;
- visible colour is structured across temperature, precipitation and soil dimensions and may integrate abiotic and biotic processes;
- gross outline is comparatively weak after spatial control;
- fine involucral contour architecture shows a coherent temperature-seasonality association.

These patterns generate mechanistic hypotheses for field or experimental work; they do not demonstrate the mechanisms themselves.

## Claim boundary

Azami Chapter 1 does **not** claim demonstrated:

- phenotypic plasticity;
- local adaptation or selection;
- pollinator causation;
- antagonist-mediated defence;
- rain-protection adaptation;
- evolutionary rate differences;
- adaptive radiation;
- a resolved *Cirsium* species tree;
- direct botanical spine/phyllary measurements from the auxiliary image proxies.

Image vertical is also a reproducible image reference rather than a direct inclinometer measurement of gravity.

## Azami -> EAzami boundary

Azami is the **global observational discovery layer**. It ends after freezing:

`continuous visible capitulum variation -> within/among-taxon environmental structure`

Azami does **not** own inference about:

- pollinator- or antagonist-mediated selection;
- defensive function of spine/phyllary or sticky traits;
- East Asian/Japanese phylogenetic or biogeographic history;
- repeated loss/regain of floral colour;
- repeated/parallel erect <-> nodding transitions;
- trait-specific evolutionary rates or modular evolvability;
- rapid/adaptive radiation.

Those questions belong to [`zuizui0223/EAzami`](https://github.com/zuizui0223/EAzami), the downstream **mechanism + evolutionary-history zoom**. The observation-to-mechanism handoff is fixed in [`manuscript/EAzAMI_MECHANISM_REDUCTION_BOUNDARY_2026-08-20.md`](manuscript/EAzAMI_MECHANISM_REDUCTION_BOUNDARY_2026-08-20.md).

## Start here

1. [`manuscript/SUBMISSION_MANUSCRIPT.md`](manuscript/SUBMISSION_MANUSCRIPT.md) — canonical manuscript-facing story, section order and current headline.
2. [`manuscript/00_title_abstract.md`](manuscript/00_title_abstract.md) through [`manuscript/06_references.md`](manuscript/06_references.md) — current submission text and submission-wide bibliography.
3. [`manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`](manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md) — frozen cohort names, counts and permitted analyses.
4. [`manuscript/final_claims.json`](manuscript/final_claims.json) — machine-readable numerical claim registry.
5. [`manuscript/FIGURE_TABLE_MAP.md`](manuscript/FIGURE_TABLE_MAP.md) — authoritative main/Supplement figure and table crosswalk.
6. [`manuscript/MAIN_FIGURE_CAPTIONS.md`](manuscript/MAIN_FIGURE_CAPTIONS.md) — captions for the six main figures.
7. [`analysis/ch1/pipeline.json`](analysis/ch1/pipeline.json) — active executable submission-analysis map.
8. [`analysis/ch1/run_submission.py`](analysis/ch1/run_submission.py) — structural/claim checks and status summary.
9. [`manuscript/supplement/`](manuscript/supplement/) — Supporting Information tables, figure inputs and provenance manifests.
10. [`manuscript/EXTERNAL_COMPLETION_GATES.md`](manuscript/EXTERNAL_COMPLETION_GATES.md) — completed versus genuinely external validation gates.

## Main figure logic

The current submission is organized around six main figures:

1. **image -> phenotype -> ecological scale** — actual photographs, YOLO localization and continuous measurement route;
2. **geographic sampling and analytical domain** — global sampling, filtering streams and BIO1 × BIO12 coverage;
3. **nested visible variance** — taxon -> photograph -> head decomposition;
4. **taxon-level trait architecture** — 148-taxon PCA using all nine primary endpoints;
5. **environmental effect sizes** — primary within-taxon, grouped-SPDE and high-resolution coefficient families;
6. **scale specificity + historical sensitivity** — within-taxon versus randomized-tree among-taxon climate structure.

The complete figure/table contract is in [`manuscript/FIGURE_TABLE_MAP.md`](manuscript/FIGURE_TABLE_MAP.md).

## Frozen implementation and provenance

Exact numbered files under `ch1_global/v2/` remain in place because executed workflows and artifact provenance cite those paths. Do not rename or rewrite them simply to simplify the repository tree.

The active author/reviewer interface is `analysis/ch1/`, while manuscript-facing controls live under `manuscript/`.

Large analysis products are often stored as hash-verified GitHub Actions artifacts rather than duplicated in Git. Figure manifests record frozen workflow runs, artifact IDs and SHA-256 digests.

## Current cohorts

- **balanced image-comparison atlas:** 3,725 observations/photos; 6,626 heads; 216 source-assigned taxa;
- **exhaustive detector-positive layer:** 406,582 observations; 286 taxa;
- **coordinate-usable layer:** 392,989 observations; 271 taxa;
- **<=10 km positional-accuracy layer:** 297,293 observations; 259 taxa;
- **exhaustive spatially thinned primary:** 46,276 observations; 259 taxa;
- **grouped SPDE-INLA complete cases:** 31,666–34,472 observations and 139–141 taxa per endpoint;
- **high-resolution involucre layer:** 1,443 usable heads; 1,292 observations; 210 taxa;
- **high-resolution <=10 km inferential family:** 904 observations; 165 taxa;
- **independent detector audit:** 1,000 source images; 323 species; 250 double-labelled images.

## Reference source

[`manuscript/06_references.md`](manuscript/06_references.md) is the **submission-wide bibliography**. Older Introduction-only reference notes and BibTeX files are drafting provenance and must not replace the canonical manuscript reference list.

## Remaining scientific submission gates

Only two scientific blockers remain genuinely external to repository computation:

1. adjudicated human bounding boxes for the independent detector precision/recall audit;
2. genuinely independent reference measurements for orientation, colour and outline.

Until those are complete, production overlays, detector-development metrics and horizontal-mirror repeatability are technical checks rather than biological accuracy validation.

Administrative finalization still includes final author/CRediT metadata, funding, acknowledgements, nomenclatural notes and the immutable DOI-backed release.

## Withdrawn result

The former raw absolute-slope lability relation (`rho = -0.333`) and median-split quadrants were withdrawn after a precision-confounding audit. Precision-aware reanalysis found no common cross-module coupling. The former result remains provenance/QA only and must not re-enter the biological headline.

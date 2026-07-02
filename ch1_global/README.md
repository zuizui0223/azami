# Chapter 1 — Global image-derived floral trait macroecology in *Cirsium*

This document redefines Chapter 1 from a single-trait study of head orientation
into a multi-trait macroecology study built on a standardized image ROI. It
supersedes the "Ch.1-Global = orientation classifier" framing in the root
README's chapter table; see that table for how this chapter now relates to
Chapters 2-4.

## Versioning: v1 (frozen) → v2 (this redesign)

- **v1** = the orientation-only pipeline, frozen and preserved under
  [`v1/`](./v1/) (see [`v1/README.md`](./v1/README.md)). Nothing is deleted;
  v1 is the reproducible baseline the redesign inherits working code from.
- **v2** = the multi-trait macroecology plan in this document (not yet
  implemented; roadmap in §4).
- Before writing anything up, read [`METHODS_AUDIT.md`](./METHODS_AUDIT.md):
  it grades v1's peer-review vulnerabilities (three of the five macro-analysis
  foundations currently **fail**) and maps each to the v2 remediation. The
  v2 design exists to move every foundation to "pass."

## 0. Why the redefinition

The existing pipeline (`v1/01_annotate_and_train_yolo.py` →
`v1/02_crop_heads_with_yolo.py` → `v1/03_train_head_direction_classifier.py` →
`v1/04_predict_head_direction_global.py` → `v1/05_rf_glm_gam_env_analysis.R`)
detects capitula, crops them, and classifies orientation alone (Acc 0.896,
macro F1 0.861 on the current validation split — a number that is **not
peer-review-defensible as reported**, see `METHODS_AUDIT.md` A1/A2). That
pipeline is the right backbone, but as a single-trait study it (a) cannot say
whether orientation is one part of a repeated multivariate combination of
traits or an isolated pattern, and (b) has no leave-one-species-out or
geographic-block validation, so it cannot rule out "the model is reading
species identity, not the trait."

Chapter 1 is redefined as: **what repeated combinations do image-derived
capitulum traits form across global *Cirsium*, and how do those combinations
and the individual traits relate to climate, terrain, and geography?**
Chapter 1 does not claim adaptive function or causal mechanism; it reports
macroecological pattern only, on the same non-causal terms already used for
`floral_trait_architecture_protocol.md`.

## 1. Aim, hypotheses, and paper structure

### 1.1 Central question

> Across the global range of *Cirsium*, what repeated combinations do
> image-derived capitulum traits (orientation, colour, involucral cover,
> head shape, corolla:involucre display) form, and how do those traits and
> combinations correspond to climate, terrain, and geographic space?

### 1.2 Chapter-1-specific hypotheses (pattern-level, non-adaptive)

These replace/refine H2 from the root README's hypothesis table for the
scope of this chapter. H2 in its adaptive phrasing ("nodding evolved and
persists because of X") remains a Chapter 2/3 claim; Chapter 1 only
establishes the descriptive pattern that a later chapter can test causally.

| ID | Hypothesis | Falsifiable prediction |
|----|---|---|
| H1-Ch1a | Image-derived capitulum traits do not vary independently; they form a small number of recurring multivariate combinations. | PCA/FAMD on the trait set has a small number of axes explaining most variance, and/or clustering finds well-separated trait syndromes rather than a diffuse cloud. |
| H1-Ch1b | Trait combinations (not just single traits) show non-random structure with climate, terrain, and geography. | RF/GLM/GAM of syndrome axes or cluster membership against CHELSA + topography + soil + space outperform null/space-only models; effect directions are consistent across the sensitivity levels in §3.4. |
| H1-Ch1c | Models read the trait itself, not species identity. | Leave-one-species-out validation (per trait) holds up materially above the majority-class/mean baseline on held-out species; large per-species accuracy collapses are reported, not hidden in an aggregate score. |
| H1-Ch1d | Involucral cover and corolla:involucre display are related but not interchangeable axes. | Their correlation is reported explicitly (§2, `corolla_involucre_display_ratio` notes); if r > ~0.9 across the global sample this is stated as a limitation, not silently merged. |

### 1.3 Proposed paper structure

1. Introduction — capitulum trait diversity in *Cirsium*, gap: prior work
   (including this project's own earlier orientation-only analysis) treats
   traits singly; motivate a joint, image-derived, non-causal macroecological
   description as the necessary first step before phylogenetic (Ch.2) and
   functional (Ch.4) claims.
2. Methods
   1. Image sourcing and flowering filter (existing iNaturalist pipeline)
   2. ROI-A / ROI-B detection and standardization (§2)
   3. Per-trait extraction models (§3.1-3.3)
   4. Human-label validation (§3.2)
   5. Group cross-validation and leave-one-species-out design (§3.3)
   6. Environmental/geographic analysis: per-trait RF/GLM/GAM, trait-syndrome
      PCA/clustering, sensitivity analysis (§3.4)
3. Results
   1. Detection and extraction performance (incl. LOSO and per-species
      breakdown, not just pooled accuracy)
   2. Trait syndrome structure (PCA/cluster axes, what traits load together)
   3. Environmental/geographic correspondence of individual traits and of
      syndromes
   4. Sensitivity: photo-level vs species-level vs species x grid-level
4. Discussion — explicitly bounded: macroecological pattern only; state what
   would be required (Ch.2 phylogeny, Ch.3 controlled field measurement,
   Ch.4 manipulation) to move from pattern to mechanism.
5. Limitations — citizen-science sampling bias, uncalibrated colour, camera
   tilt/viewing-angle confounds, ROI detection error propagating into
   downstream trait estimates.

## 2. ROI design

Two ROI types are extracted from every YOLO-detected capitulum, replacing the
current single "capitulum crop":

- **ROI-A (standard head ROI)**: capitulum + involucral bracts + pedicel
  base, cropped from the multi-class detector (capitulum, involucre,
  pedicel-base classes; see `capitulum_bbox`/`involucre_bbox` in the trait
  dictionary) with the existing padding logic (`PAD_RATIO`) applied to the
  union box. ROI-A is the input for orientation, cover, shape, and display
  traits.
- **ROI-B (corolla-only colour ROI)**: derived from ROI-A by masking out the
  involucre/background classes from the segmentation model (§3.1), leaving
  only corolla pixels. ROI-B exists only to prevent green bract/leaf pixels
  from biasing colour extraction; it is not used for shape or cover.

Full annotation definitions, output types, reference frames, and validation
metrics for every trait are in
[`data/ch1_image_traits/cirsium_ch1_image_trait_dictionary.csv`](../data/ch1_image_traits/cirsium_ch1_image_trait_dictionary.csv).
Summary of the five target traits and their operational definitions:

1. **Head orientation** — `head_orientation_class` (upright / inclined /
   nodding, thresholds at 30 deg and 60 deg from vertical) plus continuous
   `head_angle_deg` from a pedicel-to-apex keypoint vector, estimated in the
   image's gravity frame. Photos with a visibly tilted horizon/plant stem are
   flagged at annotation time and treated as a sensitivity subset, not
   silently included.
2. **Flower colour** — `flower_colour_class` (white / pale / dark) and
   `flower_colour_continuous` (CIELAB summary), both from ROI-B only, with an
   explicit caveat that citizen-science photos are not colour-calibrated
   across images (per-image white balance only; treat as relative/ordinal).
3. **Involucral cover / corolla exposure** — `involucral_cover_ratio`, the
   proportion of the head's apex-facing area that is involucre vs exposed
   corolla, from a 3-class segmentation (involucre / corolla / background) on
   ROI-A. Requires a view-angle qualifier (lateral vs frontal/oblique),
   because the same head shows different apparent cover depending on framing.
4. **Head shape** — `head_shape_aspect_ratio` (width:height on a head-axis
   rotated bounding shape, not the raw image axis) and a derived
   `head_shape_class` (globose / depressed / cylindrical) with thresholds
   fixed on a calibration subset before being applied to the full dataset.
5. **Relative corolla:involucre display** — `corolla_involucre_display_ratio`,
   computed from the same segmentation masks as (3) but asking a different
   question ("how much corolla visually dominates the head") rather than
   "how enclosed is the corolla." Its correlation with `involucral_cover_ratio`
   is reported explicitly rather than assumed (H1-Ch1d).

Traits intentionally excluded from Chapter 1 (mucilage, involucral spines,
leaf spines, pollination/herbivory/reproductive outcomes) are addressed in
§5.

## 3. AI pipeline and validation design

### 3.1 Detection and extraction stages

1. **Multi-class YOLO detector** (extends `models/best.pt`): capitulum,
   involucre, pedicel-base classes, optionally with a keypoint head
   (apex + pedicel-attachment points) for `head_angle_deg`. Requires a new
   annotation round beyond the current single-class capitulum boxes.
2. **ROI-A / ROI-B extraction**: standardized crop + padding from the
   multi-class detection (extends `v1/02_crop_heads_with_yolo.py`); ROI-B is
   produced by the segmentation model below rather than the detector alone.
3. **Segmentation model** (involucre / corolla / background) trained on
   manually annotated masks: feeds `involucral_cover_ratio`,
   `head_shape_aspect_ratio`, `head_shape_class`, and
   `corolla_involucre_display_ratio`.
4. **Per-trait extraction**:
   - Orientation: CNN classifier (extends
     `v1/03_train_head_direction_classifier.py`) + geometric angle estimate from
     keypoints.
   - Colour: deterministic CIELAB pixel statistics on ROI-B, with a CNN
     classifier fallback only if the deterministic rule underperforms human
     agreement.
   - Cover, shape, display: deterministic geometric/area computations on
     segmentation masks (no separate learned model beyond segmentation
     itself).

### 3.2 Human-label validation

Every extraction model is checked against independent human labels before
being trusted for the environmental analysis:

- Classification traits: accuracy, macro-F1, and Cohen's/weighted kappa
  against an adjudicated two-annotator label (disagreements resolved by a
  third annotator), not against a single annotator's opinion.
- Continuous traits: MAE and Spearman/Pearson correlation against a
  calibration subset with an independent measurement (protractor angle for
  orientation; ranking task for colour).
- Segmentation: Dice/IoU against manually drawn polygons on a labelled
  subset, plus correlation of the derived ratio traits against an
  independent human ordinal score (not just mask overlap, since the
  downstream quantity is the ratio, not the mask itself).

### 3.3 Group cross-validation and leave-one-species-out

Random image-level splits (the current approach in
`v1/03_train_head_direction_classifier.py`) are replaced with nested grouping,
run as separate, explicitly reported validation regimes rather than a single
number:

1. **Random split (baseline)** — kept only as a reference for how much the
   other regimes degrade performance, not reported as the headline result.
2. **Group by source photo/observation** — multiple crops or photos from the
   same iNaturalist observation must stay entirely in train or entirely in
   test, so the model cannot exploit repeated views of the same individual
   plant.
3. **Group by species** — GroupKFold over species, so within-fold species
   sets don't leak across train/test.
4. **Geographic block CV** — spatial blocking (e.g., grid-cell or buffered
   leave-block-out) so nearby, environmentally similar records aren't split
   across train/test.
5. **Leave-one-species-out (LOSO)** — for every species with at least a
   pre-declared minimum sample size, train on all other species and predict
   the held-out species; report the full per-species distribution of
   accuracy/correlation, not only the mean. Large per-species drops relative
   to (3) are treated as evidence the model may be reading species gestalt
   rather than the trait, and are discussed as a named limitation for that
   species/trait combination rather than averaged away.

### 3.4 Environmental and geographic analysis

Extends the existing `v1/05_rf_glm_gam_env_analysis.R` / `v1/05b_rf_only_env_analysis.R`
machinery (RF importance with repeats, correlation-cluster + VIF variable
selection, common-N weighted GLM, GAM) from a single binary response to:

- **Per-trait models** — one RF/GLM/GAM per trait (classification traits as
  binary/multinomial, continuous traits as Gaussian/Beta as appropriate)
  against CHELSA climate, topography, soil, and space, reusing the same
  variable-selection pipeline.
- **Trait-syndrome analysis** — PCA (continuous traits) or FAMD/MCA (mixed
  continuous + categorical) across the trait set at the observation level,
  followed by clustering (hierarchical or k-means on PCA/FAMD scores, or
  Gower-distance clustering) to define discrete syndrome types; syndrome
  axes/cluster membership are then modelled against the same environmental
  variable set as the per-trait models.
- **Sensitivity analysis across aggregation level** — every environmental
  model is re-run at (a) photo/observation level, (b) species-mean level, and
  (c) species x geographic-grid-cell level, to check that patterns are not
  artifacts of which species/regions iNaturalist happens to oversample.

## 4. Implementation roadmap (mapped to existing code)

Ordered by what blocks the rest of the pipeline, not by file number:

| Phase | Task | Relation to existing code |
|---|---|---|
| 0 | New annotation round: multi-class boxes (capitulum/involucre/pedicel), keypoints, involucre/corolla/background segmentation masks, and human trait labels (orientation, colour, cover, shape, display) on a shared calibration subset. | Blocks everything else; nothing downstream can be trained or validated without it. |
| 1 | Extend the detector to multi-class (+ optional keypoints). | Extends `v1/01_annotate_and_train_yolo.py`. |
| 2 | Build ROI-A/ROI-B extraction from the multi-class detector. | Extends `v1/02_crop_heads_with_yolo.py`. |
| 3 | Refactor the orientation classifier as one trait module and add the continuous angle estimator. | Reuses `v1/03_train_head_direction_classifier.py` / `v1/03b_..._legacy.py` and `v1/04_predict_head_direction_global.py` / `v1/04b_..._legacy.py` largely as-is; this is the lowest-risk phase since it already validates at Acc 0.896 on a random split — the new work here is regrouping the split (§3.3), not re-deriving the model. |
| 4 | Train the involucre/corolla/background segmentation model. | New; no existing counterpart. Highest modelling effort in the pipeline. |
| 5 | Colour extraction from ROI-B. | New, but mostly deterministic image processing; can proceed in parallel with Phase 4 once ROI-B exists in prototype form. |
| 6 | Shape and relative-display feature computation from segmentation masks. | New; depends on Phase 4 output. |
| 7 | Group-CV / LOSO validation harness applied to every trait model. | New; generalizes the ad hoc `validation_predictions.csv` / `validation_report.txt` pattern already used for orientation. |
| 8 | Extend `v1/05_rf_glm_gam_env_analysis.R` to per-trait models + add the trait-syndrome PCA/clustering script + the aggregation-level sensitivity script. | Extends `05_rf_glm_gam_env_analysis.R` / `05b_rf_only_env_analysis.R`; the variable-selection/RF/GLM/GAM machinery is reused, only the response side changes. |

Phase 3 can start immediately and ship an interim orientation-only result with
proper Group CV/LOSO validation while Phases 4-6 (the harder segmentation-
dependent traits) are still in development — this avoids blocking on the
hardest part of the pipeline before any validation-design improvement lands.

## 5. Boundary with Chapters 2-4 (avoiding overlap)

- **Mucilage** is not measured in Chapter 1 under any circumstance — image
  detection of mucilage is unreliable at citizen-science image quality. It
  remains exclusively a Chapter 3 field/near-macro/specimen trait
  (`floral_mucilage` in `cirsium_floral_trait_dictionary.csv`).
- **Leaf spines** are outside ROI-A/ROI-B entirely and are not introduced in
  Chapter 1. Combining leaf and floral defence requires an explicit
  cross-organ hypothesis, which Chapter 1 does not make (consistent with
  `ch3_trait_architecture/README.md`).
- **Pollination, floral herbivory/antagonism, and reproductive outcome
  effects** are reserved for Chapter 4's field manipulation experiments.
  Chapter 1 does not use interaction or fitness data even opportunistically
  (e.g., from iNaturalist comments); if cited at all, it is background
  motivation in the Introduction/Discussion, never an analysed variable.
- **Traits that share a name with a Chapter 3 field trait**
  (`involucral_cover`, `flower_colour`) are not the same measurement and
  should not be read as duplicating Chapter 3's contribution:
  - Chapter 1 measures them broadly (many species, uncontrolled
    citizen-science photography, image-derived, error-tolerant) to describe
    a **global macroecological pattern**, with no phylogenetic correction and
    no causal claim.
  - Chapter 3 remeasures the same concepts narrowly (fewer sites, controlled
    standardized photography/calipers/spectrometry) to test **trait
    syndrome structure after accounting for phylogeny**, per
    `floral_trait_architecture_protocol.md`.
  - `head_shape_aspect_ratio` (Chapter 1, relative, no scale bar) is likewise
    distinct from `capitulum_diameter` (Chapter 3, calibrated absolute mm);
    the two should never be merged into one variable in either chapter's
    dataset.
- **Chapter 2** (phylogeny, ancestral states, repeated evolution) consumes
  Chapter 1's trait table as input (species-level summaries) but performs no
  image-derived measurement itself.
- **`ch1_japan` (GBIF x GloBI pollinator SDM)** no longer represents part of
  "Chapter 1" under this redefinition, because it analyses pollinator
  assemblages and orientation together — exactly the interaction-effect
  scope now reserved for Chapter 4. It is recommended (not yet executed) to
  either (a) relabel it as a regional pilot feeding Chapter 4's pollinator
  rationale, or (b) fold its pollinator-SDM machinery into Chapter 4's field
  design as background/motivation material. This is a naming/placement
  decision for the author to confirm before the folder is moved.

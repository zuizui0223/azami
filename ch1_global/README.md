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
| H1-Ch1e | Capitulum traits are phenotypically integrated and some trait combinations are under-represented ("forbidden"), i.e. floral trait space is structured by coordination/constraint rather than freely combinatorial. | At the species level (§3.5): the trait correlation matrix departs from independence (integration index > null); trait space (PCA/FAMD) has empty or sparse regions relative to a null occupying the marginal ranges; putative modules (e.g. a display module {colour, size, exposure} vs a protection module {cover, orientation}) are recoverable. These become correlated-evolution hypotheses handed to Ch.2 (§6). |

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
4. Discussion — explicitly bounded: macroecological pattern only. Lay out the
   **candidate mechanisms** the pattern is consistent with — abiotic
   (rain/wind exposure), **pollinator-mediated selection**, floral-antagonist
   defence, and developmental/genetic constraint — as non-exclusive hypotheses
   *without testing them here*, and state what each later chapter contributes
   to testing them (Ch.2 phylogeny/correlated evolution, Ch.3 controlled field
   measurement, Ch.4 manipulation of pollination / antagonism / exposure). The
   pollinator hypothesis is named here as the explicit hook into Ch.4; no
   pollinator variable is analysed (§3.4, §5).
5. Limitations — citizen-science sampling bias, uncalibrated colour, camera
   tilt/viewing-angle confounds, ROI detection error propagating into
   downstream trait estimates.

### 1.4 Where the novelty comes from

Single-trait "nodding ~ climate" is a modest, largely known story and is not
the contribution. Chapter 1's novelty rests on three legs, in order of
strength:

1. **Multi-trait floral phenomics from citizen-science images, validated.** A
   reproducible pipeline that extracts a *suite* of capitulum traits (not one)
   from unstructured public photos and is validated with leave-one-species-out
   and human-agreement checks is a methodological contribution in its own
   right — most macro-scale floral-trait work relies on herbarium/literature
   scoring, not validated computer vision at this breadth.
2. **Global floral trait syndromes, integration, and constraint (§3.5).**
   Describing how capitulum traits co-vary, integrate into modules, and leave
   "forbidden" combinations empty, across a whole genus at global scale, is the
   biological contribution. It reframes the question from "where is nodding"
   to "how is the flowering head organized as a coordinated phenotype, and how
   does that organization map onto the planet."
3. **A phylogeny-free macroecological baseline for an intractable clade.**
   Precisely because *Cirsium* polyploidy/hybridization make a chloroplast
   phylogeny hard (see §3.4, §6), a rigorously space- and taxon-controlled
   *pattern* description is both necessary and novel: it is the honest first
   layer that a later, harder phylogenetic analysis (Ch.2) is built on, and it
   states its hypotheses of correlated evolution up front.

The framing to avoid: presenting the classifier accuracy or a single
climate–orientation curve as the result. The result is the trait-syndrome /
integration structure and its geographic-environmental correspondence, with
the pipeline as enabling method.

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
selection, common-N weighted GLM, GAM) from a single binary response to the
following. Note the design constraint that drives the structure of this
section: ***Cirsium* is polyploid with extensive hybridization and a
chloroplast phylogeny that is hard to resolve, so a reliable tree to correct
with may not exist even in Ch.2.** Chapter 1 therefore does not promise
"phylogenetic correction later" as its answer to non-independence (this
revises Audit finding B2); instead it owns a **phylogeny-free** strategy in
which space and taxonomic rank carry the load.

- **Non-independence without a reliable phylogeny.** Handle pseudoreplication
  with, in combination: (i) **spatial partitioning as a primary tool, not a
  robustness afterthought** — spatial block/cluster CV, spatial thinning, and
  a spatial term or spatial random effect (e.g. a smooth `s(lon,lat)` in the
  GAM, or an SPDE/CAR effect) so geography absorbs both sampling clustering and
  unmeasured spatially structured drivers; (ii) a **taxonomic random effect**
  (species, and section/subgenus where a classification exists) — a taxonomic
  hierarchy is a weaker substitute for a tree but still prevents one
  data-rich clade from dominating; (iii) **species and species×grid
  aggregation** as the headline unit (below); (iv) a **clade/section and
  single-species jackknife** — refit dropping each dominant species/section to
  show no one lineage produces the pattern. This replaces the phylogenetic
  control we cannot reliably build with structure we can.
- **Broadened predictor set beyond climate.** Because the strongest,
  most-novel signal is expected in trait *combinations* rather than in any one
  climate axis, and because climate alone risks a thin story, add non-climate
  layers to CHELSA + topography + soil: productivity/greenness (NDVI/EVI),
  **cloud- and fog-frequency and rain-day metrics** (directly relevant to the
  wet-exposure ideas around cover/orientation — e.g. EarthEnv cloud, CHELSA
  precipitation-timing variables), photoperiod/seasonality (latitude-derived),
  and land cover / human-footprint (which doubles as a sampling-bias covariate,
  §3.4 caveat). All enter the same correlation-cluster + VIF selection so
  collinear layers don't inflate the model.
- **Two-tier design (global + Japan) — with an explicit, reviewer-proof
  justification for the Japan tier.** The existing global/Japan code split
  becomes a deliberate two-tier structure: a **global tier** (all image traits
  × climate/terrain/soil/broadened layers × space, maximum breadth) and a
  **Japan deep-dive tier** (the same image traits at finer grain in a denser
  region). A reviewer *will* ask "why single out Japan when you already have a
  global analysis?" — and "I work there" / "more data" are not acceptable
  answers (Europe and North America are also densely sampled). The defensible
  reasons, which must be stated in the paper, are:
  1. **Japan is a *Cirsium* diversity hotspot** (many species, high endemism),
     so it lets us test trait–environment structure *within* a rich regional
     flora, partly decoupled from the broad between-region climate gradient
     that dominates the global tier — a distinct analytical question, not a
     smaller copy of the global one.
  2. **Scale-consistency test:** whether the global-scale trait–environment
     pattern replicates at fine grain in an independent, denser dataset is a
     genuine, phylogeny-free robustness argument; divergence flags a coarse-
     sampling artifact.
  3. **System continuity:** the downstream phylogenetic (Ch.2) and field
     (Ch.4) work is concentrated in Japan, so the Japan tier is the explicit
     bridge that sets up those chapters in the same system.

  Honest caveat: with the pollinator proxy deferred (below), the Japan tier's
  *distinctive* value rests on reasons 1 and 3, not on data density alone. If a
  reviewer still balks, the tier can be **demoted from a co-equal tier to a
  regional validation / sensitivity subsection** without loss to the global
  result; this fallback is deliberately kept open.
- **Pollinator hypothesis vs pollinator variable — keep the hypothesis, defer
  the variable.** Two different things must not be conflated:
  - The **pollinator *hypothesis*** (that pollinator-mediated selection is one
    candidate driver of the trait syndromes) *is* proposed in Chapter 1, as one
    of several non-exclusive candidate mechanisms in the framing and Discussion
    (§1.4, §5), and is the explicit hook into Chapter 4. This is the
    connectivity to later chapters and is wanted.
  - A **pollinator-availability *variable*** (a modelled proxy entered into the
    models) is **DEFERRED and not used in any Chapter 1 model.** Beyond the
    scope boundary, there is a decisive *methodological* reason: a stacked-SDM
    pollinator-richness surface is itself built from climate, so it is largely a
    nonlinear transform of predictors already in the model. It is therefore
    (i) collinear with climate and (ii) unable, even in principle, to separate a
    "pollinator effect" from a "climate effect" — precisely the separation the
    variable would be included to make. It also conflates total richness with
    the specific guilds that visit *Cirsium*, and availability with realized
    visitation. So the proxy cannot do the job the hypothesis needs; that job is
    Chapter 4's, with observed interaction data. If ever revisited, a
    functionally specific, regional (Japan-tier) availability surface from the
    `ch1_japan` GloBI + SDM machinery is the only honest form — raw GBIF
    pollinator density (confounded with human sampling effort) is never
    acceptable.
- **Per-trait models** — one RF/GLM/GAM per trait (classification traits as
  binary/multinomial, continuous as Gaussian/Beta), reusing the same
  variable-selection pipeline over the broadened predictor set.
- **Trait-syndrome models** — PCA (continuous) or FAMD/MCA (mixed) across the
  trait set, then clustering (hierarchical/k-means on scores, or Gower-distance
  clustering) to define discrete syndrome types; syndrome axes / cluster
  membership are modelled against the same predictors. Integration and
  constraint in that trait space are analysed in §3.5.
- **Aggregation-level sensitivity** — every model is re-run at (a)
  photo/observation level, (b) species-mean level, and (c) species×grid-cell
  level. Given the phylogeny problem, the **species-level (and species×grid)
  model is the headline** for any environmental claim; the per-observation
  model is shown only to demonstrate the pattern is not created by aggregation.
- **Inference philosophy (pre-declared).** With many traits × many predictors ×
  several tiers and aggregation levels, primary traits/hypotheses (§1.2) are
  declared confirmatory in advance and everything else is labelled
  exploratory; evidence is judged by effect size + CI and by **consistency
  across tiers and aggregation levels**, not by counting p-values (Audit D2).
- **Errors-in-variables.** Predicted trait values carry classifier/regression
  error into these models; carry predicted probabilities/continuous estimates
  (not hard labels) where possible, use a high-confidence subset as one
  sensitivity tier, and report how conclusions move (Audit D1).

### 3.5 Trait coordination, integration, and constraint

Beyond "do traits cluster into syndromes," Chapter 1 asks the sharper,
more-novel question the syndrome idea implies: **how coordinated and how
constrained is capitulum trait space?** This is the phenotypic-integration
framing and is computed at the **species level** (photo-level correlations are
inflated by within-individual pseudoreplication and by shared measurement
error, so they are not the unit for integration).

- **Trait covariation** — species-level correlation/partial-correlation matrix
  among all traits, with the measurement-error caveat (D1) stated; report which
  pairs covary beyond what shared size/allometry explains.
- **Integration magnitude** — an integration index (e.g. variance of the
  eigenvalues of the trait correlation matrix, or a relative eigenvalue
  standard deviation) tested against a null of independent traits; optionally
  whether integration strengthens in more stressful/seasonal climates.
- **Modularity** — test whether the capitulum behaves as one integrated module
  or as sub-modules (candidate a-priori split: a *display* module {colour,
  size, corolla exposure} vs a *protection/exposure* module {involucral cover,
  orientation}); report support for the modular hypothesis rather than assuming
  it.
- **Constraint / forbidden combinations** — quantify occupancy of trait space
  (convex-hull volume, kernel density, or empty-quadrant tests in PCA/FAMD
  space) against a null that fills the marginal ranges independently; sparsely
  occupied regions are candidate developmental/functional constraints.
- **Allometry** — trait–size scaling (e.g. colour or cover vs head size), since
  apparent integration can be size-driven and must be separated from it.

These outputs are not interpreted as adaptation in Chapter 1; each recovered
module, strong covariance, or forbidden combination becomes an explicit
**hypothesis of correlated evolution** passed to Chapter 2 (§6).

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
- **Pollinators: hypothesis in, variable out.** The pollinator *hypothesis* is
  deliberately proposed in Chapter 1 as one of several non-exclusive candidate
  mechanisms for the observed pattern (Discussion), because it is the explicit
  bridge to Chapter 4 — advancing it is *wanted* for chapter connectivity.
  What stays out of Chapter 1 is the **analysed pollinator data**: observed
  visits/contact/fitness (Chapter 4's field data) and any modelled
  pollinator-availability *variable* in the models (deferred; §3.4 gives the
  methodological reason — a stacked-SDM proxy is climate-derived and cannot
  separate pollinator from climate effects). So Chapter 1 *names* the
  hypothesis and hands it to Chapter 4 to *test*; it does not itself analyse a
  pollinator variable.
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
  Chapter 1's trait table and its §3.5 hypotheses as input; see §6.
- **`ch1_japan`** contributes the **Japan deep-dive tier** of Chapter 1
  (§3.4) as a denser-data region for a cross-tier consistency check — using its
  occurrence/environmental machinery, **not** its pollinator SDM, which is
  deferred with the pollinator proxy above. Its original pollinator-assemblage-
  vs-orientation interaction framing stays Chapter 4. No folder move is
  required; the same code serves both under clearly separated
  claims.

## 6. Handoff to Chapter 2 (and the phylogeny problem)

Chapter 1 is deliberately built so that its outputs are exactly Chapter 2's
inputs, and so that it does not over-promise given that a resolved tree may not
be available:

- **What Ch.1 hands to Ch.2:** a validated species-level trait matrix
  (orientation, colour, cover, shape, display + continuous companions), the
  syndrome/cluster definitions (§3.4), and — most importantly — the §3.5
  **integration/modularity/constraint results reframed as explicit hypotheses
  of correlated evolution** (e.g. "display and protection modules are
  decoupled," "cover and nodding co-occur beyond chance," "combination Z is
  forbidden"). These are the concrete propositions Ch.2 tests on a tree.
- **The phylogeny caveat, stated up front:** *Cirsium* polyploidy and
  hybridization make chloroplast phylogeny estimation hard, so Ch.2 may need
  multiple candidate trees, nuclear/target-capture data, or network rather
  than strictly bifurcating-tree methods, and some correlated-evolution tests
  may remain uncertain. Chapter 1 does **not** assume this will be solved: its
  own claims stand on space + taxonomic-rank control (§3.4) without requiring
  a tree. This makes Ch.1 a self-contained, phylogeny-free macroecological
  result and a scaffold for Ch.2, rather than a claim hostage to a phylogeny
  that may not resolve.
- **Bridge to Ch.3/Ch.4:** the syndromes and forbidden combinations that look
  environmentally structured in Ch.1 become the specific trait combinations
  Ch.3 remeasures under controlled conditions and Ch.4 manipulates in the
  field; §1.4 leg 2 is what those chapters turn from pattern into mechanism.

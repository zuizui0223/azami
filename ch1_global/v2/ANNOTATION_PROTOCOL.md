# Ch.1 v2 — annotation foundation and paper-lock protocol

## What this phase is for

This phase creates image-level biological ground truth before any global environmental model is re-run. It replaces the v1 practice of inheriting an orientation label from a species description. Each accepted image/head is labelled from the image itself, and annotators must not see taxon name, coordinates, model prediction, or previous annotator labels.

The unit for the trait label is a **detected capitulum in a source observation**. The unit for detector and segmentation annotation is the **full source image**, because a photo can contain more than one capitulum.

## The Ch.1 paper claim that this phase enables

> Across *Cirsium*, validated image-derived capitulum traits form repeatable phenotypic combinations, and the occupancy of those combinations varies across geographic and environmental space.

The paper does **not** claim that a trait evolved in response to climate, that a sparse region of trait space proves a developmental constraint, or that a pollinator variable explains the pattern. Those are Chapters 2–4 questions.

Use **phenotypic integration**, **under-represented / forbidden combinations**, and **environmental/geographic correspondence** as the Ch.1 language. Reserve evolutionary transition and causal-adaptation language for the phylogenetic and experimental chapters.

## Annotation tracks

| Track | Source shown to annotator | Labels/geometry | Purpose |
|---|---|---|---|
| Detector | full source photo | all visible flowering-head boxes | detection precision/recall and multi-head ROI-A source |
| Keypoints | ROI-A plus source photo | pedicel attachment + capitulum apex | continuous head angle |
| Trait labels | ROI-A; colour judgement additionally receives corolla-only crop when available | orientation, camera tilt, view angle, colour, cover, shape, display | biological ground truth and human-agreement ceiling |
| Segmentation | ROI-A | involucre, corolla, background polygons | cover, head shape, and display calculations |

Definitions are fixed in `data/ch1_image_traits/cirsium_ch1_image_trait_dictionary.csv`. Do not add or rename traits during annotation without a versioned dictionary update.

## Image eligibility and exclusion flags

Annotators make an eligibility decision before assigning a trait:

- `flowering_qc`: `usable`, `not_flowering`, `too_small_or_blurred`, `occluded`, `wrong_taxon_or_non-Cirsium`, `uncertain`.
- `camera_tilt_flag`: `vertical_reference_clear`, `vertical_reference_uncertain`, `tilted_or_rotated`, `not_assessable`.
- `view_angle_class`: `lateral`, `oblique`, `frontal`, `not_assessable`.

A photo can remain usable for colour or segmentation while being unusable for orientation. Cover and display analyses have a **lateral-view primary subset**; oblique/frontal images remain only in a declared sensitivity analysis. The main orientation subset includes only `vertical_reference_clear` images. This avoids treating the top edge of an arbitrary camera frame as biological gravity.

## Trait instructions

### Orientation

Label the head relative to the gravity direction inferred from a visibly upright plant stem or a reliable horizon. Use `upright` (0–30° from vertical), `inclined` (>30–60°), and `nodding` (>60°). When a reliable vertical reference is not visible, record `not_assessable`; do not infer orientation from species identity, flower posture stereotypes, or the photo frame alone.

### Colour

Label corolla pixels only as `white`, `pale`, or `dark`. This is an ordinal, image-relative assessment, not an absolute reflectance measurement. Mark `illumination_or_white_balance_unusable` when extreme shadow, colour cast, flash, or clipping makes the judgement unreliable.

### Cover, shape, and display

For a lateral view, score cover and display separately: cover asks how enclosed the corolla is by involucre; display asks how much corolla visually dominates. Shape concerns the head itself after mentally aligning it to its own long axis. Do not infer any of these from the taxon name.

## Annotation design and locked data partitions

1. **Calibration subset:** 720 observation-deduplicated heads from at least 36 species whenever coverage permits. Two independent annotators label every calibration item, then a third annotator adjudicates disagreements. This subset establishes human agreement and calibrates class definitions; it is never selected because a model found it easy.
2. **Main trait/detector pool:** target 3,000 observation-deduplicated heads. Cap contribution of common species, preserve coverage across spatial blocks, and record every exclusion.
3. **Segmentation subset:** target 1,000 heads drawn across species and spatial blocks. Mask annotators use source photo context when boundaries are ambiguous.
4. **Evaluation:** no random image-level result is headline evidence. Save `observation_group`, `species_cv_fold`, and `spatial_cv_fold` from the manifest. Report observation-group CV, species-group CV, spatial-block CV, and LOSO per species. A species, observation, or local spatial cluster must never cross a corresponding train/test boundary.

## Phase gates before global modelling

Do not expand to a global trait–environment result until the following are produced and inspected:

1. Detector precision/recall on held-out manually boxed source images, including errors by species and image quality.
2. Double-annotation agreement and adjudicated image-level labels for every retained trait.
3. A documented decision on each trait's usable image subset: orientation (vertical reference), colour (illumination), cover/display (view angle).
4. Grouped and LOSO performance, with the full species-level distribution, not only a pooled accuracy.
5. A versioned environmental-data extraction/provenance script and coordinate-QC rule before environmental inference.

## How Japan enters the first paper

Japan is a **regional replication and scale-consistency tier**, not a second paper hidden inside this one. The global tier discovers trait combinations and their broad correspondence with environment/geography. The Japan tier tests whether the same associations remain at finer scale in a species-rich study system that continues into the phylogenetic and field chapters.

The later Japanese field programme can then test the mechanisms that Ch.1 deliberately leaves open: rain/wind exposure, floral antagonism, pollinator behaviour, and reproductive consequences. The existing flower-colour, mucilage, and pollinator branches are therefore not lost; they become the mechanism-facing empirical stage rather than being forced into an observational macroecology paper.

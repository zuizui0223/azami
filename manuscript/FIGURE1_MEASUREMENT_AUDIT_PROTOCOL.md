# Figure 1 measurement and detector audit protocol

## Purpose

Figure 1 is a credibility display, not a new trait analysis. It must show how public photographs become inspectable continuous measurements while preserving detector failure, trait-specific unassessability, provenance and image licence restrictions.

## Required outputs

### Detector evaluation

Use an independent source-image subset that did not contribute to detector training or threshold selection. The audit set must contain:

- images with zero, one and multiple visible capitula;
- distant, partly occluded and complex-background heads;
- buds, anthesis-stage and post-anthesis heads where the object definition permits them;
- broad taxonomic and geographic coverage;
- explicit unassessable source images rather than silently excluding them.

Each visible target is boxed once. Images with no target retain one row with an empty `object_id` and empty coordinates. Run:

```bash
python ch1_global/v2/79_evaluate_capitulum_detector.py \
  --annotations detector_audit_annotations.csv \
  --predictions yolo_crop_metadata.csv \
  --out-dir detector_audit \
  --iou-threshold 0.5
```

Report precision, recall, F1, object counts and image counts. Supplementary material must also show missed-head and false-positive categories. Trait QC cannot diagnose a capitulum that the detector missed.

## Detector annotation schema

Required columns:

- `image_id`: immutable queue/source-image identifier matching the prediction table;
- `object_id`: unique target identifier within image; empty for a true no-head image;
- `assessable_for_detection`: whether a human can reliably audit all target objects;
- `x1`, `y1`, `x2`, `y2`: source-pixel coordinates;
- `missed_head_category`: optional manual category used for false negatives;
- `audit_notes`: optional free text.

Recommended missed-head categories are `small_or_distant`, `occluded`, `edge_truncated`, `bud_or_postanthesis`, `unusual_view`, `crowded_heads`, and `other`.

## Figure 1 panel selection

Selection is manual but reproducible. It must include:

1. source image with detector bounding box;
2. tight crop;
3. wider context crop;
4. orientation overlay;
5. colour-pixel overlay;
6. shape contour, convex hull and major-axis overlay;
7. representative assessable and trait-specific unassessable cases;
8. at least one horizontal-flip technical replicate in the assembled figure or supplement.

Run:

```bash
python ch1_global/v2/78_build_figure1_measurement_audit.py \
  --selection figure1_panel_selection.csv \
  --measurements primary_trait_continuous_head_measurements.csv \
  --out-dir figure1_measurement_audit
```

The overlay script imports the frozen v2 production measurement functions. Figure-only substitute segmentation is not allowed.

## Trait-specific assessability

Do not label a whole photograph simply assessable or unassessable. Record separate states for:

- orientation;
- colour;
- shape.

One image may support colour and shape while failing orientation. `unassessable` means that the photograph does not support the measurement; it is never biological absence.

Recommended panel roles:

- `all_traits_assessable`;
- `orientation_unassessable`;
- `colour_unassessable`;
- `shape_unassessable`;
- `occlusion_failure`;
- `background_leakage_failure`;
- `camera_tilt_example`;
- `horizontal_flip_replicate`;
- `low_value`, `middle_value`, `high_value` for each continuous endpoint.

## Licence gate

Source photographs may be analysed without redistribution, but a manuscript figure republishes an image. Every selected panel must retain:

- observation and photo identifiers;
- taxon;
- photographer attribution;
- exact photo licence code;
- original observation URL;
- access date;
- whether cropping or overlays were applied;
- explicit `licence_verified=true`.

The default overlay builder accepts CC0, CC BY and CC BY-SA variants. ND-licensed images must not be cropped or overlaid. NC and unrecognised licences require explicit publisher/legal review and must not bypass the default gate casually.

## Figure/Supplement division

Main Figure 1 should remain compact:

- workflow;
- one source-to-crop example;
- one overlay example per module;
- a small set of assessable/unassessable examples;
- the transition from repeated observations to a species distribution.

Place full detector metrics, error categories, low/median/high endpoint examples, horizontal-flip pairs and additional QC failures in Supplementary figures.

## Freeze boundary

This work does not change detector weights, trait thresholds or frozen numerical results. If the audit demonstrates a material detector or measurement failure, create a new analysis version rather than overwriting the submitted result.

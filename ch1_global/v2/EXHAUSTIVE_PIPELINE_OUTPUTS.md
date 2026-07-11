# Exhaustive all-photo pipeline outputs

## Source-of-truth outputs

- `exhaustive_photo_queue_merged.csv` — every queued public photo;
- `exhaustive_screening_results.csv` — one terminal YOLO status per photo;
- `exhaustive_photo_detection_audit.csv` — photo metadata plus inclusion/exclusion reason;
- `exhaustive_species_detection_audit.csv` — species-level head-photo availability;
- `exhaustive_yolo_crop_metadata.csv` — detector-positive head boxes;
- `exhaustive_continuous_head_level.csv` — corrected continuous measurements;
- `exhaustive_continuous_photo_level.csv` — usable-head medians by photo;
- `exhaustive_continuous_observation_level.csv` — medians across photos by observation.

Leaf-only and other no-head photographs appear in the first four outputs with
`screen_status=no_detection`. They never appear in the continuous-trait tables.

## Derived cohorts

- `all_detected_observations.csv`;
- `coordinate_usable_observations.csv`;
- `strict_10km_observations.csv`;
- `strict_spatial_thinned_observations.csv`;
- `between_species_balanced_40_observations.csv`.

The derived cohorts are disposable analysis views. They never replace the
source-of-truth outputs.

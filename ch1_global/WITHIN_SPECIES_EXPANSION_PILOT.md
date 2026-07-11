# Exhaustive all-photo within-species expansion

## Decision

The 40-observation cap belongs only to the balanced between-species atlas. It
must not be applied before machine prediction when evaluating within-species
variation.

The high-density pipeline therefore follows this order:

```text
all eligible Cirsium photos
        |
YOLO visible-capitulum detection on every photo
        |
leaf-only / no-head photos retained as no_detection audit rows
        |
corrected continuous colour, orientation and outline measurement
        |
photo -> observation aggregation
        |
coordinate-quality and spatial sensitivity cohorts
        |
optional 40-observation balanced derivative for between-species comparison
```

The unthinned prediction tables are the source of truth. Spatial thinning never
overwrites them.

## Pre-prediction inclusion

The exhaustive queue applies only filters needed to define the intended public
photo population:

- species-level *Cirsium* identification;
- selected iNaturalist quality grade, normally research grade;
- non-captive observation;
- retrievable photo URL;
- optional licence or coordinate policy.

It does **not** apply:

- a maximum number of photos or observations per species;
- one-photo-per-observation selection;
- spatial or environmental balancing;
- flower/head annotations;
- any measured or predicted trait value.

Every photo from an eligible observation is queued. Multiple photographs from
one observation are technical views of the same observation and are aggregated
after head measurement.

## Leaf-only and no-head photographs

All photos are screened by the frozen YOLO detector. Photos with no detected
visible capitulum, including leaf-only images, are:

- labelled `no_detection`;
- excluded from colour, orientation and outline tables;
- retained in `exhaustive_photo_detection_audit.csv`;
- summarized by species to quantify differential availability of head-bearing
  photographs.

Thus a leaf-only photo cannot receive a floral trait value, but it also does not
vanish from the sampling audit.

## Correct analysis units

- detector box: one visible-capitulum candidate;
- photo: median across usable detected heads in one photograph;
- observation: median across all detector-positive photos for one `obs_id`;
- species: derived only after observation aggregation.

Multiple photos and multiple detected heads do not become independent plants.

## Post-prediction cohorts

`74_build_postprediction_analysis_cohorts.py` creates explicit derived tables:

1. `all_detected_observations.csv` — every observation with at least one detected
   head;
2. `coordinate_usable_observations.csv` — public coordinates usable for
   environment extraction;
3. `strict_10km_observations.csv` — positional accuracy <=10 km;
4. `strict_spatial_thinned_observations.csv` — one observation per species and
   configurable spatial cell;
5. `between_species_balanced_40_observations.csv` — a separate balanced
   derivative, never the within-species source.

Within-species models should compare the unthinned, strict-coordinate and
spatially thinned cohorts. The 40-observation derivative exists only for a fair
comparison with the original between-species atlas.

## Workflows

- `ch1-build-exhaustive-photo-queue.yml` builds the full queue and partitions
  every photo exactly once across 1–100 shards.
- `ch1-run-merge-exhaustive-continuous.yml` runs detector and corrected
  continuous measurements, merges all shards, preserves leaf/no-head audit rows
  and builds post-prediction cohorts.

Large runs are manual and resumable by shard. Result artifacts contain metadata,
screening status and numeric traits, not redistributed source photographs.

## Interpretation

This expansion can materially improve within-species resolution for common
species. It does not automatically remove:

- observer and road-access bias;
- non-random flowering photographs;
- detector false negatives;
- image-angle and colour-calibration error;
- spatial autocorrelation;
- uneven climatic ranges among species.

Those limitations are evaluated after prediction rather than used as a reason
to discard most photos beforehand.

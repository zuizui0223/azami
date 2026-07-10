# Continuous automated primary traits

## Decision

The Chapter 1 image analysis does not require a human to assign every head to
one discrete category. The main automated measurements are continuous:

- **orientation** — signed longitudinal head-axis angle relative to the
  EXIF-oriented image vertical: 0° upward, 90° horizontal and 180° downward;
- **corolla colour** — visible-floret fraction, Lab lightness and chroma,
  circular hue coordinates, and mutually exclusive colour-family fractions;
- **capitulum outline** — aspect ratio, circularity, solidity and width-profile
  variables.

The earlier CLIP classes remain useful as exploratory descriptions, but they do
not define ground truth and are not required by this pipeline.

## Version 2 corrections

The first complete run was retained as a diagnostic and revealed two issues
before ecological modelling:

1. white/cream and yellow masks could overlap, allowing visible-corolla coverage
   to exceed one;
2. the original orientation value was a head–peduncle bend angle and therefore
   could not reliably distinguish an upward head from a downward head.

Production runs now pass through `56_run_primary_traits_continuous_v2.py`.
Colour masks are mutually exclusive and checked to sum to one. The sign of the
head PCA axis is selected from the visible floral-pixel centroid, after which
the angle is measured against image vertical. `55_run_primary_traits_continuous.py`
remains the workflow entrypoint but routes execution to v2.

## Why continuous measurements

Boundaries such as upright versus ascending, globose versus hemispherical, or
white versus pale pink are partly conventional. Asking one person or one
zero-shot model to force every photograph into these bins creates avoidable
measurement error. Continuous image measurements preserve more information and
allow uncertainty to be propagated into the ecological models.

## Image processing

The pipeline starts from the existing YOLO-detected head and context crops. It:

1. estimates the central foreground from colour contrast with the crop border;
2. extracts mutually exclusive colour masks and outline measurements;
3. estimates the undirected longitudinal head axis by PCA;
4. uses the visible floral-pixel centroid to select the floral end of that axis;
5. calculates the signed axis angle relative to image vertical;
6. repeats every measurement after horizontal mirroring;
7. marks a trait `usable` only when image quality and flip repeatability pass the
   predeclared thresholds.

The flip is a technical replicate. A biological measurement should not change
when the photograph is mirrored.

## Orientation interpretation

`orientation_angle_degrees` is image-referenced rather than a direct
inclinometer measurement of gravity. iNaturalist images are read after EXIF
orientation, so image vertical is a reproducible proxy, but tilted photographs
remain a source of error. Analyses must retain only `orientation_status=usable`,
report this limitation, and use species/observation aggregation to reduce
single-image noise.

## Analysis units

`53_aggregate_primary_trait_continuous.py` calculates medians and median absolute
deviations:

- detected head → observation;
- observation → species.

Only rows with the corresponding `*_status=usable` contribute to a trait
summary. Numbers of usable heads and observations are retained rather than
silently treating missing measurements as biological absence.

## Recommended model variables

- `orientation_angle_degrees_median`
- `corolla_lab_lightness_median`
- `corolla_lab_chroma_median`
- `corolla_hue_sin_median` and `corolla_hue_cos_median`
- `shape_aspect_ratio_median`
- `shape_circularity_median`
- `shape_solidity_median`
- `shape_width_cv_median`

Hue must be modelled circularly through sine and cosine, not as a linear
0–360-degree variable.

## Categorical figures

Categories may still be created for maps or descriptive plots, but thresholds
must be declared after the continuous measurement system is frozen. Those bins
are display conventions and must not be described as independently learned
biological truth.

## Human audit

The completed human audit is retained as a diagnostic of ambiguity. It is not
used to fit colour, shape or orientation thresholds. Therefore disagreement
between the annotator and the old CLIP classifier does not become training noise
in the automated measurements.

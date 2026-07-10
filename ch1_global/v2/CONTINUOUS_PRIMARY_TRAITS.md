# Continuous automated primary traits

## Decision

The Chapter 1 image analysis no longer requires a human to assign every head to
one discrete category. The main automated measurements are continuous:

- **orientation** — angle between the recovered peduncle direction and the
  longitudinal head axis, from 0 to 180 degrees;
- **corolla colour** — visible-floret fraction, Lab lightness and chroma, circular
  hue coordinates, and colour-family pixel fractions;
- **capitulum outline** — aspect ratio, circularity, solidity and width-profile
  variables.

The earlier CLIP classes remain useful as exploratory descriptions, but they do
not define ground truth and are not required by this pipeline.

## Why

Boundaries such as upright versus ascending, globose versus hemispherical, or
white versus pale pink are partly conventional. Asking one person or one
zero-shot model to force every photograph into these bins creates avoidable
measurement error. Continuous image measurements preserve more information and
allow uncertainty to be propagated into the ecological models.

## Image processing

`52_measure_primary_traits_continuous.py` starts from the existing YOLO-detected
head and context crops. It:

1. estimates the central foreground from colour contrast with the crop border;
2. extracts colour and outline measurements from the head crop;
3. locates the head crop inside the context crop and recovers a nearby peduncle
   line for the orientation angle;
4. repeats every measurement after a horizontal flip;
5. marks a trait `usable` only when image quality and flip repeatability pass the
   predeclared thresholds.

The flip is a technical replicate. A biological measurement should not change
when the photograph is mirrored.

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
used to fit colour, shape or orientation thresholds in this continuous pipeline.
Therefore disagreement between the annotator and the old CLIP classifier does
not become training noise in the new measurements.

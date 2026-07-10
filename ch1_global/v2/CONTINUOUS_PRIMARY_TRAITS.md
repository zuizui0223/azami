# Continuous automated primary traits

## Decision

Chapter 1 does not require a human or zero-shot model to force every head into a
discrete category. The primary automated measurements are continuous:

- **orientation** — signed longitudinal head-axis angle relative to the
  EXIF-oriented image vertical: 0° upward, 90° horizontal and 180° downward;
- **corolla colour** — visible-floret fraction, Lab lightness and chroma,
  circular hue coordinates, and mutually exclusive colour-family fractions;
- **capitulum outline** — aspect ratio, circularity, solidity and width-profile
  variables.

Earlier CLIP classes remain diagnostic descriptions only. They do not define
ground truth and are not used by the submission pipeline.

## Production version

The first complete continuous run revealed two issues before ecological
modelling:

1. white/cream and yellow masks could overlap, allowing visible-corolla coverage
   to exceed one;
2. the original orientation value was a head–peduncle bend angle and could not
   reliably distinguish upward from downward heads.

The frozen production entrypoint is
`56_run_primary_traits_continuous_v2.py`. Colour masks are mutually exclusive and
checked to sum to one. The sign of the head PCA axis is selected from the
visible floral-pixel centroid, after which the angle is measured against image
vertical.

`52_measure_primary_traits_continuous.py` supplies shared image-processing
functions. `55_run_primary_traits_continuous.py` is retained only as the tested
OpenCV compatibility layer imported by the frozen v2 implementation. New
workflows must invoke script 56 directly.

## Why continuous measurements

Boundaries such as upright versus ascending, globose versus hemispherical, or
white versus pale pink are partly conventional. Forced bins introduce avoidable
measurement error. Continuous image measurements preserve more information and
make QC and uncertainty explicit.

## Image processing

The pipeline starts from YOLO-detected head and context crops. It:

1. estimates the central foreground from colour contrast with the crop border;
2. extracts mutually exclusive colour masks and outline measurements;
3. estimates the undirected longitudinal head axis by PCA;
4. uses the visible floral-pixel centroid to select the floral end of that axis;
5. calculates the signed axis angle relative to image vertical;
6. repeats every measurement after horizontal mirroring;
7. marks a trait `usable` only when image quality and flip repeatability pass the
   frozen thresholds.

The mirrored image is a technical replicate. A biological measurement should
not change when the photograph is mirrored.

## Orientation interpretation

`orientation_angle_degrees` is image-referenced rather than a direct
inclinometer measurement of gravity. Images are read after EXIF orientation, so
image vertical is a reproducible proxy, but tilted photographs remain a source
of error. Analyses retain only `orientation_status=usable`, report this
limitation and aggregate across heads and observations to reduce single-image
noise.

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

Hue is circular and must be represented jointly by sine and cosine. A
significant sine or cosine coefficient alone is not a complete biological colour
interpretation.

## Categorical figures

Categories may be created for maps or descriptive figures only. Thresholds must
be declared after the continuous measurement system is frozen. Display bins are
not independently learned biological states.

## Human-category audit

The earlier human-category exercise is retained in project history as a
diagnostic of ambiguity. It is not used to fit colour, shape or orientation
thresholds. Disagreement between a person and the old CLIP classifier therefore
does not become training noise in the continuous measurements.

# Continuous trait universe redesign

## Decision

Chapter 1 is not an involucre-only extension and is not a catalogue of visual
categories. Its measurement claim is that heterogeneous public photographs can
yield repeated **continuous, trait-specific and explicitly assessable numerical
phenotypes**. This contract applies equally to orientation, colour, display,
outline, involucral architecture, armature and high-resolution surface signals.

The botanical categories in `ch1_trait_ontology.csv` remain annotation and
calibration language. They are forbidden as inferential outcomes in the
continuous pipeline.

## Problem fixed

The previous repository contained three incompatible layers:

1. categorical CLIP candidates;
2. nine whole-capitulum numerical model columns;
3. several high-resolution contour columns, including two differently named
   columns with the same maximum-projection formula.

`ontology/ch1_continuous_trait_contract.csv` is now the only registry for
continuous inference. It distinguishes a biological construct from the actual
image measurement, gives units and bounds, records image requirements and states
what each number must not be interpreted as. Formula identifiers must be unique,
so a duplicated statistic cannot masquerade as another trait.

## Why hundreds of thousands of photographs do not equal one analysis sample

The denominators refer to different units and must not be collapsed into one
three-line explanation. The frozen exhaustive queue contains 777,766 photographs
from 460,036 observations. Detection retained 406,582 observations with at least
one visible capitulum. Coordinate availability retained 392,989 observations;
the predeclared positional-accuracy threshold retained 297,293; and one
observation per taxon × 0.25-degree cell retained 46,276 observations from 259
taxa. This last reduction controls spatial pseudoreplication and occurs only
after detection and continuous prediction; it is not a trait-based choice.

Extended endpoints then have separate resolution, view, segmentation, focus and
mirror-repeatability gates. Those gates create different usable counts for each
endpoint, not one universal "usable photograph" count. The workflow therefore
exports the full attrition chain and endpoint-specific coverage. The old
3,725-photograph atlas remains a balanced development and variance-partitioning
resource; it is not the sampling frame for the exhaustive extended analysis.

## Measurement modules

### Orientation

The current numerical endpoint is the signed PCA head axis relative to
EXIF-oriented image vertical. It is retained because it is reproducible, but the
contract explicitly prevents it from being described as an inclinometer angle or
head-peduncle deflection. A future local-peduncle reference requires continuous
keypoint validation and must be added as a separate endpoint rather than silently
changing this one.

### Colour

Lightness, chroma and circular hue remain continuous. Hue is one biological
dimension represented by a sine/cosine pair and receives one joint permutation
test. Mutually exclusive colour-family pixel fractions are a compositional
description and cannot enter separate univariate primary models.

### Floral display

Visible floret fraction is a continuous image phenotype. It is a candidate
display/flowering-exposure measurement, not a direct anthesis score or a
whole-plant display measure until calibrated against continuous stage annotation.

### Capitulum outline

Aspect ratio, circularity, solidity and longitudinal width-profile variation
remain continuous primary endpoints. They quantify the visible two-dimensional
outline and are not three-dimensional head shape.

### Involucral architecture and armature

The extended high-resolution measurement now returns:

- involucre length/width ratio;
- apical and basal taper ratios;
- radial projection roughness, 95th percentile and maximum;
- fraction of outward-projecting radial bins;
- projection-peak density;
- bilateral projection asymmetry.

The old `spine_relative_length_max_proxy` alias is not a separate endpoint. The
canonical value is `involucre_projection_max`. These quantities remain
contour-derived until bract base/tip reference annotations establish their
relationship to free-tip length, spreading angle and recurvature.

### Indumentum and visible exudate signal

When the detected head is at least 300 pixels across and passes focus, mask and
mirror-repeatability gates, the pipeline additionally records edge density, LBP
texture entropy, normalized high-frequency energy and a low-saturation
high-value specular fraction. These are validation-only numerical image signals.
They are not automatically named hair density, gland density or mucilage.

## Analysis tiers

- `primary`: frozen whole-capitulum endpoints with technical repeatability.
- `candidate`: continuous endpoints that may enter a separately corrected
  exploratory family after coverage checks.
- `validation_only`: numeric surface measurements requiring botanical reference
  calibration before trait-environment inference.
- `descriptive_only`: constrained compositions or support quantities that must
  not be treated as independent univariate traits.

Missing or unstable measurements remain missing. No tier converts
`unassessable` to absence or zero.

## Required execution order

1. Validate the contract with `87_validate_continuous_trait_contract.py`.
2. Run the frozen v2 whole-capitulum measurements on the exhaustive detector
   output.
3. Intersect exhaustive crop metadata with the frozen strict spatial cohort,
   apply only predeclared resolution and URL gates, and select one deterministic
   highest-resolution head per observation. Trait values and category labels are
   forbidden from this selection.
4. Run `89_measure_extended_continuous_traits.py` on those crops.
5. Build the category-free long tables with
   `88_build_continuous_trait_universe.py`.
6. Join the same spatially thinned observation cohort and environment table.
7. Run `90_run_continuous_trait_universe_climate.py`. Primary and candidate
   families receive separate BH correction; circular hue receives one joint
   within-taxon permutation test.
8. Promote a candidate only after double-annotated continuous reference data
   validate direction, scale and error against viewpoint and resolution.

## Continuous reference annotation

The biological calibration packet must sample low, middle and high values of
every candidate metric plus failures, across taxa, image resolutions and views.
Two annotators record continuous references rather than whole-head classes:

- head-axis and visible peduncle vectors;
- involucre polygon and longitudinal axis;
- bract base/tip keypoints for visible outer, middle and inner bracts;
- free-tip length divided by involucre width;
- bract angle relative to the local involucre surface;
- signed curvature along the free bract tip;
- visible hair/filament area and gland-point count per annotated surface area;
- visible exudate area, with a separate `not visible enough` flag.

Required validation summaries are continuous agreement (ICC or repeatability,
rank correlation, MAE and calibration slope), mirror repeatability, resolution
strata, viewpoint strata and taxon-held-out performance. Category accuracy alone
cannot validate this pipeline.

## Manuscript claim boundary

The defensible novelty is the end-to-end measurement framework: repeated public
photographs are retained as continuous image phenotypes with trait-specific
assessability, not reduced to species categories. The present data can establish
technical reproducibility and ecological association. Botanical identity,
plasticity, adaptation, defence and secretion mechanisms require the reference
and experimental layers stated above.

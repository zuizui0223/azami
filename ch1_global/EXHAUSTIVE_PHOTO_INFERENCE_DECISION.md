# Exhaustive photo inference before thinning

Decision date: 2026-07-11

## Correction to the original sampling design

The existing 3,725-observation dataset was intentionally capped for balanced
between-species comparison. That cap is not an appropriate source dataset for
high-resolution within-species inference.

Chapter 1 will therefore retain two distinct products:

1. **Exhaustive prediction layer** — every eligible public *Cirsium* photo is
   passed through the frozen head detector, and every detector-positive head is
   measured with the corrected continuous-trait pipeline.
2. **Derived analysis cohorts** — coordinate-quality, spatially thinned and
   species-balanced tables created only after all predictions are archived.

No thinned table may replace or overwrite the exhaustive source predictions.

## Leaf-only photographs

Leaf-only and other no-head photographs are expected in the public-image pool.
They are processed by YOLO and retained as `no_detection` rows. They do not enter
colour, orientation or capitulum-outline measurement.

This creates two useful outcomes:

- floral traits cannot be assigned to a photograph without a detected visible
  capitulum;
- the probability that an eligible public photograph contains a detectable
  capitulum can be audited by species, region and environment.

## Revised role of the current 3,725 observations

The current balanced dataset remains valid for the existing between-species
atlas and historical-sensitivity analysis. Its within-species null results are
provisional because the data were not built to maximize within-species density.

The exhaustive layer will determine whether high-density within-species analysis
is feasible for each trait and species. The final manuscript wording must not
claim absence of within-species effects until this expanded analysis is complete.

## Analysis-unit hierarchy

- head: one YOLO detector box;
- photo: median across usable heads in one photo;
- observation: median across all detector-positive photos sharing one iNaturalist
  observation ID;
- species: summary across observations.

Multiple photos and multiple heads remain repeated views, not independent
plants.

## Post-prediction cohorts

The pipeline creates:

- all detector-positive observations;
- coordinate-usable observations;
- positional-accuracy <=10 km observations;
- one observation per species and 0.25-degree spatial cell;
- a maximum-40 species-balanced derivative for direct comparison with the
  original atlas.

The first four are available for within-species modelling. The fifth is for
between-species balancing only.

## Remaining bias controls

Exhaustive prediction improves sample size but does not eliminate:

- observer and road-access bias;
- seasonal and flowering-photo bias;
- detector false negatives;
- repeated observations of the same local population;
- spatial autocorrelation;
- camera and viewpoint effects.

These must be modelled or tested after prediction through detection audits,
spatial cohorts, observation aggregation and regional sensitivity analyses.

# Main Figure 2 geography release manifest

This manifest defines the submission-facing geographic visualization layer for Chapter 1. It is a visualization and sampling-bias audit only; it does not create a new inferential family or alter any frozen biological result.

## Canonical products

- **Main Figure 2** `Figure_2_geographic_sampling_and_analysis_domain.*`
  - A: 1° density of detector-positive observations before spatial filtering. The frozen detector-positive cohort contains 406,582 observations; 406,178 have finite longitude/latitude values that can be rendered geographically, while the full 406,582 count remains the cohort total.
  - B: coordinates of the 46,276-observation spatially thinned primary cohort.
  - C: the two executed analysis streams and their derived cohorts. The balanced image atlas leads to nested variance/PCA/historical analyses and the high-resolution involucre subset; the exhaustive detector-positive stream passes through coordinate and positional-accuracy filtering and taxon × 0.25° thinning before the primary and grouped SPDE analyses.
  - D: CHELSA BIO1 × BIO12 environmental domain of the primary cohort.
- **Figure S6** `Figure_S6_sampling_geography_across_filters.*`
  - 1° sampling density for detector-positive, coordinate-usable, ≤10 km and primary-thinned stages on a common log-count scale. The first panel labels both the frozen detector-positive total and the number of mappable coordinates.
- **Figure S7** `Figure_S7_geographic_trait_assessability.*`
  - 2° cell-level usable fractions for orientation, colour, outline and the mean across all nine primary endpoints; cells require at least 20 coordinate-usable observations.

Each figure is released as SVG, PNG and PDF. The release artifact also contains the aggregated plotting tables, a summary JSON and SHA-256 inventories.

## Frozen source runs

### Exhaustive post-prediction cohorts

- workflow run: `29216617585`
- artifact ID: `8269246732`
- artifact name: `ch1-exhaustive-continuous-merged-all_photos_20260711_v1-recovery`
- artifact digest: `sha256:5f18b42d18cfcb81691c38ce0f04bcef754e6a67382025ea90110dbc50ae194b`
- source tables used:
  - `all_detected_observations.csv`
  - `coordinate_usable_observations.csv`
  - `strict_10km_observations.csv`
  - `strict_spatial_thinned_observations.csv`

Frozen cohort checks:

| Stage | Observations | Taxa |
|---|---:|---:|
| Detector-positive | 406,582 | 286 |
| Coordinate usable | 392,989 | 271 |
| ≤10 km positional accuracy | 297,293 | 259 |
| Primary spatially thinned | 46,276 | 259 |

The geographic generator additionally records the number of rows with finite coordinates at every stage. This rendering check is separate from the frozen cohort definition and never changes cohort membership.

### Enriched primary environment

- workflow run: `29306454759`
- artifact ID: `8314947270`
- artifact name: `ch1-exhaustive-enriched-environment-20260713`
- artifact digest: `sha256:b1252258fff85b397de8a58f322e24e90731d55f07f53ff7de16a169c5a4dfe3`
- source table: `strict_spatial_thinned_with_climate_topography_soil.csv`
- rows: 46,276; taxa: 259.

CHELSA BIO1 is converted from the frozen raster storage transform (`scale = 0.1`, `offset = -273.15`) to degrees Celsius for display. BIO12 is displayed in the frozen precipitation units (mm).

## Map geometry

Natural Earth admin-0 50 m country polygons are downloaded during the workflow only as geographic context. They are not used for any biological inference.

## Reproduction

- generator: `analysis/build_main_figure2_geography.py`
- workflow: `.github/workflows/ch1-main-figure2-geography-ci.yml`

The workflow downloads the frozen source artifacts directly, validates exact cohort counts, generates the three figures and their aggregated plotting data, records source and release hashes, and uploads one release artifact. Large observation-level source tables are not duplicated in Git.

## Interpretation boundary

The maps visualize where public-image observations occur, how filtering changes spatial representation, where the primary analytical cohort lies in geographic/environmental space, and whether image-trait assessability varies geographically. They do not demonstrate sampling representativeness of all global *Cirsium*, causal environmental effects, adaptation, or absence of residual observation bias. Figure S7 shows measurement availability in images; unassessable images are not biological absences.

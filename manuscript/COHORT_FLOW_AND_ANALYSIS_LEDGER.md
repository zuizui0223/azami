# Chapter 1 cohort flow and analysis ledger

This ledger replaces the ambiguous use of *strict*, *expanded*, *balanced* and *primary*. Each label below refers to one immutable table and one analysis purpose.

## Exhaustive all-photo inference stream

The exhaustive stream began with all eligible public photographs before any trait-based thinning. Cohorts were then derived using coordinate quality, positional accuracy, spatial cells and stable hashes only.

| Canonical cohort name | Frozen file | Observations | Taxa | Selection | Analyses permitted |
|---|---|---:|---:|---|---|
| Exhaustive detected | `all_detected_observations.csv` | 406,582 | 286 | At least one detector-positive capitulum | Detection and availability accounting only |
| Exhaustive coordinate-usable | `coordinate_usable_observations.csv` | 392,989 | 271 | Public usable coordinates | Geographic coverage summaries |
| Exhaustive ≤10 km | `strict_10km_observations.csv` | 297,293 | 259 | Positional accuracy ≤10 km | Pre-thinning sensitivity only |
| **Exhaustive spatially thinned primary** | `strict_spatial_thinned_observations.csv` | **46,276** | **259** | One observation per taxon × 0.25° cell after ≤10 km restriction | Primary within-species climate coefficients and revised precision-aware lability analysis |
| Exhaustive between-species balance | `between_species_balanced_40_observations.csv` | 3,723 | 259 | Stable-hash maximum 40 observations per taxon | Exhaustive-stream between-species sensitivity only |

The exhaustive merge contained 777,766 photographs from 460,036 source observations, 637,745 detector-positive photographs, 406,582 observations with detected heads and 1,255,791 detected heads. These counts describe the execution stream, not the smaller image-comparison atlas.

## Balanced image-comparison atlas

A separately executed balanced image layer contained 6,626 detected capitula from 3,725 observations and 216 accepted image-analysis taxa. It supports visible variance partitioning, species-level trait PCA, among-species trait summaries and historical-placement sensitivities.

Its earlier observation-level climate screens are retained as sensitivity analyses:

| Sensitivity analysis | Endpoint rows | Typical usable rows per endpoint | BH-supported main endpoint–predictor rows |
|---|---:|---:|---:|
| Balanced atlas, all coordinate-eligible observations | 36 | 2,834–3,327 | 2 |
| Balanced atlas, positional accuracy ≤10 km | 36 | 2,029–2,388 | 0 |

These balanced-atlas screens must not be called the 46,276-observation primary cohort.

## Primary exhaustive within-species coefficients

The 46,276-observation spatially thinned cohort produced 36 component-wise models: nine endpoints × four CHELSA predictors. Eight component rows passed BH correction. Four were ordinary linear endpoints—orientation versus BIO1, chroma versus BIO12, and aspect ratio versus BIO4 and BIO12. Four additional rows were hue sine/cosine components. Because hue components are not independent biological endpoints, they are not counted as four additional colour conclusions; a joint circular test is required.

Accordingly, manuscript wording is:

> In the exhaustive spatially thinned primary cohort, four non-circular linear endpoint–climate associations passed BH correction and were small in standardized magnitude. Four hue-component rows also passed component-wise correction but were not interpreted separately.

## Revised precision-aware lability cohort

The original 102-taxon lability cohort allowed differing endpoint and predictor sets and used the RMS magnitude of unpooled slope estimates. That score was withdrawn after reviewer-style audit showed near-deterministic dependence on slope sample size and standard error.

The revised primary cohort contains 101 taxa and requires:

- all seven linear endpoints;
- all four predeclared CHELSA predictors for every endpoint;
- at least 10 observations for every species × endpoint × predictor slope;
- equal weighting of orientation, colour and shape modules;
- explicit use of every archived slope standard error.

Circular hue remains in the colour and PCA analyses but is excluded from this precision-corrected cross-species test because the archived species-level joint hue vectors do not contain component standard errors.

## Naming rules

- **Image-comparison atlas** = 6,626 heads, 3,725 observations, 216 taxa.
- **Exhaustive detected layer** = 406,582 observations, 286 taxa.
- **Exhaustive spatially thinned primary cohort** = 46,276 observations, 259 taxa.
- **Revised precision-aware lability cohort** = 101 fully complete taxa.
- Do not use *expanded pooled cohort* without a filename and count.
- Do not report an FDR count without the cohort, endpoint family and number of tests.

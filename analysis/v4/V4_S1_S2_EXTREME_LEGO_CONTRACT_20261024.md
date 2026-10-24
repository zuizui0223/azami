# azami V4 S1–S2: extreme bioclim + phenotypic LEGO (2026-10-24)

**Active target:** V4 PR #144 branch `codex/v4-native-integration-workflow`,
not V3 main. This is an **explicitly post-hoc, experimental V4 branch**.
Do not rewrite the frozen V3 result or call V3's 33/36 a V4 result.

## Scientific question

Are *Cirsium* floral colour, image-referenced capitulum presentation, and
two-dimensional capitulum outline fixed as an indivisible syndrome, or are
distinct visible trait combinations present both within and among taxa?
Does changing annual rain to warm-season/dry-quarter rain, adjusting for
temperature extremes and topography, alter those observable associations?
Is the native/introduced contrast within *C. arvense* and *C. vulgare*
different from the native between-taxon pattern?

## Reuse real V4 original source without changing its denominators

S1: 313,350 authority/WCVP-native original source observations;
170,066 orientation-assessable observations and 244 keys are the original
angle-layer statistics, **not** guaranteed complete for all 8 LEGO constructs.
Start from the original encrypted source, independently decrypt with the
existing Actions secret, verify exact SHA, and preserve head/QC eligibility.

S2: 297,622 WCVP-classified wild source observations in the two
original taxa *C. arvense* and *C. vulgare*, native+introduced combined.
Do not call other range states introduced. Same source-native images overlap
between S1/S2, so neither counts as independent validation of the other.

S1: use the actual 8 V4 constructs/28 pairs from
`analysis/v4/native_integration.py`:
orientation, lightness, chroma, joint circular hue,
aspect, circularity, solidity, width-profile CV.
These are **not** V3's 9 biological constructs/36 pairs.
Original S2 measurement columns are audited before deciding if
multi-trait native/introduced pairing is estimable; never fill missing
head shape from species means or photos without source evidence.

## New outcome-blind environmental candidates

Original V4 annual basis:
BIO01, BIO04, BIO12, BIO15 and RSDS.

Added *source-pixel original CHELSA V2.1 (1981–2010)*:
BIO05 warmest-month maximum; BIO06 coldest-month minimum;
BIO13 wettest-month rain; BIO14 driest-month rain;
BIO17 driest-quarter **mean monthly** rain;
BIO18 warmest-quarter **mean monthly** rain.
These are **NOT** independently measured rainfall extremes at flowering.

Add source-pixel **WorldClim V2.1 30arc-second elevation** (metres), a
different independent terrain source, for diagnostic/alternative adjustment.
Do not infer elevation effects from coordinates or from 5-degree cell means.
WorldClim elevation does not indicate the site-level timing of flowering.

S1 + S2 locations are combined only to avoid re-downloading identical
source pixels, with per-source obs_id checks and coordinate identity checks.
No user/observation/point identifiers are exported.

## Collinearity before fitting any new traits

On the original native environment rows, before any trait filter:
Pearson and Spearman matrices and VIF for 12 total candidates at three
levels (original observations, equal-taxon between means and taxon-centred
within deviations). Keep a fixed matched all-12 environmental support for
comparisons.

For each new interpretation, protect one focal precipitation scale, one
predeclared thermal scalar and RSDS; prune the same additional BIO04,
BIO15 and elevation nuisance terms by |r| > .8 and/or VIF > 5,
*without looking at any measured head trait*.
If the focal variables themselves are redundant, mark
UNRESOLVED and do not force a joint model.

Explicit candidate comparisons: annual BIO12 vs warmest-quarter BIO18,
driest-quarter BIO17, driest-month BIO14, wettest-month BIO13;
with/without elevation and alternate hot/cold BIO05/BIO06, using exactly
the same complete-12 × eight-trait S1 observations.

**V4 old model remains frozen.** New P values have a separate exploratory
multiplicity family; no retroactive original Holm adjustment.

## LEGO prediction checks against V4 phenotypes

**F (rigid package):** one dominating head configuration should imply
strong cross-module co-occurrence, few independent head presentations at
fixed colour/shape and little within-taxon reconfiguration.

**M (modular phenotype):** colour, presentation and outline may exhibit
different within-/between-taxon association geometry. Strong same-module
relations can coexist with weaker cross-module relations, with measurable
head configurations that cannot be represented by one fixed package.

**R (reconfigurable genetic LEGO):** needs one-plant functional G inventory,
independently validated hereditary module origins, compatibility, local
ancestry/cytotype, and actual ecological fitness. **No public image analysis
can prove R**. The V4 test can establish or falsify only visual phenotypic
predictions, and inform aza3's molecular design.

We compare all **28 V4 pairs** on an exact common sample at within and
among scales using RV, sign for scalar pairs, equal-species weighting,
environment-adjusted relations, and equal-n resampling. The 28 raw-RV
direction counts have finite-sample bias (as V3's 33/36 did). Do not
declare genetic independent modules, evolutionary transitions or adaptation
from them.

A three-way image configuration test for angle × chroma × aspect, if
admissible after QC, is a phenotype contrast only. It is *not* the anto
higher-order joint functional resource probability G123.

## Mandatory falsifiers and stop gates

- Identical three-candidate climate gradients unsupported or VIF >5:
  do not interpret independent annual vs summer rainfall effects.
- Joint 8-construct V4 cohort too small or missing an exposure:
  mark not identifiable and preserve trait/environment coverage.
- New S2 shape/colour original per-record columns unavailable:
  quantify missingness, do not create any simulated observed pair.
- Pairwise RV differences disappear under equal-n sampling, or are
  explainable by same-photo camera/colour/outline construction:
  no claim of strong scale-dependent module assembly.
- Geographic/camera/phylogeny adjustment not passed:
  no rain-caused adaptive architecture claim.
- Functional gene modules, joint genomic resources and viable achenes
  never measured here: no claim of **genetic** LEGO proof.

CHELSA units: https://www.chelsa-climate.org/datasets/chelsa_bioclim
Elevation source: https://geodata.ucdavis.edu/climate/worldclim/2_1/base/
aza3 functional definition:
https://github.com/zuizui0223/aza3/blob/main/docs/ANTO_HIGHER_ORDER_LEGO_MEASUREMENT_PLAN_V2.md

## Reproducible action

`analysis/v4/s1_s2_extreme_environment.py` downloads exact original
source rasters, records URLs, SHA256, spatial affine, scale, offset,
nodata and coverage; private extraction uses one original source pixel
per site, not 5-degree climate cell averages.

`analysis/v4/s1_s2_extreme_lego.py` retains all 12 candidate levels,
checks VIF, fits outcome-blind alternative environmental models, reports
V4 28-edge raw/conditional RV with equal-n diagnostics, and checks
S2 source-linked trait availability/range niches.

`.github/workflows/v4-s1-s2-extreme-lego.yml` guards private archived
source checksums. It publishes only aggregate JSON after genuine completion.

**Until an actual successful Action with seven new source layers, any new
climate/elevation coefficients or LEGO outcomes are UNCOMPUTED, not zero.**


## Explicit reduction to identify a conditional elevation sensitivity (24 October)

The original full-adjustment annual BIO12/elevation and BIO18/elevation
candidate models excluded elevation for VIF >5. Before interpreting that
as no elevation effect, diagnose a distinct and deliberately **reduced**
four-predictor sensitivity: original-pixel precipitation time-scale
(BIO12 **or** BIO18), BIO05 warm-season heat, RSDS, elevation. On the
all-12-original-climate source-only correlation matrices, both
four-variable designs pass VIF<5 at raw, equal-taxon and within-taxon
levels (raw maximum VIF approximately 4.98 and 4.96 respectively).
BIO4 and BIO15 are not adjusted in this sensitivity; the resulting
coefficient has a different, less conditional meaning. The samples,
bio05/rsds/elevation controls and 8 response constructs remain identical
between these two newly declared alternative rain time scales.
These are additional exploratory comparisons, not a revised original
S1 model, an elevation causality test, or proof that BIO12 beats BIO18.
S2 uses the same four adjusted source variables in both ranges of each
focal species, if its measured photo and geographic-cell support passes.

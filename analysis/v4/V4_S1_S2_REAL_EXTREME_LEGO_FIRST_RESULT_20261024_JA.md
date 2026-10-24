# V4 S1–S2 real extreme-climate and phenotypic LEGO results (24 October 2026)

**EXECUTED original-source exploratory results; not a V3 replay.** Original
V4 PR #144 branch is the parent; see [source-verified successful Actions
40057060469](https://github.com/zuizui0223/azami/actions/runs/40057060469)
and numerical artifact `12546132092` (original ZIP SHA256
`ccc9125f70c2f35b4fab89622969dda92e6035a7aa872b067e408b6fc68efc35`).
This note describes the **first fully successful** expanded-environment run,
before later V4 S2 climate regressions or candidate bract-module extensions
have been audited. Those later tests must not inherit this success statement.

## Original environmental source and covariate support

S1 original native = **313,350 observations**. S2 source-authenticated
wild native+introduced in two species = **297,622 observations**.
Identical original native records overlap by 141,014; exact ID+coordinate
deduplication leaves **469,958 genuine distinct observation locations**.

CHELSA V2.1 original 1981–2010 source grids:
BIO05, BIO06, BIO13, BIO14, BIO17 and BIO18 all sampled and finite at
**469,958/469,958 locations**; 30-arcsecond WorldClim V2.1 elevation
finite at **466,546/469,958 (99.27%)**. Source URLs, raster download SHA,
dtype, affine transform, scale/offset, nodata, row count and extraction
coverage are preserved in the original Actions source receipt. No
5-degree average, latitude-derived pseudo-altitude, or imputation is used.

S1 all-12-environment complete source **311,585 rows**, 303 source taxa.
Important actual Pearson correlations from the *source-only* census:

| Predictors | Raw-r | Equal-taxon mean r | Within-taxon r |
|---|---:|---:|---:|
| BIO14 versus BIO17 | +0.9909 | +0.9966 | +0.9844 |
| BIO12 versus BIO13 | +0.9229 | +0.9350 | +0.8848 |
| BIO1 versus BIO6 | +0.9025 | +0.9245 | +0.8782 |
| BIO12 versus BIO18 | +0.7981 | +0.8956 | +0.6720 |

All 12 predictors jointly have raw VIF up to **104.98 (BIO6)** and
equal-taxon mean VIF up to **218.56 (BIO17)**. They MUST NOT be simultaneously
entered as independent regressors. Deterministic *outcome-blind*
alternative-model selection retains annual BIO12 OR warmest-quarter BIO18.
The originally attempted elevation nuisance is **removed by the VIF rule**
from both primary annual and warm-quarter models; this is NOT a biological
null test for elevation. With BIO6, BIO18 and the alternative thermal
combination, elevation can remain while BIO4 is removed; these are
**different estimands**.

## V4 image phenotypes: 156,683 observations / 174 taxa

The original 8 image constructs and their QC yield
**157,654 observations / 174 taxa** before the new all-12 environmental
complete requirement and **156,683 / 174** after same-observation
trait+climate matching (minimum five complete observations per taxon).

V4 uses **8 constructs = 28 pairs**, NOT V3's
9 constructs / 36 pairs. Exact same-photo/construct sample throughout
the following annual vs seasonal comparisons.

| V4 same-cohort statistic | Unadjusted source phenotypes | Annual BIO12-adjusted | BIO18-adjusted |
|---|---:|---:|---:|
| Pairs with among-RV > within-RV, raw estimator | **21/28** | **23/28** | **23/28** |
| Median within-taxon RV | 0.003481 | 0.003459 | 0.003386 |
| Median among-taxon RV | 0.022713 | 0.019054 | 0.021316 |
| Median equal-n among-minus-within RV, n=taxa matched | **+0.014037** | **+0.010230** | **+0.012460** |
| Equal-n 500-resampling 2.5–97.5% range | **+0.003826 to +0.019002** | **+0.000296 to +0.015383** | **+0.002511 to +0.017690** |
| Median equal-n count with stronger among | **19/28** | **19/28** | **19/28** |

Do not report the raw 21/28 or 23/28 as a sample-size-independent
statistical success rate. The equal-n ranges conditional on original
taxa and fixed among estimates are **not population confidence intervals**.
This is an observational shift in image phenotypes, not evolutionary
acceleration or independently derived modules.

Among the actual raw edges, image lightness × hue has RV within **0.1042**,
among **0.4273**; chroma × hue has **0.0601 / 0.3482**;
orientation × 2D elongation has **0.0048 / 0.0722**.
The colour pair coherences involve shared pixel sources, and the orientation
× elongation combination may reflect camera perspective.

The original a priori module grouping shows a descriptive pattern:
- Within taxa: median **same-domain RV 0.08498**, cross-domain **0.00140**.
- Among taxa: median **same-domain RV 0.22384**, cross-domain **0.01062**.

Independent algebraic enumeration of **280 equal-size label partitions**
(1 presentation / 3 colour / 4 outline) gives upper-tail label-reference
P **1/280 = 0.00357** for within, **2/280 = 0.00714** among.
This is a **posthoc technical module-label reference**, not a phylogenetic
permutation or evidence of independent genetic LEGO. It needs its separate
CI result to be elevated beyond independent aggregate re-calculation.

## A three-piece LEGO test beyond pairwise correlations

The prespecified visible triad is image-relative angle × corolla chroma ×
2D outline elongation. With *within-taxon-median* dichotomization, **124/135**
well-observed taxa have all 8 visible low/high combinations. This is a
positive **phenotypic-configuration observability** result; high photo
sample sizes alone can make all eight appear, and no independent heritable
genotypes or dispersal event histories were measured.

The actual third-order test must distinguish it from pairwise information.
Under *taxon-mean*-centred binary configurations, fixed single and pairwise
margins imply a maximum-entropy reference of **0.08995** for the within-taxon
high/high/high frequency; observed is **0.08696**, difference **−0.00299**
(taxon-bootstrap resampling range **−0.00663 to +0.00078**).
Among taxa, observed **0.08046** and expected **0.08046**, difference
effectively **0** (range about **−0.019 to +0.017**).
The separate within-taxon-*median* parity statistic is +0.03849 with a
positive taxon-bootstrap interval, i.e. **the inferred higher-order excess
depends on dichotomization and test definition**. Do not claim stable
third-order organisation or functional G123 co-retention from this discrepancy.

## S2: real warmest-quarter climate expansion in two introduced taxa

Using all-12-climate complete original WCVP wild source locations in each
range (not the previous BIO12-only denominator):

| Species | Introduced minus native BIO18 median (mean monthly precipitation in warmest quarter) | Introduced minus native BIO12 median |
|---|---:|---:|
| C. arvense | **+48.2 mm/month** | **+263 mm/year** |
| C. vulgare | **+49.7 mm/month** | **+298 mm/year** |

These two cases are positive **occupied climate distribution** differences,
not within-individual plant evolution; the earlier original BIO12-only
comparison gave +217.4 and +293 mm with different sample coverage.
The actual S2 original image measurements for angle, colour and outline
are present in both species and both ranges (tens of thousands of
measurements per primary coordinate). New S2 same-cohort fitted slopes and
uncertainty should be quoted **only after their own source-verified run**,
not manufactured from the above between-range medians.

## Scope of aza3 LEGO bridge

azami V4 only tests *visible phenotypic packages* in authentic original
source heads; aza3's stronger genetic/developmental/functional modularity
requires genomic resource identities and same-individual developmental,
compatibility, phylogeny/cytotype and viable-achene tests. The protective
phyllary envelope is only an image projection proxy in the current V4
low-support SI cohort, NOT measured true spine, gland, adhesion or defense.
A physically unbreakable fixed whole-capitulum package is not established;
observed reconfigurations motivate the aza3 experiment, but genetic LEGO
is **not yet proven**.

**Reproduce source:** original aggregate artifact `12546132092` and
`analysis/v4/lego_module_label_null.py` for the 280-label post-hoc diagnostic.

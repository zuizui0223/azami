# Chapter 1 minimal v3 revision

This directory is intentionally small. The revision does **not** replace the frozen v2 ecological analysis with a new full-source model.

## Canonical path

1. **Source provenance** — reconcile the original metadata/API archives so the v2 processing universe is traceable upstream.
2. **Historical measurement reuse** — recount the existing all-photo processing archive and recover the five head fields that were measured but omitted from the historical observation aggregation.
3. **v2 ecology** — reuse the existing v2 ecological scripts and the frozen 46,276-observation / 259-taxon primary cohort.
4. **Output** — rebuild figures/tables/manuscript-facing summaries with corrected claim boundaries.

Run:

```bash
python analysis/v3/minimal_preflight.py
```

## What v3 fixes

- Restores the upstream source ledger rather than treating the v2 analysis table as the acquisition universe.
- Recounts the historical queue/screen/head pipeline without downloading or remeasuring images.
- Recovers five already-measured display/colour-composition fields that the historical aggregator omitted.
- Keeps the old responsiveness-versus-variation headline out of the canonical claim set.
- Treats below-species-mean image variation as visible image variation, not genetic variance.
- Treats the v2 environment coefficients as marginal associations along correlated observational gradients.

## What v3 no longer requires

The following are preserved only on `archive/pr92-preclean-20260909` and do not block the canonical revision:

- Wave B/C image remeasurement;
- independent assessability research by taxon/region/environment;
- alternative uniform-chroma qualification;
- tail/repair/full-family/final calibration branches;
- new 36-slot calibration expansions;
- hypervolume/breadth synthesis;
- all-27 exploratory results as a route-admission rule;
- a new 319,244-observation primary ecological cohort.

The principle is: **reuse v2 wherever it is still valid, patch only the parts that v2 could not support cleanly.**

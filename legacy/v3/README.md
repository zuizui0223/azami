# Frozen Chapter 1 V3 — archived 25 October 2026

**Status: legacy only.** The current Research Article is V4 S1–S2, linked from [root README](../../README.md). The V3 study is **not superseded numerically** or claimed invalid; it has a different construct family, sampling cohort and inferential scope. Its former `main` was commit `b463dbc938e454baa54b2e480f7f1d8badd7795d` (28 September 2026).

## Location and original path mapping

| Historical path | Archived path | Treatment |
|---|---|---|
| `analysis/v3/**` | `legacy/v3/analysis/v3/**` | moved exactly; original tree blob identities preserved |
| `reproducibility/**` | `legacy/v3/reproducibility/**` | moved exactly; includes historical figure and v2 reference outputs |
| `.github/workflows/**` (14 workflows at archive point) | `legacy/v3/workflows/active_as_of_20260928/**` | moved exactly; no longer executable Actions on `main` |
| `legacy/workflows/**` (4 previously archived V3 workflows) | `legacy/v3/workflows/previously_archived/**` | moved exactly |
| Dedicated V3 `tests/test_*.py` | `legacy/v3/tests/test_*.py` | historical checks; their source revision must be checked out for replay |

Shared `analysis/ch1/`, `src/azami_ch1/` and V2 `analysis_outputs/` remain in their earlier locations because V2/V3 rely on some of them. No v3 statistical numbers are migrated or recalculated.

## Reproduce historical V3

```bash
git checkout b463dbc938e454baa54b2e480f7f1d8badd7795d
# Reproduce using that version's reproducibility/CURRENT_ANALYSIS.md
```

From the current checkout, **do not** run old `python -m analysis.v3...` commands or historical Actions definitions: the archive retains their original relative imports and paths. This is intentionally a historically replayable snapshot, not an active pipeline. Original legacy README may describe V3 as "current" because it is preserved as an exact historical receipt.

V3 snapshot described nine constructs / 36 pair relations on a strict 1,734-observation / 42-taxon common cohort and a distinct 46,276-observation / 259-taxon atlas. **These are not the V4 eight constructs / 28 pairs on 156,683 photographs / 174 taxa.**

The historical Zenodo v3 publishing workflow must not be manually copied back to active `.github/workflows` without a separate explicit authorization and new publication-safety review.

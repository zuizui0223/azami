# Chapter 1 v3 — continuous capitulum traits along environmental gradients

## Active environment-first extension of PR93

The main question is how the constituent biological traits of a capitulum covary along environmental gradients, together or separately. Reuse PR93 at `4deae0815880baf5db3e2795bd79369344e8bace`; retain its nine non-surface constructs as the biological core. Surface measurements remain supporting diagnostics. Within- and among-taxon results remain separate estimands, but are presented under the same environmental gradients rather than as competing paper narratives.

The workflow is continuous measurement -> biological trait definition -> complete environment atlas with inherited FDR -> sampling/spatial/among-taxon placement sensitivity -> whole-capitulum synthesis. No new variance-partitioning model, native-only filter, image acquisition, or categorical-superiority test is added. Hypervolume plots can illustrate distributions but do not replace environmental tests.

The first executed step reorganizes the existing PR93 numerical results into 162 core rows (9 constructs x 9 predictors x 2 scales). Of these, 21 retain the original FDR support; 3 among-taxon rows pass the full declared sensitivity sequence: presentation angle–BIO12, floral chroma–radiation, and floral chroma–NPP. NPP remains exploratory, not a new frozen-v2 headline. No within-taxon core row passes the entire sampling/spatial sequence. These are inherited results, not independent replication. Correlated environmental predictors remain marginal associations.

Whole-capitulum interpretation combines the gradient response map with PR93's separately established common-cohort integration; counting supported traits alone does not establish coordinated environmental response or a syndrome. A direct joint response test is not yet added.

Numerical inputs were downloaded from Actions run `34418904597`, artifact `10130210432` (`ch1-v3-biological-axes-34418904597`). Public aggregate source tables and their hashes are saved with the generated outputs at `analysis_outputs/environment_first_20260910/`; no raw coordinates or reviewer replies are included. Scientific claims remain subject to the measurement and placement limitations documented below.

```bash
python analysis/v3/build_environment_first_synthesis.py --source analysis_outputs/environment_first_20260910/source --out local_data/environment_first_replay
python -m pytest tests/test_environment_first_synthesis.py -q
```

The following sections preserve the PR93 baseline and its provenance.

## Canonical boundary

The frozen Chapter 1 **v2 remains the canonical endpoint-level analysis**. This directory does not replace its endpoint atlas, multiplicity families, public reproduction bundle, or two manuscript-level headline ecological conclusions.

The v3 work asks a higher-level question using the same frozen measured endpoints:

> How are visible capitulum traits biologically organized within taxa and among taxa, how does that organization change with scale, and do the two frozen-v2 ecological anchors remain robust inside that broader reorganization?

## What v3 adds

1. **Biological reaggregation** — the 22 measured v2 endpoints are reorganized into named biological constructs without using environmental outcomes to choose the grouping.
2. **Scale-dependent integration** — construct–construct integration is compared directly within taxa and among taxa.
3. **Common-cohort robustness** — all nine complete-18 constructs and all 36 relations are forced onto the same 1,734-observation / 42-taxon cohort.
4. **Uncertainty and modularity** — taxon bootstrap and module-label permutations quantify partial cross-scale conservation and biological module cohesion.
5. **Direct scale contrast** — among-taxon visible-phenotype integration is stronger overall than within-taxon integration across the common-cohort taxon bootstrap.
6. **Ecological-anchor preservation** — floral chroma × shortwave radiation and presentation angle × annual precipitation survive biological reaggregation plus sampling, spatial/residual and 52-tree historical sensitivity gates.
7. **Image-based alternative audits** — outcome-blind assessability and an explicit frozen-summary-calibrated orientation stress test narrow selected measurement/selection explanations without claiming physical-trait accuracy.

## Current core results

- common-cohort within–among integration alignment: **rho = 0.439125**, QAP **P = 0.0041**;
- taxon-bootstrap alignment: median **0.31210**, 95% interval **0.01708–0.53517**, **97.9%** positive;
- module cohesion: within-taxon **P = 0.0013**, among-taxon **P = 0.0365**;
- observed median RV: **0.002238 within** vs **0.043212 among**;
- **33/36** relations are stronger among taxa;
- 1,000-bootstrap difference of median RV (among − within): median **+0.065928**, 95% interval **+0.038197 to +0.119048**, positive in **100%** of replicates;
- construct-level frozen anchors remain floral chroma × radiation (**beta = -0.345372**) and presentation angle × annual precipitation (**beta = +0.304359**), and both pass the full construct-level robustness chain.

## Result documents

- `BIOLOGICAL_AXIS_REANALYSIS_PLAN_20260910.md` — construct definitions and analysis plan.
- `BIOLOGICAL_AXIS_RESULT_SUMMARY_20260910.md` — 22-endpoint → biological-construct reanalysis.
- `BIOLOGICAL_AXIS_SENSITIVITY_CHAIN_RESULT_20260910.md` — sampling → spatial/residual → historical sensitivity at construct level.
- `SCALE_DEPENDENT_INTEGRATION_RESULT_20260910.md` — pairwise-cohort within/among integration synthesis.
- `CONSTRUCT_SCALE_UPGRADE_RESULT_20260910.md` — exact complete-18 common cohort, bootstrap, module cohesion and six environmental blocks.
- `DIRECT_SCALE_CONTRAST_RESULT_20260910.md` — direct uncertainty-qualified among-vs-within integration-strength contrast.
- `ASSESSABILITY_AND_TECHNICAL_STRESS_RESULT_20260910.md` — outcome-blind availability audit and frozen-summary-calibrated orientation stress.

## Interpretation boundary

Supported language:

> The capitulum is neither one undifferentiated syndrome nor a set of unrelated image measurements. Biological modules are detectable at both scales, but their detailed integration geometry is only partly conserved, and visible-phenotype integration is stronger overall among taxa. Against this broader scale-dependent reorganization, the two frozen-v2 among-taxon ecological anchors remain robust.

Do **not** infer from these analyses alone:

- genetic or developmental modularity;
- plasticity from within-taxon spatial covariance;
- causal environmental effects;
- gravity-referenced accuracy of image vertical;
- calibrated physical colour accuracy;
- unbiased image sampling or absence of value-dependent missingness.

Colour error propagation remains unresolved because the frozen technical audit preserves no numeric per-image colour-error distribution. Orientation random-stress results do not address systematic camera roll or environment-dependent measurement error.

## Reproducibility receipts

- construct reanalysis / sensitivity chain: see the receipts linked in the corresponding result documents;
- common-cohort upgrade: run `34432661967`, artifact `10135139012`;
- assessability / technical stress: run `34434467062`, artifact `10135679053`;
- direct scale-contrast rerun: run `34435743710`, artifact `10136229131`, digest `sha256:7b4c5703b048ff108cd794921394ae22d1423461c9aaf330b47029c5d1f4908f`.

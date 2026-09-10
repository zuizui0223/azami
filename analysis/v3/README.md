# Chapter 1 v3 — shortest reuse route

## Active PR92 exploration since 2026-09-10

PR92 is now the full-range capitulum distribution / hypervolume exploratory branch, not a native-only primary reanalysis. Reuse the 46,276-observation v2 measurement universe, without new image acquisition or a native-range filter. PR93's construct integration analyses are a separate workstream.

The first bounded pilot is `run_distribution_breadth_pilot.py`: equal-sample density breadth versus environmental/geographic breadth, followed by spatial-cell-held-out within-taxon prediction. Public aggregate outputs are in `analysis_outputs/pr92_distribution_breadth_pilot_20260910/`. This pilot uses 50 observations per taxon in 30 rarefactions; its KDE volumes describe common projections, not full biological trait space. At most three axes are retained for this pilot, with retained variance reported explicitly. Involucre and whole-capitulum projections have too few eligible taxa and too much discarded variance for whole-trait-volume inference.

Breadth correlations have not been corrected for spatial/phylogenetic dependence or multiplicity. The geographically held-out prediction comparison is a descriptive screening result, not a significance test; coordinate calibration uses the full trait sample and is not a nested independent measurement validation. Do not turn positive predictive gain into proof of plasticity. Environmental breadth refers to sampled precipitation, radiation, VPD and wind exposures, not a complete ecological niche. Historical runs below are preserved, not the active specification.

Run with existing numerical inputs:

```bash
python analysis/v3/run_distribution_breadth_pilot.py --traits continuous_trait_universe_observation_long.csv --environment strict_spatial_chelsa_full9.csv --out local_data/distribution_breadth_pilot
python -m pytest tests/test_distribution_breadth_pilot.py -q
```

The follow-on `run_distribution_overlap_pilot.py` reuses exactly these shared bases, input hashes and 50-observation eligibility. It estimates density overlap (integral of the smaller of two KDE densities) and exact equal-weight empirical Wasserstein W2. W2 squared is separated into centroid distance squared and centred-distribution distance squared. Thirty rarefaction repetitions are summarized; their percentile ranges are not confidence intervals. Pairwise values are descriptive, not independent hypothesis tests. Overlap depends on fixed smoothing, and finite samples create positive W2 even for identical underlying populations. Neither overlap nor distance demonstrates functional equivalence or convergence.

The completed local pilot and numerical checks are recorded in `analysis_outputs/pr92_distribution_overlap_pilot_20260910/`: 4,318 module-specific taxon pairs and 129,540 pair-repetitions. Full replicate arrays remain in the local run directory and can be regenerated from the script; public tables contain aggregate pair results. This is not a new independent image-validation dataset.

```bash
python analysis/v3/run_distribution_overlap_pilot.py --traits continuous_trait_universe_observation_long.csv --environment strict_spatial_chelsa_full9.csv --basis-dir analysis_outputs/pr92_distribution_breadth_pilot_20260910 --out local_data/distribution_overlap_pilot
python -m pytest tests/test_distribution_overlap_pilot.py tests/test_distribution_breadth_pilot.py -q
```

## Full-range figure and finite-sample follow-up

`build_distribution_pilot_figures.py` produces the breadth scatter and rule-selected distribution examples, plus bounded sampling and smoothing checks. The saved products are in `analysis_outputs/pr92_distribution_figures_20260910/`.

Same-taxon comparisons use 30 disjoint 50+50 draws, so only taxa with at least 100 eligible observations enter. Each existing between-taxon median W2 is compared with the larger of its two taxon-specific same-taxon 95th percentiles. This is a conservative descriptive reference, not a calibrated hypothesis test or biological variance partition. The centred-distance reference is computed separately.

The earlier 84% median centred share for outline is not, by itself, evidence of strong taxon differentiation: only 47/561 outline pairs exceed the centred same-taxon reference, compared with 13/528 orientation and 349/528 colour pairs. These fractions are not estimates of statistically significant pairs and pairs are not independent. Photography and measurement variation remain inseparable from phenotype variation.

Bandwidth checks use an outcome-blind seeded subset of 100 pairs per core module, ten matched rarefactions, and multipliers 0.75/1/1.25. Rank correlations exceed 0.97 against the baseline multiplier, but absolute overlaps change. Do not present absolute overlap percentages as estimator-independent results. Figure examples are explicitly outcome-selected illustrations, not independent validation.

```bash
python analysis/v3/build_distribution_pilot_figures.py --traits continuous_trait_universe_observation_long.csv --environment strict_spatial_chelsa_full9.csv --breadth-dir analysis_outputs/pr92_distribution_breadth_pilot_20260910 --overlap-dir analysis_outputs/pr92_distribution_overlap_pilot_20260910 --out local_data/distribution_figures
python -m pytest tests/test_distribution_pilot_figures.py tests/test_distribution_overlap_pilot.py tests/test_distribution_breadth_pilot.py -q
```

## Historical native-only route

The earlier v3 route was deliberately small. It repaired upstream provenance and aggregation, then reused the frozen v2 ecological machinery rather than rebuilding a second analysis system.

## Current path

`frozen source provenance -> frozen strict-spatial cohort -> native-range restriction -> historical measurement reuse -> phenotype-blind VIFstep<=10 predictor filter -> frozen v2 within/among estimators -> existing broad-spatial sensitivity -> existing 52-tree historical-placement sensitivity`

### Cohort

- Frozen v2 strict-spatial universe: **46,276 observations / 259 source-assigned taxa**.
- Current native-only reanalysis: **27,066 observations / 241 taxa with environment data**.
- Native membership is deterministically regenerated from frozen v2 auxiliary tables because the historical native-status Git LFS object is no longer byte-recoverable.

### Traits

- All **27 registered endpoints** are retained.
- Hue sine/cosine remain one joint circular inferential unit, giving **26 inferential units**.
- Five fields already measured historically at head level but omitted from the old observation aggregation are restored; no new image measurement is performed.
- Closed colour fractions remain descriptive and are not narrated as independent primary biological discoveries.

### Within and among scales

This distinction is inherited from v2, not introduced by v3.

- **Within taxon:** taxon-demeaned standardized marginal slope, taxon-clustered uncertainty; circular traits use the registered joint test and within-taxon predictor permutations.
- **Among taxon:** taxon-level trait medians and environment medians, standardized marginal slope, predictor permutation across taxa; min5 primary and min2 sensitivity scopes are retained.
- BH-FDR families remain separate for within, among-min5 and among-min2.

### Environment

The frozen v2 nine-predictor atlas is retained as the historical benchmark. A user-requested comparative lane applies phenotype-blind VIFstep<=10 before rerunning the same v2 marginal models.

For the native cohort, VIFstep removes **BIO1 only** at both scales. The retained eight predictors are BIO4, BIO12, BIO15, shortwave radiation (`rsds`), VPD, surface wind, GSP and NPP.

Final maximum VIF is **1.961** for within-taxon demeaned environment and **6.977** for among-taxon median environment. VIF is used only to filter the family of marginal tests; it does not convert the analysis to a simultaneous multivariable regression.

## Current shortest-route result

VIFstep<=10 native reanalysis:

- within FDR: **37**
- among-min5 FDR: **27**
- among-min2 FDR: **16**
- within broad-spatial passes: **6**
- among broad-spatial passes: **4**
- among pairs passing all 52 historical-placement trees: **4/4**

The two primary among-taxon survivors are `corolla_lab_chroma x chelsa_rsds_mean` (negative) and `corolla_lab_chroma x chelsa_npp` (positive). Two closed-composition descriptive rows also survive: `corolla_yellow_pixel_fraction x chelsa_gsp` and `corolla_yellow_pixel_fraction x chelsa_npp`.

`orientation_image_vertical_angle x chelsa_bio12` remains in the retained predictor family but does **not** pass the existing broad-spatial sensitivity because residual Moran structure remains detectable. VIF filtering does not rescue it. `orientation x BIO1` is not tested in this lane because BIO1 is the phenotype-blind VIFstep exclusion.

## Receipts

- `native_primary_spatial_historical_receipt_20260909.json` — native-only nine-predictor benchmark.
- `native_primary_vifstep10_receipt_20260909.json` — current VIFstep<=10 comparative result and exact workflow/artifact provenance.

## What v3 no longer requires

The following are preserved only on `archive/pr92-preclean-20260909` and do not block the current route: Wave B/C image remeasurement, independent assessability research, alternative uniform-chroma qualification, tail/repair/calibration branches, hypervolume/breadth synthesis, or a new 319,244-observation primary ecological cohort.

## Boundary

This remains a pattern-first observational analysis. Broad-spatial sensitivity is not a full SPDE process model, the 52 trees are historical-placement sensitivity rather than a resolved *Cirsium* species tree, and no surviving association establishes causation, adaptation, selection or mechanism.

The principle is: **reuse v2 wherever it is still valid, patch only the parts that v2 could not support cleanly.**

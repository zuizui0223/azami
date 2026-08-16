# azami — global image-derived capitulum traits in *Cirsium*

This repository contains a multi-chapter project on the ecology and evolution of *Cirsium* floral architecture. Chapter 1 is the submission-focused component: a global analysis of continuous capitulum traits measured from public biodiversity photographs.

## Chapter 1 in one sentence

Species means conceal hierarchically structured visible phenotype distributions, but no common coupling between visible dispersion and within-species spatial climatic association remains after slope uncertainty is modelled.

The analysis is observational and non-causal. Spatial climatic association is not proof of temporal response, local adaptation, phenotypic plasticity, pollinator selection or evolutionary rate.

## Role in the wider research program

Chapter 1 is the **global macro-scale hypothesis-generation layer**. It deliberately views *Cirsium* from far away using a very large public-image dataset before a resolved species history is imposed.

The next evolutionary-resolution stage is maintained in `zuizui0223/EAzami`:

```text
azami / Chapter 1
Global public-image macro screen
        ↓ hypotheses and trait distributions
EAzami
East Asian nuclear phylogeny + explicit transition histories
        ↓ replicated focal transitions
population / mechanism studies
Ancestry + expression + pigment + biological interaction + fitness
```

Therefore Chapter 1 does **not** need to absorb definitive ancestral-state reconstruction, transition counts, pollinator causation or molecular regain mechanisms before submission. Those are downstream tests. The explicit boundary and handoff are recorded in `manuscript/EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md`.

Chapter 1 already contains an **exploratory auxiliary involucre layer** derived from high-resolution head crops. It measures image-geometry proxies including involucre projection roughness, involucre spread fraction, spine-like peak count and maximum relative spine-like projection length, and screens their within-species and between-species climatic associations. These are auxiliary image proxies, not botanically validated phyllary-angle or spine measurements. Downstream work should carry these existing results into EAzami, then test whether they track explicit phyllary spreading/recurvature and spine architecture on the resolved nuclear history. Truly new/refined traits include botanically validated phyllary angle, actual spine length/orientation, visible involucral stickiness/glandularity and display architecture.

## Canonical executed cohorts

Two data streams answer different questions.

| Cohort | Observations | Taxa | Purpose |
|---|---:|---:|---|
| Balanced image-comparison atlas | 3,725 | 216 | 6,626-head nested visible variance, species PCA, among-species summaries and historical sensitivity |
| Exhaustive detector-positive layer | 406,582 | 286 | Post-detection source layer |
| Exhaustive spatially thinned primary | 46,276 | 259 | Primary within-species climate coefficients |
| Revised precision-aware lability cohort | 101 | 101 | Seven linear endpoints × four predictors with all slope standard errors |

The full flow from 777,766 photographs to each derived table is recorded in `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md`.

## Nested image hierarchy correction

The balanced atlas retains one photograph per public observation and may contain multiple detected heads per photograph. The earlier two-level variance partition combined these levels. A revised exact decomposition separates:

- differences among assigned species means;
- differences among photographs within assigned species; and
- differences among heads within photographs.

Across nine endpoints, the combined below-species fraction is 0.589–0.931. The among-photograph component is 0.440–0.691 and the among-head-within-photo component is 0.143–0.379. One-head-per-photo sensitivities retain fractions of 0.582–0.899. Equal-replication sensitivities using 10 photographs per eligible species retain median fractions of 0.528–0.879.

These are visible image sums-of-squares components, not genetic or standardized phenotypic variance. The correction is implemented in `analysis/decompose_nested_visible_variance.py` and `.github/workflows/ch1-nested-visible-variance.yml`.

## Reviewer-driven statistical correction

The previous lability analysis used the RMS magnitude of unpooled absolute species-specific slopes. That score was strongly confounded with median slope sample size (Spearman rho = -0.919) and median slope standard error (rho = 0.914). The previously reported rho = -0.333 and median-split quadrants are withdrawn.

The replacement analysis:

- uses a common 101-taxon cohort with all seven linear endpoints and all four climate predictors;
- weights orientation, colour and shape equally;
- subtracts slope sampling variance through `beta^2 - SE^2` for a transparent association-energy summary;
- fits a hierarchical variance meta-regression to all 2,828 slope estimates with their standard errors.

No common coupling was detected:

- noise-adjusted axis correlation: rho = -0.057, P = 0.572, 95% species-bootstrap CI -0.265 to 0.155;
- hierarchical log-variance coefficient: -0.042, 95% profile CI -0.274 to 0.200, likelihood-ratio P = 0.732.

The correction is implemented in `analysis/reanalyze_lability_precision.py` and rerun by `.github/workflows/ch1-reviewer-precision-reanalysis.yml`.

## Other retained results

Species-level PCA remains multidimensional: PC1 explains 32.9%, PCs 1–2 explain 56.1% and PCs 1–3 explain 69.3%.

In the exhaustive 46,276-observation primary cohort, eight of 36 component rows pass BH correction. Four are non-circular linear associations and all are small in standardized magnitude. Four additional hue sine/cosine rows require joint circular interpretation. The earlier balanced-atlas ≤10 km sensitivity has zero BH-supported main rows; it is a different dataset and is not the primary exhaustive FDR result.

Among-species environmental sorting remains strongest for orientation and visible colour and weaker for most outline traits. Historical conclusions remain sensitive to grafted taxon placement because only 54 of 216 image-analysis taxa occur directly in the dated backbone.

## Canonical manuscript path

Start here:

1. `manuscript/SUBMISSION_MANUSCRIPT.md` — submission-facing section order and status;
2. `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — immutable cohort names, counts and analysis permissions;
3. `manuscript/final_claims.json` — reviewer-revised machine-readable claim registry;
4. `manuscript/results/nested_visible_variance_summary.csv` — image-hierarchy decomposition;
5. `manuscript/results/reviewer_precision_summary.json` — numerical precision audit;
6. `analysis/ch1/pipeline.json` — canonical executable stages;
7. `manuscript/EXTERNAL_COMPLETION_GATES.md` — detector, measurement, taxonomy, spatial and archive gates still requiring external input;
8. `manuscript/EAzami_HANDOFF_AND_CHAPTER_BOUNDARY_2026-08-16.md` — explicit macro→phylogeny→mechanism handoff and Chapter 1 stop rules.

The validator `analysis/validate_final_claims.py` and workflow `.github/workflows/ch1-final-manuscript-claims.yml` enforce the revised claim set.

## Primary traits

- orientation angle relative to EXIF-oriented image vertical;
- corolla Lab lightness and chroma;
- circular hue sine/cosine;
- capitulum aspect ratio, circularity, solidity and width-profile variation.

### Exploratory auxiliary involucre traits already analysed

- involucre projection roughness;
- involucre spread fraction;
- spine-like peak count proxy;
- maximum relative spine-like projection length proxy.

These auxiliary endpoints are measured only on QC-passing high-resolution crops and are explicitly treated as image-geometry proxies. They have been included in continuous environment screening, including within-species and between-species analyses, with auxiliary historical sensitivity retained as exploratory rather than promoted to the Chapter 1 headline.

Circular hue remains a joint endpoint. It is retained in colour and PCA analyses but excluded from the precision-corrected cross-species lability test because archived species-level joint hue vectors lack component standard errors.

## Reproducibility policy

- Every FDR count must name its cohort, endpoint family and number of tests.
- Multiple heads from one photograph must remain a nested level or be reduced by one-head-per-photo sensitivity.
- Raw absolute species-slope RMS and median-split quadrants are legacy provenance only.
- The revised primary lability test requires all seven linear endpoints, four predictors per endpoint and n ≥ 10 for every slope.
- Final artifacts must include CSVs, figures, metadata, package versions, commit SHA, workflow run IDs and SHA-256 checksums.
- GitHub Actions artifacts are temporary; the submission bundle must be archived durably before submission.

## Interpretation limits

Citizen-science photographs are not colour calibrated, image vertical is only a proxy for gravity and outline measures remain viewpoint dependent. The auxiliary involucre/spine variables are image-geometry proxies rather than direct botanical measurements. Detector and continuous-measurement accuracy require independent validation. Taxonomic decisions and residual spatial/region-omission diagnostics remain submission gates. The grafted mega-tree is a historical sensitivity analysis, not a resolved evolutionary history.

# azami — global continuous capitulum phenomics in thistles

This repository is the submission workspace for a global analysis that converts public biodiversity photographs into repeated **continuous head-level trait measurements** rather than one categorical state or one mean value per taxon.

## Submission in one sentence

Global public photographs recover continuous within-taxon distributions of capitulum orientation, visible colour and outline that are largely hidden by taxon means, with trait- and scale-specific environmental structure strongest for orientation and visible colour.

The analysis is observational and non-causal. Spatial climatic associations are not proof of temporal response, local adaptation, phenotypic plasticity, pollinator selection, evolutionary rate or adaptive radiation.

## Canonical path

Start here:

1. `manuscript/SUBMISSION_MANUSCRIPT.md` — manuscript-facing story and section order;
2. `manuscript/COHORT_FLOW_AND_ANALYSIS_LEDGER.md` — immutable cohort names, counts and analysis permissions;
3. `manuscript/final_claims.json` — machine-readable frozen claims;
4. `analysis/ch1/pipeline.json` — canonical executable stages and rerun policy;
5. `manuscript/supplement/` — Supplementary Information and machine-readable tables;
6. `manuscript/EXTERNAL_COMPLETION_GATES.md` — genuinely external remaining validation.

Executed numbered scripts under `ch1_global/v2/` are retained at their original paths for provenance. They are implementation targets, not the navigation interface.

## Methodological contribution

```text
public biodiversity photograph
        ↓
YOLO localization of visible capitula
        ↓
deterministic continuous measurement per assessable head
        ↓
orientation + visible colour + outline geometry
        ↓
repeated within-taxon trait distributions
        ↓
nested variance + within-taxon environment association + among-taxon sorting
```

Inferential traits are continuous numerical measurements. Failed or unassessable measurements remain missing rather than being coded as biological absence.

## Canonical executed cohorts

| Cohort | Observations | Taxa | Purpose |
|---|---:|---:|---|
| Balanced image-comparison atlas | 3,725 observations / 6,626 heads | 216 | Nested visible variance, PCA, among-taxon summaries |
| Exhaustive detector-positive layer | 406,582 | 286 | Post-detection source layer |
| Exhaustive spatially thinned primary | 46,276 | 259 | Primary within-taxon climate coefficients |
| Grouped SPDE-INLA complete cases | 31,666–34,472 per endpoint | 139–141 | Spatially explicit within-taxon models |
| High-resolution involucre subset | 1,292 observations / 1,443 usable heads | 210 | Auxiliary contour proxies |
| Independent detector audit | 1,000 source images | 323 | Human validation; annotation pending |

## Main retained results

- Across nine primary endpoints, 0.589–0.931 of visible image variance occurs below source-assigned taxon means.
- The 46,276-observation primary climate family has 8/36 BH-supported component rows: four non-circular linear rows and four circular-hue components requiring joint interpretation.
- Grouped SPDE-INLA most consistently retains environmental structure for orientation and visible colour; gross-outline effects are weaker/model-dependent.
- In the high-resolution auxiliary layer, involucre roughness, outward spread and maximum spine-like projection increase with temperature seasonality; these remain 2-D image proxies rather than demonstrated defensive traits.
- Among-taxon environmental sorting exceeds permutation-null expectation mainly for orientation and visible colour.
- Residual Moran's I is near zero for all nine primary endpoints, broad-region omission retains taxon rankings strongly, and simultaneous collapse of eight WCVP synonym conflicts changes neither primary coefficient signs nor FDR decisions.
- Historical/PGLS results remain sensitivity analyses because only 54/216 atlas taxa are direct dated-backbone tips.

The former raw absolute-slope lability relationship and median-split quadrants are retained only as statistical provenance/QA and are not part of the current biological story.

## Primary continuous endpoints

- orientation angle relative to EXIF-oriented image vertical;
- CIELAB lightness and chroma;
- circular hue sine/cosine;
- aspect ratio, circularity, solidity and width-profile variation.

The auxiliary high-resolution layer retains involucre projection roughness, outward spread fraction and maximum relative spine-like projection. These are image-geometry proxies, not direct botanical phyllary/spine measurements.

## Current submission gates

The computational taxonomic, niche-null and spatial-robustness gates are complete. Two scientific blockers remain external to repository computation:

1. adjudicated human boxes for the independent detector audit;
2. genuinely independent reference measurements for orientation, colour and outline.

Nomenclatural notes, authorship/administrative metadata and the durable archive DOI are finalization tasks after the measurement gates close.

## Reproducibility controls

```bash
python analysis/ch1/run_submission.py check
python analysis/ch1/run_submission.py claims
python analysis/ch1/run_submission.py summary
```

Every FDR count must name its cohort and multiplicity family. Circular hue components are interpreted jointly. Public photographs are not colour calibrated; image vertical is not gravity; outline measures remain viewpoint dependent; and grafted trees are historical sensitivity devices rather than a resolved species tree.

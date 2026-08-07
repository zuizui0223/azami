# External completion gates for Chapter 1

The remaining work is now split into **completed reproducible diagnostics** and **genuinely external/human validation**. Completed diagnostics must not be rerun merely to change manuscript framing; they are rerun only if a frozen input, taxonomic unit, measurement definition or accepted model changes.

## 1. Independent detector audit — MANUAL ANNOTATION REMAINS

The previous 1,000-image detector-development packet cannot serve as independent validation because it supplied the open-vocabulary proposals and 270 pseudo-labelled training/validation images for the recovered YOLO11n detector.

The leakage-free replacement packet and evaluation infrastructure are complete. The workflow:

- excludes every prior proposal/training `photo_id` and `obs_id`;
- freezes 1,000 taxonomically and geographically spread source images before detector evaluation;
- blinds taxonomy, coordinates and detector predictions;
- assigns 250 images to a second annotator;
- fixes the production weight SHA-256 and confidence threshold 0.25;
- generates hidden predictions down to confidence 0.01.

Precision, recall and F1 remain **unreported** until the two independent human annotation files are completed and disagreements are adjudicated. No repository-only analysis can substitute for those boxes.

## 2. Independent orientation, colour and outline validity — HUMAN/REFERENCE DATA REMAIN

Horizontal-mirror stability and production overlays are technical checks, not independent accuracy validation. The remaining gate requires reference measurements that are genuinely external to the production extraction:

- orientation: human landmarks/axes and camera-roll assessment;
- colour: calibration-supported references or repeat-image error with an explicitly justified reference design;
- outline: independent manual masks/contours with overlap and continuous-trait agreement.

Any material failure requires an analysis-version increment and propagation of measurement uncertainty into headline sensitivities. This gate cannot be closed by a second automated interpretation of the same images.

## 3. Authority-backed taxonomic freeze — AUTOMATED REVIEW PREP COMPLETE; 24 NAMES REQUIRE DECISION

A reproducible World Checklist of Vascular Plants (WCVP) review workflow is implemented in `analysis/build_wcvp_taxonomic_review_candidates.py` and `.github/workflows/ch1-wcvp-taxonomic-candidates.yml`.

Successful run `31152400481` reviewed the union of 259 frozen source names:

- 243 names had one exact WCVP canonical-name match;
- 16 names had multiple exact canonical-name matches, typically because accepted and homonymous/synonym records share the same author-free canonical name;
- 235 were routine accepted-name candidates;
- 8 were WCVP synonym candidates;
- 24 rows were flagged for high-priority human review;
- no source name was unmatched.

The eight WCVP synonym candidates map onto names that are already active analysis units, so blindly accepting those synonymies would merge analysis units and trigger a new analysis version. They must therefore be reviewed as taxonomic decisions rather than silently applied. The manuscript should continue to use *assigned taxon/species* language until this freeze is final.

## 4. Residual spatial and broad-region diagnostics — COMPLETED

The original compact SPDE workflow did not archive observation-level residuals, so a diagnostic-only refit/export was added without replacing accepted coefficient tables. Run `31152400475`, artifact `8983877726`, reproduced the frozen best-WAIC specifications closely and exported observation-level residuals.

Global k-nearest-neighbour Moran diagnostics now use spherical 3-D coordinates and a fixed neighbour graph across 999 permutations. Across all nine primary endpoints:

- Moran's I ranged from -0.0096 to 0.0030;
- no endpoint had permutation P < 0.05;
- minimum leave-one-broad-region-out species-rank Spearman rho ranged from 0.856 to 0.972.

Natural Earth continent assignments covered 46,274 of 46,276 source observations; two unmapped island observations were retained as an explicit diagnostic stratum. Full provenance and endpoint results are in `manuscript/results/SPATIAL_ROBUSTNESS_RESULTS_2026-08-07.md`.

**Decision:** the residual-spatial and broad-region robustness gate is closed for the frozen analysis. A new spatial diagnostic is required only if the accepted cohort, model specification or taxonomic units change.

## 5. Environmental-niche permutation and null tests — COMPLETED

Run `31152400487` executed 10,000 trait-label permutations while preserving the same complete species, environmental PCA and environmental availability. BH correction was applied across the nine primary traits within each metric and scope.

For all 148 complete species, significant environmental sorting was concentrated in orientation and visible colour; width-profile CV retained centroid-distance but not overlap support. In the 49-species direct-backbone sensitivity, support narrowed, with hue cosine most consistently retained.

Full results are frozen in `manuscript/results/NICHE_PERMUTATION_RESULTS_2026-08-07.md`.

**Decision:** the niche-null gate is closed for the current taxonomic units. It must be rerun if the taxonomic freeze merges or removes active units.

## 6. Authorship and administrative metadata — TEAM CONFIRMATION REMAINS

The final author order, affiliations, CRediT roles, acknowledgements, funding numbers and corresponding-author details require confirmation from the research team. These fields are intentionally not inferred from repository metadata.

## 7. Durable release — FINAL STEP AFTER HUMAN SCIENTIFIC GATES

After detector/measurement validation and taxonomic decisions are final, regenerate any sensitivity required by those decisions, freeze the manuscript, create an immutable release tag, produce SHA-256 manifests and deposit the licence-safe submission bundle in a durable repository with a persistent DOI.

## Current submission-blocking items

Only the following scientific items still require information that cannot be created honestly from repository evidence alone:

1. human boxes for the independent detector audit;
2. genuinely independent reference measurements for orientation, colour and outline;
3. review of the 24 high-priority WCVP taxonomic rows, especially the eight synonym conflicts that would merge active analysis units.

All other previously listed computational credibility gates are now executable and, for the frozen analysis, completed.

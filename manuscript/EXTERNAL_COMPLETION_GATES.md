# External completion gates for Chapter 1

The remaining work is split into **completed reproducible diagnostics** and **genuinely external/human validation**. Completed diagnostics are not rerun merely to change manuscript framing; they are rerun only if a frozen input, operational taxonomic unit, measurement definition or accepted model changes.

The current submission-facing story is global **continuous within-taxon image phenomics**: public photographs are converted into repeated numerical capitulum-trait measurements rather than categorical trait labels. The remaining gates therefore focus on whether those numerical measurements are independently valid.

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

Horizontal-mirror stability and production overlays are technical checks, not independent accuracy validation. The remaining gate requires reference measurements genuinely external to the production extraction:

- orientation: human landmarks/axes and camera-roll assessment;
- colour: calibration-supported references or repeat-image error with an explicitly justified reference design;
- outline: independent manual masks/contours with overlap and continuous-trait agreement.

This is the most important remaining scientific gate because the manuscript's methodological novelty depends on treating each assessable head as a numerical phenotype observation. Any material failure requires an analysis-version increment and propagation of measurement uncertainty into headline sensitivities. This gate cannot be closed by a second automated interpretation of the same images.

## 3. Taxonomic robustness — COMPLETED FOR THE OPERATIONAL-UNIT ANALYSIS

Chapter 1 uses exact source-platform taxon assignments as **operational taxonomic units** rather than claiming to resolve genus-wide species limits. A reproducible WCVP authority audit reviewed the union of 259 frozen names in run `31152400481`:

- 243 unique exact canonical-name matches;
- 16 multiple exact canonical-name matches;
- 235 routine accepted-name candidates;
- 8 WCVP synonym candidates;
- 24 high-priority nomenclatural-review rows;
- 0 unmatched names.

All eight synonym candidates map to names that are already separate active source-platform units. Run `31153385658` therefore collapsed all eight simultaneously as a sensitivity analysis rather than silently rewriting the primary taxonomy.

The collapse changed the balanced atlas from 216 to 211 units and the exhaustive primary from 259 to 251 units, but:

- all 8 original BH-supported primary component rows remained BH-supported;
- no coefficient changed sign;
- maximum absolute standardized-beta change was 0.000385;
- the minimum below-unit visible-variance fraction increased from 0.5886 to 0.5970;
- the largest absolute endpoint change in below-unit fraction was 0.0154.

**Decision:** taxonomic uncertainty does not explain the headline ecological result families. The scientific taxonomic-robustness gate is closed for the source-assigned operational-unit analysis. The 24 high-priority rows still require nomenclatural notes in the supplement, but they are no longer an untested source of the primary conclusions.

## 4. Residual spatial and broad-region diagnostics — COMPLETED

The original compact SPDE workflow did not archive observation-level residuals, so a diagnostic-only refit/export was added without replacing accepted coefficient tables. Run `31152400475`, artifact `8983877726`, reproduced the frozen best-WAIC specifications closely and exported observation-level residuals.

Global k-nearest-neighbour Moran diagnostics use spherical 3-D coordinates and a fixed neighbour graph across 999 permutations. Across all nine primary endpoints:

- Moran's I ranged from -0.0096 to 0.0030;
- no endpoint had permutation P < 0.05;
- minimum leave-one-broad-region-out taxon-rank Spearman rho ranged from 0.856 to 0.972.

Natural Earth continent assignments covered 46,274 of 46,276 source observations; two unmapped island observations were retained as an explicit diagnostic stratum.

**Decision:** the residual-spatial and broad-region robustness gate is closed for the frozen analysis. A new spatial diagnostic is required only if the accepted cohort or model specification changes materially.

## 5. Environmental-niche permutation and null tests — COMPLETED

Run `31152400487` executed 10,000 trait-label permutations while preserving the same complete taxa, environmental PCA and environmental availability. BH correction was applied across the nine primary traits within each metric and scope.

For all 148 complete taxa, significant environmental sorting was concentrated in orientation and visible colour; width-profile CV retained centroid-distance but not overlap support. In the 49-taxon direct-backbone sensitivity, support narrowed, with hue cosine most consistently retained.

**Decision:** the niche-null gate is closed for the current operational taxonomic units.

## 6. High-resolution involucre/spine-like layer — ANALYSIS COMPLETE, BOTANICAL VALIDATION DOWNSTREAM

The final manuscript uses three inferential auxiliary contour proxies:

- involucre projection roughness;
- involucre spread fraction;
- maximum spine-like projection relative to head radius.

All three are positively associated with BIO4 temperature seasonality in the ≤10 km auxiliary family after BH correction. `spine_peak_count_proxy` remains in integrated derived tables for provenance/downstream exploration but is not one of the three final manuscript auxiliary inferential endpoints.

These are image-geometry proxies, not direct botanical measurements of phyllary angle, actual spine length, orientation or stiffness. Their ecological interpretation is therefore bounded in Chapter 1. Botanical proxy→character validation and resolved phylogenetic mapping are handed to EAzami/focal studies rather than treated as additional Chapter 1 submission blockers.

## 7. Authorship and administrative metadata — TEAM CONFIRMATION REMAINS

The final author order, affiliations, CRediT roles, acknowledgements, funding numbers and corresponding-author details require confirmation from the research team. These fields are intentionally not inferred from repository metadata.

## 8. Durable release — FINAL STEP AFTER HUMAN MEASUREMENT GATES

After detector and continuous-measurement validation are final, regenerate any sensitivity required by a material audit failure, freeze the manuscript, create an immutable release tag, produce SHA-256 manifests and deposit the licence-safe submission bundle in a durable repository with a persistent DOI.

## Current submission-blocking scientific items

Only two scientific items still require information that cannot be created honestly from repository evidence alone:

1. human boxes for the independent detector audit;
2. genuinely independent reference measurements for orientation, colour and outline.

Nomenclatural annotations for the 24 high-priority WCVP rows and authorship/administrative fields still require human confirmation, but the computational spatial, niche-null and taxonomic-robustness gates are completed for the frozen operational-unit analysis.

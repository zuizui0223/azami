# External completion gates for Chapter 1

The remaining work is split into **completed reproducible diagnostics**, **image remeasurement**, and **genuinely external/human validation**. Chapter 1 is on submission hold. Completed diagnostics are not rerun merely to change manuscript framing; they are rerun only if a frozen input, operational taxonomic unit, measurement definition or accepted model changes.

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

## 6. High-resolution involucre/spine-like layer — HEADLINE WITHDRAWN

The final manuscript uses three inferential auxiliary contour proxies:

- involucre projection roughness;
- involucre spread fraction;
- maximum spine-like projection relative to head radius.

All three were positively associated with BIO4 temperature seasonality in the original ≤10-km auxiliary family. After the predeclared minimum-dimension and sharpness adjustment, none passed BH correction (q = 0.0696–0.0730), and projection roughness and maximum spine-like projection reversed sign in the 150–199-pixel stratum. The three rows are withdrawn from the biological headline and retained only as provenance.

These are image-geometry proxies, not direct botanical measurements of phyllary angle, actual spine length, orientation or stiffness. A later direct botanical study may revisit them, but Chapter 1 does not claim an involucral seasonality result.

## 7. Developmental stage — IMAGE LABELLING REMAINS

All 46,276 strict-cohort rows have parseable observation dates. Raw-calendar taxon-specific sine/cosine adjustment retained all four non-circular primary rows, but three dominant taxa were sampled in both hemispheres. A second contract was therefore frozen before outcome inspection. It added a half-cycle to Southern Hemisphere dates or fitted separate taxon-specific cyclic curves by hemisphere. All four rows retained their frozen sign, BH support and omission-sign stability under both definitions. These checks resolve the hemispheric calendar issue but do not distinguish bud, anthesis and post-anthesis heads. A stage-labelled or anthesis-restricted sensitivity remains submission-blocking for the retained orientation and outline claims.

## 8. Visible-colour negative control — IMAGE REMEASUREMENT REMAINS

The paired flower-versus-background `region × BIO12` implementation is tested and fixed. It measures outer-context and nonfloral-green CIELAB controls and does not infer specificity from a nonsignificant background slope. A self-contained 1,413-head context-crop artifact remains available, but the full production workflow removed temporary context crops before uploading its numerical outputs. The main 46,276-observation sensitivity therefore requires inserting the fixed control measurement into the existing image-reconstruction workflow. Chroma–BIO12 was independently withdrawn from the current headline after failing native-only BH correction; this pass remains necessary for broader corolla-specific colour interpretation.

## 9. Repeat photographs — COHORT FROZEN, MEASUREMENTS REMAIN

An outcome-blind metadata preflight identified 20,073 strict-cohort observations from 236 taxa with at least two public photographs: 58,748 photo rows and 38,675 photographs beyond the first. Trait remeasurement is pending. Even after remeasurement, between-photo variance includes camera, illumination, viewpoint and possibly subject differences unless manual review verifies the same individual.

## 10. Native-range sensitivity — COMPLETE, TWO OF FOUR ROWS RETAINED

The locked analysis pinned WCVP checklist DOI 10.15468/6h8ucr and the exact TDWG level-3 geometry commit and content hash. It resolved 245 of 259 source taxa and classified 27,066 of 46,276 observations as native; ambiguous accepted-name joins, unlisted units and unmapped coordinates were excluded fail-closed. Native-only refits retained orientation–BIO1 and aspect ratio–BIO4. Chroma–BIO12 and aspect ratio–BIO12 retained their primary signs but failed native-only BH correction and were withdrawn from the current headline. The analysis does not infer introduction history or adaptation.

## 11. Authorship and administrative metadata — TEAM CONFIRMATION REMAINS

The final author order, affiliations, CRediT roles, acknowledgements, funding numbers and corresponding-author details require confirmation from the research team. These fields are intentionally not inferred from repository metadata.

## 12. Durable release — FINAL STEP AFTER SCIENTIFIC GATES

After the scientific gates are final, regenerate any sensitivity required by a material audit failure, freeze the manuscript, select and record an appropriate code/data licence, create an immutable release tag, produce SHA-256 manifests and deposit the licence-safe submission bundle in a durable repository with a persistent DOI. GitHub Actions artifacts are time-limited and are not the durable archive.

## Current submission-blocking scientific items

The current scientific blockers are:

1. human boxes for the independent detector audit;
2. genuinely independent reference measurements for orientation, colour and outline;
3. developmental-stage labels or an anthesis-restricted sensitivity;
4. paired visible-colour negative controls;
5. repeat-photo trait remeasurement and bounded variance decomposition.

Nomenclatural annotations for the 24 high-priority WCVP rows, licence selection and authorship/administrative fields still require human confirmation. The computational spatial, niche-null, taxonomic, raw-calendar and hemisphere-aware seasonal, dominant-taxon, native-range and involucre image-quality audits are complete for the frozen operational-unit analysis.

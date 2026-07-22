# External completion gates for Chapter 1

The repository now contains the reviewer-driven statistical corrections and the reproducible infrastructure for the remaining submission gates. Gates must be closed in the order below because later analyses depend on the accepted measurement and taxonomic units.

## 1. Independent detector audit — packet and evaluation implemented; manual boxes pending

The previous 1,000-image detector-development packet cannot serve as independent validation because it supplied the open-vocabulary proposals and 270 pseudo-labelled training/validation images for the recovered YOLO11n detector.

The replacement workflow `.github/workflows/ch1-independent-detector-audit-packet.yml`:

- excludes every prior proposal/training `photo_id` and `obs_id`;
- selects a new 1,000-image taxonomically and geographically spread source-image sample before detector inference;
- blinds taxonomy, coordinates and detector predictions;
- assigns 25% of images to a second annotator;
- fixes the production weight SHA-256 and production confidence threshold 0.25;
- generates hidden predictions down to confidence 0.01 for later diagnostics.

The local annotation, adjudication and evaluation programs are implemented in `analysis/`. Precision, recall and F1 remain **unreported** until two human annotation files and third-party adjudication of disagreements are complete. See `DETECTOR_INDEPENDENT_AUDIT_PROTOCOL.md`.

## 2. Independent orientation, colour and outline validity

Horizontal-mirror stability and production overlays are technical checks, not accuracy validation. The next gate requires independent reference measurements for each module:

- orientation: human landmarks/axes and camera-roll assessment;
- colour: controlled or calibration-supported reference subsets and repeat-image error;
- outline: independent manual masks/contours with overlap and continuous-trait agreement.

The validation sample must remain separate from threshold development. Error distributions must be propagated to headline variance and association sensitivities.

## 3. Authority-backed taxonomic freeze

The accepted-name decision table must be populated and externally reviewed. Automated matching alone is insufficient for synonyms, infraspecific ranks, hybrids, approximate matches and taxonomically difficult Cardueae names.

Required input: one citable authority-backed decision per frozen source name, including the number of records affected. A decision that merges, removes or reassigns an active analysis unit triggers a new analysis version and rerun.

## 4. Residual spatial and broad-region diagnostics

The audit and preparation code exists, but the grouped-SPDE result bundle still needs an observation-level fitted/residual export and a reviewed observation-to-region lookup.

Required outputs:

- endpoint-specific residual Moran's I;
- regional observation and taxon coverage;
- leave-one-region-out stability;
- dominant-region and spatial-mesh/prior sensitivity where supported by the frozen model outputs.

## 5. Environmental-niche permutation and null tests

The trait-extreme niche contrasts must be compared with predeclared null distributions that preserve the relevant taxonomic sample size and geographic availability. Report uncertainty, multiplicity correction and sensitivity to direct-backbone/taxonomic restrictions. Descriptive centroid distance or overlap alone is not inferential support.

## 6. Authorship and administrative metadata

The final author order, affiliations, CRediT roles, acknowledgements, funding numbers and corresponding-author details require confirmation from the research team.

## 7. Durable release

After all scientific gates and manuscript revisions are frozen, deposit the checksummed submission bundle and add the immutable release tag, final commit and archive DOI.

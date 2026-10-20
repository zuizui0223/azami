# Original-source v4 execution — 2026-10-20

## What was actually executed

The historical pre-detector archive was recovered from its SHA-identical owner backup and scanned: **1,122,854 photo rows / 665,115 observations** (2026-07-03 snapshot; not the current entire iNaturalist database). No v3 phenotype table or detection-positive eligibility was used.

A redistributable diagnostic lane contains **131,637 photographs / 72,602 observations / 214 source taxa** after species-rank, noncaptive and CC0/BY/BY-SA restrictions. Interleaved geographic cells and a cap of ten observations per taxon give **1,409 queued photos**. The licence and cap restrictions are diagnostic design choices, not a final ecological universe. The original master frame is retained, and native-primary status is a separate upstream join.

## Actual original photos, not thumbnails

GitHub source-acquisition run **39371939888** succeeded and preserved artifact **12241906850** (ZIP SHA256 `5600dbe9db39dddbf0f74a5ba9534d2abfd865be0c3071f781c1c0a4efc7aece`). It retrieved original bytes for 64 photos from 32 source taxa. Initially 63 decoded, with one ICC transform failure.

The remaining file was a real **monochrome photograph with a GRAY ICC profile**, not an unavailable photo. The decoder had converted it to RGB before applying the GRAY transform. The corrected decoder keeps it available for geometry and explicitly disallows colour inference. Replay of the exact bytes now yields **64 geometry-decodable images, 63 source-colour candidates, one monochrome image**. This is not a statement that all 63 show assessable flowers. No altered byte or replacement photo is substituted.

## Detector rerun (local CPU)

The recovered frozen model weight was independently checksum-verified (`4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed`). A fresh run on all 64 originals, Ultralytics 8.4.183 / torch 2.10.0+cpu, imgsz 960, confidence 0.25, IoU 0.7, returned **161 candidate boxes in 61 photographs**. All no-detection cases remain in the ledger. These are localization outputs, **not 161 verified capitula and not 61/64 recall**. No v3 bounding boxes were reused. This detector run was local; remote acquisition/re-measurement CI does not imply independent detector accuracy.

## Direct-pixel/anatomical proof of execution

Five convenience examples were annotated after inspection, labelled **assistant_visual_proposal** throughout. Four contain floral display patches measured directly in image Lab; the monochrome example is explicitly colour-unassessable. Three have proposed head/peduncle landmarks, and two have explicit involucre-body masks. Missing landmarks and unannotated phyllary tips are unavailable, not biological zeros.

The colour masks are **sample patches of visible floral display**, not complete botanically resolved corolla segmentation. They may include styles or other floral organs. Their role is to test the new measurement contract and expose the underlying pixels, not support between-species colour comparisons. Proposed annotation and image hashes, acquisition attribution and measured pixels travel together. Image L*a*b* is not calibrated reflectance or anthocyanin content.

The source and anatomy tests cover white-pixel retention, unavailable hue, original URL checks, deterministic sampling, wrong source hashes, malformed URLs, invalid ICC, monochrome exclusion, absent anatomy and exclusion masks. The new GitHub anatomical workflow verifies the original source ZIP and re-runs the measurements; no frozen v3 phenotype CSV is supplied.

## Scientific state

- Complete: metadata recovery/scan; bounded original acquisition; local detector rerun; direct-pixel/geometric measurement implementation and five-example execution.
- Still required: independently checked anatomical segmentation at scale; repeated-image and exposure errors; native-primary linkage; new environment extraction and spatially held-out multivariate models; within/among covariance from those **new** observations.
- The earlier 42-taxon R²/importance numbers remain downstream exploratory diagnostics, not newly validated upstream results.
- Original source artifact expiry is recorded by GitHub (2026-11-19). Source photos/model also have separate owner-side provenance; this new source bundle has not yet been claimed as a permanent public archive.

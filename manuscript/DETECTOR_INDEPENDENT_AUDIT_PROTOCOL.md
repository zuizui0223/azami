# Independent visible-capitulum detector audit protocol

## Decision and scope

The production detector is evaluated as a single-class `visible_capitulum` object detector. The audit concerns only whether every human-auditable visible capitulum is detected and whether detector boxes correspond to real capitula. It does not validate orientation, colour or outline measurements; those are the next submission gate.

The frozen production inference threshold is **0.25** and the predeclared object-matching threshold is **IoU = 0.50**. Detector predictions are generated down to confidence 0.01 so the production threshold can be evaluated without rerunning the model and a diagnostic precision–recall curve can be reported without selecting a new threshold on the audit labels.

## Why the earlier 1,000-image packet is not independent

The earlier detector-audit source packet (artifact `8068122589`) supplied the images used by the open-vocabulary proposal workflow. The recovered detector package (artifact `8076736948`) shows that 270 of those images became pseudo-labelled training/validation images for the YOLO11n bootstrap detector. Its training provenance explicitly states that the reported training metrics quantify agreement with pseudo-labels and are not independent accuracy estimates.

Accordingly, the earlier 1,000 images are retained as detector-development provenance but are prohibited from the independent evaluation. The new selector excludes **both `photo_id` and `obs_id`** from every prior proposal/training source queue. Excluding observations as well as photos prevents another photograph from the same public observation entering the evaluation.

## Frozen detector identity

- recovery artifact: `8076736948`;
- initialization: `yolo11n.pt`;
- target recorded by training workflow: `visible_capitulum`;
- training labels: open-vocabulary pseudo-label proposals;
- pseudo-labelled images: 270 (211 train, 59 validation);
- input size: 640 pixels;
- training epochs: 40;
- production weight SHA-256: `4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed`;
- production confidence threshold: 0.25.

The audit workflow refuses a detector package whose weight digest differs from this value.

## Independent source-image sample

`analysis/build_independent_detector_audit_packet.py` selects the audit before detector inference from the frozen merged metadata artifact `8066010557`.

Eligibility:

- species-rank public observations;
- non-captive records;
- at least one downloadable image rendition;
- one photograph per observation;
- no `photo_id` or `obs_id` appearing in the prior proposal/training source set.

The default sample contains 1,000 source images. Every eligible remaining species is represented before additional images are allocated. Within species, distinct 10-degree spatial blocks are prioritized. The deterministic selection seed is `20260722`; 25% of images are independently annotated twice.

The selector writes two strictly separated products:

1. **blinded annotation packet** — audit ID, image filename and annotation assignment only;
2. **private key** — taxon, coordinates, observation/photo IDs, spatial block, latitude band, user and licence metadata.

Taxon, coordinates and detector predictions are never shown to annotators.

## Biological object definition

A `visible_capitulum` is one discrete thistle capitulum that a human can bound in the source image. Anthesis-stage, bud and post-anthesis capitula are recorded separately for stratified performance summaries but use the same detector target when their boundary is visually identifiable.

Annotators must:

- draw one tight source-pixel box around every visible target;
- record `no_target` when no capitulum is visible;
- use `unassessable` only when the image cannot be exhaustively audited;
- record life stage, occlusion, edge truncation and overall image quality;
- avoid using taxon identity or detector suggestions.

## Human annotation and adjudication

The local application `analysis/detector_audit_annotation_app.py` displays only blinded source images. Annotator 1 completes all 1,000 images; annotator 2 completes the preassigned 25% subset.

`analysis/finalize_detector_audit_annotations.py` validates all assignments and compares double annotations using one-to-one matching at IoU 0.50. A double-labelled image is accepted automatically only when:

- assessability states agree;
- object counts agree; and
- every box can be paired at IoU >= 0.50.

All disagreements are exported to `detector_adjudication_required.csv` and require a third blinded adjudicator. The canonical ground-truth table is not written as complete while unresolved images remain.

## Detector evaluation

`analysis/evaluate_independent_detector_audit.py` evaluates the frozen predictions at confidence 0.25 and IoU 0.50. Main outputs are:

- object-level precision, recall and F1;
- Wilson 95% intervals for precision and recall;
- source-image bootstrap 95% intervals;
- true no-target image specificity;
- positive-image sensitivity;
- false-positive and false-negative review queues;
- performance by ground-truth box size, life stage, occlusion, edge truncation and image quality;
- post-unblinding summaries by latitude band and spatial block, with taxon summaries only above a declared minimum denominator;
- a diagnostic confidence curve from 0.05 to 0.95.

The diagnostic curve is descriptive. The production threshold remains 0.25 unless a new detector version is trained and the complete downstream image analysis is rerun.

## Interpretation and stop rule

No detector precision, recall or F1 value may enter the manuscript until the independent annotations are complete and adjudicated. The production-output Figure 1 is a provenance demonstration, not detector validation.

If the audit shows material failure overall or in biologically relevant strata—particularly small/distant, occluded or unusual-view capitula—the existing downstream dataset is not silently repaired. The detector, measurement layer and all dependent analyses must receive a new analysis version and be rerun from source photographs.

## Reproducible entry points

1. Generate the leakage-free packet and hidden predictions:
   `.github/workflows/ch1-independent-detector-audit-packet.yml`.
2. Annotate the downloaded blinded packet twice.
3. Finalize and adjudicate:

```bash
python analysis/finalize_detector_audit_annotations.py \
  --manifest detector_independent_audit_blinded_manifest.csv \
  --assignments detector_independent_audit_assignments.csv \
  --annotations annotator_1.csv \
  --annotations annotator_2.csv \
  --adjudication adjudicator.csv \
  --out-dir detector_ground_truth
```

4. Evaluate the frozen predictions:

```bash
python analysis/evaluate_independent_detector_audit.py \
  --annotations detector_ground_truth/detector_independent_audit_annotations.csv \
  --predictions yolo_crop_metadata.csv \
  --private-key detector_independent_audit_private_key.csv \
  --production-confidence 0.25 \
  --minimum-prediction-confidence 0.01 \
  --iou-threshold 0.50 \
  --bootstrap-repeats 2000 \
  --out-dir detector_independent_evaluation
```

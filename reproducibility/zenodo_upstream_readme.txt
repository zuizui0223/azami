Azami: source-to-analysis preservation map

The starting point is public photographs, not a finished trait table.
This record preserves historical processing artifacts as well as current inputs
and results. It is a draft. Upload completion is not permission to publish.

1. OBSERVATIONS AND PHOTOGRAPH SOURCE
8066010557: full iNaturalist metadata inventory; retains original photo/source
identifiers and source metadata. 8066675131: original screening queue.
8068122589: original medium-image screening/QC package. This package was used in
detector development and must NOT be treated as independent detector validation.
The complete mutable original-image collection is not preserved by these files.

2. PSEUDO-LABELS AND YOLO TRAINING
8069610715 / 8071529579: earlier Grounding DINO proposal and bootstrap dataset
artifacts. 8076736948: recovered frozen production detector package including
best.pt, last.pt, dataset manifests, pseudo-box proposals, training arguments,
results and recovery provenance. best.pt SHA-256:
4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed
The archived training manifest records 270 pseudo-labelled images: 211 training,
59 validation. Agreement with pseudo-labels is not independent accuracy.
The model package alone does not contain all training image/YOLO label files.
Reconstruction from the QC image package and pseudo-box/manifest records has
not yet been executed and verified. Do not replace that test with file presence.

3. HISTORICAL CATEGORICAL MACHINE LEARNING
8077189280: CLIP zero-shot outputs, ontology, prompts and model provenance.
8099953404: merged historical ensemble outputs and crop metadata.
Recorded ensemble revisions:
openai/clip-vit-base-patch32 3d74acf9a28c67741b2f4f2ea7635f0aaf6f0268
openai/clip-vit-base-patch16 57c216476eefef5ab752ec549e440a49ae4ae5f3
Grounding DINO is identified as IDEA-Research/grounding-dino-base in the proposal
provenance. External pretrained weights themselves are not archived here.
These categorical model outputs are historical; the current numerical analysis
does not substitute their class scores for deterministic continuous traits.

4. DETECTION, CROPS AND CONTINUOUS MEASUREMENTS
8225059018: early continuous-trait/provenance source.
8269246732: exhaustive detector-positive merged measurement output, preserved
as its exact original ZIP. Individual source rows and processing provenance are
retained; it must not be described as a full archive of original photographs.
9612943217: continuous endpoint universe subsequently used by the current paper.

5. ENVIRONMENT AND CURRENT ECOLOGICAL ANALYSIS
The analysis-data ZIP contains the continuous universe, environment extract,
region lookup, placement trees, native-status input and current result artifacts.
ARTIFACT_CATALOG.json maps these to the current seven-stage numerical workflow.
The earlier endpoint-only baseline remains in DOI 10.5281/zenodo.22295791.
All original artifact bytes are preserved, including unsupported results and
historical fields. Current reference manifests distinguish active evidence.

6. INDEPENDENT AUDIT MATERIALS: INCOMPLETE AND NON-PUBLIC
8521924881 / 8521925441 / 8521926057: independent annotation packet, private audit
mapping/QC and hidden predictions. These are materials for an unfinished audit,
not completed independent precision/recall or biological validation.
Keep the private mapping and predictions protected while blinding matters.
Image-specific redistribution rights and audit disclosure require review before
publishing any of these files. Saving them in the owner draft does not approve
public release. Use a separate restricted record or omit protected items from a
public version if needed; do not simply publish all draft files together.

CODE AND REPRODUCTION STATUS
Current code commit: 52c207a64a09b243e021dc6e7b20598ab4ab7f02
Scientific reference: fe25abd46e7235c85e6da191f976e1c8a02d0406
Historical acquisition/training/measurement scripts and ontology are under
legacy/ch1_global/v2/ in that snapshot. Current analyses start in analysis/v3/.
The source workflow commit for each artifact is recorded in its catalog.
No retraining, fresh image acquisition or new ecological fitting is performed
by the archive-transfer job. Zero-from-scratch reproduction is NOT yet certified.

PENDING BEFORE CLAIMING COMPLETE FROM-ZERO REPRODUCTION
- Verify recovery of training images and labels and reproduce the training setup.
- Verify availability and licenses of exact external pretrained model revisions.
- Verify coverage and access of the original photograph collection; URLs can die.
- Review image, model, backbone and private audit-material distribution rights.
- After explicit publication approval, anonymously download and verify the final
  intended public/restricted data partition. A private readback is not that test.

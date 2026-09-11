Data dependencies of the current Azami Chapter 1 manuscript

SCOPE
This is the manuscript's source-to-result chain, not a repository-history backup.
Code is referenced at 52c207a64a09b243e021dc6e7b20598ab4ab7f02 on GitHub.
Scientific reference: fe25abd46e7235c85e6da191f976e1c8a02d0406.
Legacy location alone neither admits nor excludes a file: only its role in the
current Main or Supporting Information determines whether it is needed.

1. PHOTOGRAPHS USED TO BUILD THE STUDY COHORT (MAIN METHODS; SI S4)
8066010557: iNaturalist metadata inventory including photograph identifiers and
source information. 8269246732: exhaustive detector-positive merged processing
output underlying the measured universe and subsequent filtering/thinning.
These are not an archive of every original photograph. Original URL availability
and photo-specific licenses remain distinct from preserved numerical data.

2. THE ADOPTED HEAD DETECTOR (MAIN METHODS; SI S1)
8076736948: adopted YOLO11n model, best.pt/last.pt, pseudo-box proposals, training
and validation manifests, training settings, diagnostics and recovery provenance.
best.pt SHA-256:
4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed
8066675131 and 8068122589: source screening queue and image/QC package needed to
trace the adopted detector's input data. They are not independent validation.
The model's training record has 270 pseudo-labelled images (211 training, 59
validation); 40 epochs, image size 640, batch 8, CPU, seed 0, pretrained yolo11n.pt.
The package records Grounding DINO proposal generation; it does not preserve a
verified exact Grounding DINO revision or initial pretrained YOLO checkpoint.
Retraining from the archived images/manifests has not been verified.

3. HEAD EXTRACTION, CONTINUOUS MEASUREMENTS AND TECHNICAL CHECKS
8225059018: continuous-trait and figure measurement provenance.
8099953404: 6,626-head source/crop records relevant to the mirror technical check
in SI S2 and the historical figure-measurement provenance. Its accompanying old
categorical predictions are incidental retained bytes, not current evidence.
8269246732: full merged continuous measurement output.
9612943217 (inside analysis-data ZIP): endpoint universe used by the paper.
The manuscript's biological constructs derive from continuous measurements;
retired MobileNet classes and CLIP category scores are not their inputs.

4. ENVIRONMENT, SENSITIVITIES AND WHOLE-CAPITULUM RESULTS
The analysis-data ZIP supplies the exact nine-predictor environment, broad-region
lookup, placement trees, native-status sensitivity table and current result
artifacts. ARTIFACT_CATALOG.json maps source runs/commits and member checksums.
The seven-stage GitHub runner compares 15 current reference aggregates.
Unsupported rows remain present. Superseded report fields are explicitly marked
as historical and are not restored as manuscript evidence.

EXCLUDED
Unrelated projects, abandoned branches, exploratory model development, retired
categorical ML pipelines, earlier failed training exports, and incomplete private
independent-audit annotation packets/mappings/predictions are not required to
reproduce this manuscript and are excluded from this preservation pass.

INTEGRITY AND LIMITS
Original artifact ZIPs are preserved unchanged; hence a necessary artifact can
contain incidental historical files. Use the stated role, not every embedded
file, as the manuscript dependency definition. Catalogs record every member hash.
This is a data-preservation operation, not a new analysis or accuracy evaluation.
The full original photograph collection, initial external pretrained weights,
and end-to-end training replay are not certified complete. Data and image rights
must be confirmed before publication. No publication is performed by this job.

Azami Chapter 1: numerical data and GitHub Actions artifacts
Code: https://github.com/zuizui0223/azami/tree/52c207a64a09b243e021dc6e7b20598ab4ab7f02
Scientific reference: fe25abd46e7235c85e6da191f976e1c8a02d0406

PURPOSE
Preserve the exact Actions artifacts underlying the current paper independently
of Actions retention limits. Code remains on GitHub; no code snapshot, manuscript,
reviewer correspondence, or original photograph collection is added here.
The archive preserves numerical measurements, not a repeat of image acquisition,
YOLO training or an independent biological accuracy evaluation.

CONTENTS
Four input artifacts: continuous endpoint measurements (9612943217), nine-predictor
environment table (9633419268), broad-region lookup (8983877726), and phylogenetic
placement trees (8227254443). The source cohort includes native and introduced
records; native status is a sensitivity input, not the primary filter.
Four current result artifacts: biological constructs and sensitivity chain
(10130210432), scale integration (10131007603), complete-18 integration plus direct
scale contrast (10136229131), assessability and technical stress (10135679053).
One provenance artifact: the earlier complete-18 upgrade (10135139012), retained
as an earlier stage, not a competing final result.
The native-status CSV is supplied separately from its immutable Git tag.
Source contracts and the 15-file current reference manifest are metadata, not new
analyses. Every embedded ZIP is preserved unchanged with its original digest.

CURRENT VERSUS HISTORICAL CONTENT
Input artifacts can contain older results in addition to the exact inputs used
now. Their inclusion is provenance preservation, not endorsement as current
results. Use ARTIFACT_CATALOG.json and the current reference manifest to select
current evidence. Retired fields environment_signature_alignment and
integration_environment_coupling remain inside original artifact bytes but are
not part of current scientific claims. Unsupported results are not removed.
The earlier v2 deposit 10.5281/zenodo.22295791 remains unchanged and separately
identifies the earlier endpoint baseline. This collection is not every CI run,
pilot, failed job, original image, or upstream raw environmental raster.

USING THE INPUTS WITH THE GITHUB RUNNER
At the fixed code commit, copy the four input ZIPs into work/current/archives/,
renaming input_continuous to continuous, input_environment to environment,
input_spatial to spatial and input_historical to historical in their filenames.
Copy supplemental_inputs/native_status.csv into work/current/inputs/.
Install the pinned runtime using reproducibility/requirements-current.txt and
run python -m reproducibility.run_current_analysis without --download.
Numerical data are supplied here; no GitHub authentication is needed for them.
See the fixed GitHub CURRENT_ANALYSIS.md for the seven stages and validation.
The runner verifies extracted input hashes and compares 15 current aggregates.

SOURCES AND INTERPRETATION
Continuous measurements derive from public iNaturalist observations and use the
image-based definitions and units present in the measurement table. They are not
direct physical colour, pigment concentration, or gravity-referenced angles.
Environmental values are the frozen CHELSA v2.1/BIOCLIM+ 1981-2010 extract,
not newly computed exposures. Preserve stored encodings; source metadata lists
variables and upstream URLs. Do not infer physical units from column names.
Native-status metadata derives from WCVP/TDWG regional status; source contract
and exact native-status bytes are included. Unknown statuses remain unknown.
Placement trees are hypotheses derived from GBOTB.extended.LCVP/V.PhyloMaker2;
52 placements do not mean 52 independent phylogenies. S1 and S3 are identical
in the frozen inputs. Preserve those original bytes.

RIGHTS AND PUBLICATION
This is an unpublished draft until explicitly published in Zenodo.
Third-party data retain their source terms; the software MIT license does not
relicense this collection. CHELSA v2.1 source page lists CC0; Natural Earth map
data are public domain; the frozen WCVP contract records CC BY 4.0. V.PhyloMaker2
package is GPL-2, but underlying backbone-data redistribution terms and final
record-level license scope still require confirmation before publication.
Existing record license metadata must not be read as blanket new clearance.
Do not publish merely because upload verification passed.

INTEGRITY
SHA256SUMS.txt verifies the outer ZIP. MANIFEST.json verifies each payload.
ARTIFACT_CATALOG.json links artifact/run/commit, original ZIP hashes and every
member hash. Upload verification is not a new numerical or biological result.

# Chapter 1 Zenodo publication handoff — 2026-09-27

## Status

Repository-side preparation is complete.

The verified public-upload candidate is:

- file: `azami_ch1_v3_analysis_inputs.zip`
- SHA-256: `be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438`
- size: 59,461,974 bytes
- GitHub code commit paired to the data: `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f`
- build workflow: `36288310738`
- Actions artifact: `10921221306`
- durable owner copy of the workflow wrapper: Drive file `16pildr5nOL08HAb6ATEHepUPei5oEn0t`
- release manifest: `release_ready=true`, `release_gaps=[]`
- public upload performed: **no**

## Publish as

Create a **new version under the existing Zenodo concept**, not a new unrelated record and not an overwrite of v2.

Existing records:

- concept DOI: `10.5281/zenodo.22295790`
- published v2 record DOI: `10.5281/zenodo.22295791`

Recommended metadata:

- title: **Azami Chapter 1 v3 numerical analysis input package**
- resource type: **Dataset**
- creator: **ZHANG, Ruiqi**
- ORCID: leave blank unless explicitly supplied
- affiliation: leave blank unless explicitly supplied
- record license: **CC BY 4.0**
- archive strategy: **new version of the existing concept**

Description:

> SHA-locked processed numerical inputs for reproducing the current Azami Chapter 1 v3 analysis with the separately versioned GitHub code. The archive contains only analysis inputs: continuous image-derived measurements, the nine-predictor environment input, broad-region lookup, historical-placement trees, and the frozen observation native-status table. Analysis code, manuscripts, figures, fitted/reference outputs, replay receipts, original third-party photographs, and detector training material are intentionally excluded.

GitHub code reference:

- repository: `https://github.com/zuizui0223/azami`
- pinned commit: `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f`
- runbook: `reproducibility/CURRENT_ANALYSIS.md`

## Exact public payload boundary

The uploaded ZIP contains exactly ten files:

1. `inputs/artifact-9612943217-continuous.zip`
2. `inputs/artifact-9633419268-environment.zip`
3. `inputs/artifact-8983877726-spatial.zip`
4. `inputs/artifact-8227254443-historical.zip`
5. `inputs/observation_native_status.csv`
6. `input_contract.json`
7. `README.txt`
8. `checksums.json`
9. `release_manifest.json`
10. `zenodo_release_metadata.json`

Do **not** add:

- repository/code ZIPs;
- manuscript/SI files;
- figures;
- current reference outputs;
- replay receipts;
- fitted outputs;
- raw photographs;
- detector training material.

Zenodo is only the durable byte store for replay inputs. GitHub remains the executable analysis record.

## Pre-upload verification already completed

The final candidate was independently unpacked after CI.

Verified:

- outer upload-file SHA-256 matches the sidecar;
- all nine rows recorded in the embedded `checksums.json` match their actual files;
- `analysis_input_count=5`;
- `code_included=false`;
- `manuscript_files_included=false`;
- `figures_included=false`;
- `reference_outputs_included=false`;
- `replay_receipts_included=false`;
- final metadata has `release_approved=true`;
- final metadata is bound to GitHub commit `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f`.

## Required post-publication verification

Do not mark current v3 public reproduction complete merely because Zenodo accepts the upload.

After publication:

1. record the new version DOI;
2. open/download the record without owner credentials;
3. download `azami_ch1_v3_analysis_inputs.zip`;
4. verify SHA-256 equals `be57e9ba80773c97afdc200f9a3982eaa85044cfc6791b62187ccfcd41f54438`;
5. verify the embedded `checksums.json`;
6. check out GitHub commit `7483ae8e20fa5df2439ad524312d3ac9a0fbf71f`;
7. run the current eight-stage numerical replay using the downloaded inputs;
8. confirm the 16-file current-reference comparison passes;
9. record the public DOI and verification receipt in the repository.

Only after step 8 should the repository status change from **release-ready / not public** to **credential-free public reproduction verified**.

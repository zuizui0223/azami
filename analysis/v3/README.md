# Chapter 1 v3: retain the source information, then define analysis views

V3 starts with the recovered acquisition snapshot, not the 46,276 observations
already selected for v2. The current entry point is
`python -m analysis.v3.workflow inventory`.

The archived July 2026 photo metadata contain **665,115 observations and
1,122,854 photos**. These are records of photos, not proof that every image file
has been downloaded. The snapshot is not a count of all iNaturalist records today.

Executed evidence: [source inventory receipt](../../reproducibility/v3_source_inventory_20260907.json),
[original-archive recovery](../../reproducibility/v3_upstream_recovery_20260907.json), and
[historical processing recount](../../reproducibility/v3_historical_processing_20260907.json).
The full snapshot was inventoried twice with byte-identical copies of all four
outputs; no source row was removed. The 46,276 legacy observations all matched.

The later [pre-merge/API reconciliation](../../reproducibility/v3_source_reconciliation_20260907.json)
restored **47 observation-photo links** missing from the merged snapshot:
one lost during merging and 46 recoverable from archived API records. The local
source-link ledger now covers **665,139 observation IDs, 1,122,854 unique photo IDs
and 1,122,901 unique links**. These are not 47 new photos: 44 photos have links to
multiple observations. Keep these links and flag shared-photo dependence before
evaluation splits or analysis; do not count them as independent image evidence.

## One source, reversible analysis views

```text
Acquisition snapshot: every observation-photo record
  → Rights-aware image cache: exact bytes, hashes, missing/restricted states
  → Detection: all heads plus explicit negative and failed photo jobs
  → Measurement: per-head/photo features, masks, QC and uncertainty
  → Analysis views: observation-level or nested, question-specific membership
  → Shared result tables: methodological and ecological findings
```

Keep all available source records and relationships. Native range, coordinate
quality, taxonomic rank, captive status, flowering state and image quality are
annotations, not reasons to delete records from the master data. A photo can be
useful for one endpoint and not another. Coordinates are needed for exposure
alignment, not for every image measurement. Restricted locations must not be
de-obscured; metadata retention does not authorize every image download or redistribution.

Multiple photos and multiple heads are retained, but they are not independent
observation replicates. One observation can show different heads or individuals.
A same-head match is needed before interpreting repeated photos as photographic
repeat measurements. Preserve raw values before deriving summaries, and use
appropriate nesting/weighting so photo-rich observations do not get unintended
extra weight.

The first-photo queue in `05_build_image_screening_queue.py` and the later
all-photo queue in `71_build_within_species_expansion_queue.py` are different
historical routes. Do not assume the first-photo filter applied to all v2 data.
V3 records which route each source followed. Spatial thinning becomes a justified
analysis view, not irreversible source reduction.

## Implementation and present boundary

| Stage | Implemented now | Still to implement |
|---|---|---|
| Inventory | Every merged/source row and archived API link accounted for; lost links restored; source hashes, observation aggregation and v2 overlap | First chunk's raw API file is absent; complete collection-time recovery cannot be certified |
| Recover | Six original metadata archives, the unthinned processing archive and verified local image bytes with all source links | Full-source image availability and a rights-aware acquisition/retention schedule; processing records are not image files |
| Measure | Historical five-field recovery; resumable cached-image detection and all-27 baseline extraction with raw/QC values and paired context | Full-source execution, crop/resolution/photometric evaluation and propagation of measurement uncertainty |
| Analyse | Workflow specification; no v3 ecological fitting | Save exact estimands, cohort memberships, formulas, uncertainty, spatial/nesting design and multiplicity before fitting |
| Report | Aggregate inventory, historical recovery and cached-execution receipts | Complete technical-evaluation and ecological tables/figures from the same saved outputs |

The specification is [`workflow_contract.json`](workflow_contract.json).
It fixes procedures and evidential limits, **not the result direction or its
explanation**. V2 results have already been seen: this is retrospective redesign,
not preregistration. Unexpected results and alternative explanations are welcome;
new analyses prompted by them must be labelled as post-inspection exploration.

The completed inventory/reconciliation/recount receipts retain the canonical
contract hash from [their execution checkpoint](https://github.com/zuizui0223/azami/blob/4e55670373c2e2fe0b8fc55b13ee065727369703/analysis/v3/workflow_contract.json).
The current contract updates the implementation labels and known source limits;
the historical receipts and their denominators are not silently rewritten.

## Integrate controls; do not move a long list of repairs upstream

- **Common design:** source identity, dates with hemisphere-aware encoding,
  photo/observation nesting, exposure uncertainty and question-specific spatial
  structure belong in preparation and the relevant model.
- **Measurement-specific evaluation:** crop shifts for orientation/outline;
  paired flower/context colour and photometric perturbations for colour;
  pixel size, sharpness and fixed-resolution checks for fine geometry.
  Save parameters and linked results before environmental interpretation.
- **Limited sensitivity comparisons:** assumptions not covered by the main
  design, such as dominant-taxon influence, range-scope differences and propagation
  of measurement uncertainty. Report coefficient ranges and support changes,
  not only preserved signs.

Native/introduced/unknown strata remain available. The old native-status table
covers only the v2 subset; it must not assign nativeness to the rest by default.
Define an ecological target after mapping coverage and the scientific question
are explicit. Use matching observation IDs for trait and environment summaries.

Univariate atlas associations do not establish independent effects of correlated
predictors. Conditional models need collinearity diagnostics on their actual
cohort and a saved covariate rationale. Tree placements are sensitivities, not
independent datasets; report direct-tip coverage and lambda.

## Two contributions, without required positive findings

**Methodological:** a traceable photo-to-feature procedure, coverage, failures,
technical error, comparability and reproducibility. Explain detector localization
separately from deterministic measurement. Report existing training/evaluation
metrics with their actual split provenance. The fully automated route does not
require new human reference measurements, but consistency alone does not establish
detector or physical-trait accuracy.

**Ecological:** image-feature associations and uncertainty for a stated population,
scale and model. V2-selected candidates remain post-selection; they need not
survive. Image variation is not automatically biological variation, low chroma
does not determine anthocyanin content, and spatial association does not establish
adaptation. Context colour and calendar timing do not by themselves rule out
illumination or developmental-stage explanations.

Keep all 27 original endpoints in the registry. **Correction to the earlier v3
description:** the five endpoints missing from the v2 atlas were not unmeasured
because of unfinished functions. Visibility and four colour fractions were already
computed at head level; the historical aggregator did not carry them into the
observation-level atlas input. They are not five QC rejections either.
Their [versioned recovery](../../reproducibility/v3_display_composition_recovery_20260907.json)
retains all 1,255,791 head records and produces values for 347,608 observations,
using 1,053,623 colour-QC-usable heads from 530,128 photos. It does not rewrite
the frozen v2 atlas or retroactively execute its missing analyses.

The new view averages eligible heads within each photo, then photo means within
each observation, giving each eligible photo equal weight. Original values and
missingness remain available. The four fractions sum to one and form a single
composition; joint hue has two columns but one inferential unit. Do not turn these
columns into independent biological traits or replace missing fractions with zero.
Changed functions, aggregation or endpoint extensions require explicit new versions.

## Run the full-source inventory

From the repository root (standard-library Python only):

```bash
python -m analysis.v3.workflow plan
python -m analysis.v3.workflow inventory \
  --metadata /path/to/photo_metadata_merged.csv \
  --legacy-v2-native /path/to/observation_native_status.csv \
  --out-dir /path/to/new-external-v3-inventory
python -m pytest -q tests/test_v3_workflow.py
```

The optional v2 native input only annotates overlap; it never filters the source.
The source member SHA-256, archived source identity and collection time are pinned
in the contract. Recover artifact `8066010557`, named
`ch1-inat-metadata-merged-full_inventory_20260703`, from run `28659167379` or its
owner archive; verify the extracted metadata against the contract. Actions
artifacts are temporary, not a durable public v3 release.

Use a new output directory; within the repository use ignored `local_data/` or
`outputs/`. Inputs and earlier runs are not overwritten. Local outputs are:

- `source_ledger.sqlite`: every photo record with its original row number and
  source-chunk link, plus observation-level counts; duplicates/conflicts are
  visible rather than silently removed.
- `observation_ledger.csv`: observation/photo aggregation without deleting the
  linked photo rows.
- `workflow_contract.json` and `source_inventory_report.json`: exact source
  and output hashes, canonical JSON identity, software versions, counts and
  stage status. Code text hashes explicitly normalize newlines to LF.
- On failure, `incomplete_run.json`: the partial run must not be used as complete.

The source CSV itself remains the immutable reference for all original fields,
including those not duplicated in the compact ledger. Keep it locally with the
ledger. Raw metadata, user identifiers, attribution, coordinates and images are
not added to GitHub. Only aggregate receipts are public.

The first inventory receipt describes the already merged input and remains an
unchanged historical checkpoint. The later reconciliation retains the original
1,122,855 chunk rows and 640,141 archived API records separately, together with
the unique-link view. It verifies every merged row against its exact original
source version, rather than treating photo-ID deduplication as lossless.

## Recovered processing history

The all-photo processing archive was independently recounted at photo, observation
and head level. All eight source/queue/screen/head/crop identity and state checks
had zero mismatches. These are historical processing counts, not new v3 results:

| Historical stage | Photos | Observations |
|---|---:|---:|
| Merged photo-metadata snapshot | 1,122,854 | 665,115 |
| Queued for the all-photo detector pass | 777,766 | 460,036 |
| At least one detected head | 637,745 | 406,582 |
| Strict spatially thinned observation view | Not a photo-level count | 46,276 |

The detector pass recorded 1,255,791 heads, 139,797 no-detection photos and 224
missing-image jobs. Another 345,088 source photos were not in this historical
queue. Do not convert unqueued or missing jobs into detection negatives, or
detection negatives into ecological absences. The archive retains the earlier
orientation/colour/outline measurements; it does **not** establish full-22-endpoint
coverage across 406,582 observations.

There are 137,492 detector-positive observations with more than one photo.
They are candidates for automated cross-photo matching, not already established
repeated measurements of the same head. Their extra photos must not silently
receive independent observation weight. Preserve the historical head-to-observation
summaries and create explicitly versioned alternatives when changing aggregation.

The recovered processing archive contains **no image files**. It preserves exact
processing records and old image/crop paths, not a verified current image cache.

## Executed local image workspace and versioned extraction

The [image workspace](../../reproducibility/v3_image_workspace_20260907.json)
retains the full 1,122,854-photo universe and all 1,122,901 recovered source links.
An audit of the available local development, audit and perturbation caches found
2,100 file declarations, representing 2,095 photo IDs and 2,092 distinct byte
objects. All declared files were available and decoded; 2,000 matched previously
archived file hashes and 100 received their first cache-inventory hashes.
Three byte objects have multiple photo IDs. Usage pools overlap and are not new
independent v3 evaluation splits. Shared observation, photo, byte or decoded-pixel
identities must be grouped before splitting or inference.

This local cache audit leaves **1,120,759 source photo IDs without a verified
image in this workspace**. It does not show that those images never existed,
are unavailable online, or are absent from every other storage location. No new
source-photo requests were made. The workspace receipt retains its execution
contract hash from [checkpoint a4f910a](https://github.com/zuizui0223/azami/blob/a4f910a7b40dabfa593efb4ff998b2e7925d02e5/analysis/v3/workflow_contract.json).

The [cached detector pass](../../reproducibility/v3_cached_detection_20260907.json)
processed all 2,092 distinct objects: 1,759 had detections and 333 did not.
It saved 2,853 head/context pairs, with no failed jobs, invalid crops or 300-box
limit flags. This is execution and coverage evidence, **not detector accuracy**.
The pinned weights were trained against automatic pseudo-labels; historical
development metrics are not independent precision or recall against true heads.
All photo/observation links remain in the source workspace rather than expanding
shared images into independent detections.

The measurement adapter uses corrected engine `56_run_primary_traits_continuous_v2.py`
(through its 55/52 compatibility dependencies) and extended engine 89. It does
not use the obsolete head-peduncle orientation definition from engine 52 alone.
Each detected head has all 27 endpoint slots, raw original/mirror values, finite
means even when QC fails, and endpoint-specific eligibility. The registry's pixel
requirements are explicit v3 gates. Primary foreground/floral and context masks,
crop identities, image dimensions and sharpness are saved. Engine and worker
errors remain distinct from ordinary low-quality or missing measurements.

The [executed baseline measurement](../../reproducibility/v3_cached_measurement_20260907.json)
completed all 2,853 head jobs without engine or worker errors and retained 77,031
endpoint rows (2,853 times 27), including null and ineligible values. In this
historical cached sample, 1,408 heads had at least one eligible endpoint but only
25 had all 27 eligible. Eligibility was 995 heads for orientation, 1,250 for colour
and visibility, 1,136 for outline, 160 for each architecture endpoint and 53 for
each surface endpoint. Thus requiring complete data for all endpoints would discard
most otherwise usable measurements in this cache. These are endpoint-coverage
counts, not accuracy estimates or estimates for the full source population.

Paired colour diagnostics use the same Lab/hue statistics on union floral pixels,
all non-head context and green non-head context. Context masks exclude the union
of all detected head boxes, not only the focal head. The new uniform floral
statistic is separate from legacy chroma, whose engine selects redmagenta pixels
when redmagenta is dominant and the floral union otherwise. Green context is not
verified leaf tissue, and undetected flowers can remain in it. These diagnostics
do not, by themselves, establish illumination correction or flower specificity.

### Run cached-image processing

Create a separate Python 3.12 environment; do not mix the image worker's OpenCV
package with `opencv-python-headless` from the numerical environment. The
[executed package inventory](../../reproducibility/v3_image_environment_20260907.json)
and per-run code/model/runtime hashes record the local environment. Cross-platform
bitwise equivalence has not been established.

```bash
python -m pip install torch==2.11.0+cpu torchvision==0.26.0+cpu --index-url https://download.pytorch.org/whl/cpu
python -m pip install -r reproducibility/v3-image-cpu-requirements.txt
python -m analysis.v3.build_image_workspace \
  --metadata /path/to/photo_metadata_merged.csv \
  --source-links /path/to/source-reconciliation/source_reconciliation.sqlite \
  --local-materials-root /path/to/azami_ch1_method_reproducibility_20260831 \
  --perturbation-root /path/to/perturbation_n100 \
  --out-dir /path/to/new-image-workspace
python -m analysis.v3.detect_cached_images \
  --workspace /path/to/image-workspace --weights /path/to/model/weights/best.pt \
  --out-dir /path/to/new-cached-detection
python -m analysis.v3.measure_cached_heads \
  --detection /path/to/cached-detection --out-dir /path/to/new-cached-measurement
python -m analysis.v3.verify_cached_pipeline \
  --workspace /path/to/image-workspace --detection /path/to/cached-detection \
  --measurement /path/to/cached-measurement --out-dir /path/to/new-verification
python -m analysis.v3.recover_display_composition \
  --heads /path/to/exhaustive_merged/exhaustive_continuous_head_level.csv \
  --out-dir /path/to/new-display-composition-view
```

The detector accepts only the model hash recorded in its receipt; the model
comes from archived artifact `8076736948`. The local cache adapter uses explicit
photo-ID/filename joins from archived manifests; it does not guess identities
from image similarity or inspect human labels. All inputs and outputs remain
local, including licenses/attribution and private observation links.

Detector and measurement workers accept `--limit N` for a bounded invocation.
Reusing the same output directory resumes pending jobs under an identical saved
execution context. Completed/error jobs are not silently retried. Changed code,
inputs or parameters require a new versioned output. A bounded pass does not
redefine the source denominator; missing, pending and no-detection are separate.
Cached baseline extraction and mirror QC are not the full crop, resolution and
photometric evaluation, and no v3 ecological model has been fitted here.

### Recover and audit original archives locally

The archive manifest is [`upstream_sources.json`](upstream_sources.json). Recovery
needs read access to the existing GitHub Actions archives (a `GH_TOKEN` in the
environment), but makes no new source-photo requests. It preserves exact ZIPs,
verifies selected extracted members, and will not overwrite changed local files.
Use ignored `local_data/` or an external directory. For example:

```bash
python -m analysis.v3.recover_upstream --out-dir /path/to/archives
python -m analysis.v3.reconcile_sources \
  --archives /path/to/archives --metadata /path/to/photo_metadata_merged.csv \
  --out-dir /path/to/new-source-reconciliation
python -m analysis.v3.audit_processing_history \
  --archive /path/to/archives/8269246732/source.zip \
  --metadata /path/to/photo_metadata_merged.csv \
  --out-dir /path/to/new-processing-recount
```

Source reconciliation reads original photo rows and raw API observation-photo
links, including the final chunk's compressed API file. The derived unique-link
table does not replace the original records: exact archives, source row/line
locators and hashes remain available locally. Multiple observations can reference
the same photo, so deduplicating image bytes must not discard their source links.

The first chunk has no archived raw API file; collection-time losses there cannot
be fully reconstructed. Collection used a photo-bearing-record query while the
API population could change. Matching chunk ID boundaries does not turn it into
a closed census of all iNaturalist observations. Its 25,000-observation collection
checkpoint and 24,998 metadata observation IDs differ by two; the missing raw API
file prevents resolving that difference from this archive.

### Preserved numerical checkpoints

The [offline numerical audit](OFFLINE_AUDIT.md) retains the existing PR #92
arithmetic and original subset verification. It does not define full v3 coverage,
its final ecological cohort or its model plan. Interim native-only preparation
files are retained locally, not published as a competing v3 entry point.

Frozen v2 reproduction remains documented in
[the public runbook](../../reproducibility/README.md).
Prepublication manuscripts and individual comment responses remain outside GitHub.

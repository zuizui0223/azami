# Chapter 1 v3: retain the source information, then define analysis views

V3 starts with the recovered acquisition snapshot, not the 46,276 observations
already selected for v2. The active design and readiness entry point is
`python -m analysis.v3.integrated_preflight`. The source-inventory utility remains
`python -m analysis.v3.workflow inventory`.

## Active five-stage workflow

The [integrated contract](integrated_workflow_contract.json) is the current
ordering and revision authority. Its [evidence index](integrated_evidence_index.json)
distinguishes pinned historical evidence, tested code and unexecuted work.
Earlier contracts below remain historical input definitions; the integrated
contract lists exactly which model/ordering provisions it supersedes.

| Stage | Question | Current implementation boundary |
|---|---|---|
| 1. Source | Which records support the native-range ecological question? | Exact 665,139-row source and 319,244-row native cohort recovered. Dates and known dependence added without membership changes; saved authority inputs reproduce the join offline. An 872-file source snapshot was restored on Actions and its returned worker packet verified locally. Historical HTTP identity is still unverified. |
| 2. Measurement | What can each photo measure, with what uncertainty? | A local 128-observation pilot retained all 27 raw slots for 315 detected heads. Its numerical integrity is verified, not physical accuracy. Full-cohort streaming and a separately qualified uniform-floral chroma definition remain unfinished; the 13 held routes stay held. |
| 3. Assessability | Which environments and taxa lose measurement support? | `assessability` reports endpoint/module attrition against the eligible native source. `model_design_diagnostics` separates raw, within, among and nuisance-adjusted exposure diagnostics. Both are tested, not yet run on the full realized native measurement cohort. |
| 4. Ecology | Do within- and among-taxon associations agree? | Joint slope algebra is verified. Covariance-aware pooling, matched-support scale contrasts and dependence-aware inference still require a validated runner. |
| 5. Synthesis | Do observed environmental and phenotype breadths covary? | Optional secondary Hypervolume B; common axes, bandwidth and sample-size/convergence qualification must precede it. Its failure cannot block otherwise supported primary ecology. |

This keeps two contributions together: an auditable image-to-distribution method,
and a test of module-specific ecological associations across scales. It does not
require any positive result or the survival of the v2 candidates. Technical
repeatability is not physical trait accuracy, and spatial association is not
adaptation. Full-production and ecological-fit authorization remain **false**.

The frozen 9,503-row environment diagnostic is not the final model cohort.
VIF/rank must be reported on each actual endpoint/module design, including
within-taxon and among-taxon support. Partial pooling addresses unstable taxon
slopes; it does not by itself repair selective photography or recording effort.

### All-27 measurement, without v2-result selection

The [native chunk contract](measurement_chunk_contract.json) extends the
completed pilot to bounded raw acquisition. The complete source is partitioned
into 7,023 remaining chunks, retaining whole known dependence components and
excluding the already measured 128 observations. The first
[two-chunk batch](native_measurement_batch_20260908.json) contains 100 observations
and 128 request candidates. Each chunk has at most 128 observations and 64
requests. Only these two chunks are authorized by this batch, not all 7,023.
The Actions worker keeps the pinned detector, measurement functions and CPU
runtime, saves all 27 slots and bbox perturbations, and uploads numerical
transactions to the unpublished draft. A completed transaction, including a
download failure, is restored rather than re-requested. Interrupted unverified
partial files are retained separately and never counted as complete.
No environmental values are read and no ecological model is fitted in this step.

The [executed two-chunk receipt](../../reproducibility/v3_native_raw_chunks_20260908.json)
now verifies 100 additional observations, 128 successful transfers, 232 heads,
6,264 all27 slots and 5,800 bbox-condition slots. Both protected bundles were
restored locally (1,412 files). One post-upload confirmation timed out; the failed
job restored all 64 completed photo units on its next attempt with **zero new
image requests**. The successful sibling job was not re-executed. The earlier
128-observation pilot also remains reused, including its terminal transfer failure.
The [next 16-chunk wave](native_measurement_wave_20260908_b.json) authorizes only
the next source-ordered 711 observations and 1,000 requests, with four concurrent
workers and the same algorithms. Pending wave counts are not completed measurements.
For Windows verification of Linux-produced units, use a separate LF decision-file
view only after its exact SHA-256 matches the execution report. Do not relax the
unit hash checks, rewrite the original checkout or treat a newline mismatch as a
scientific result.

GitHub Actions remains the intended production compute platform. Images are
temporary inputs, not Git objects or a required permanent image archive. The
preservation requirement concerns numerical results, source links and processing
provenance: protected numerical assets, exact local recovery and hash checks
before cleanup. It does not require the user to supply another external drive.
Because this repository is public, private identifiers/coordinates cannot simply
be added to an unrestricted workflow artifact. The protected transfer described
below has passed an [executed cloud/local roundtrip](../../reproducibility/v3_protected_numerical_replay_20260908.json)
in Actions run 34197391286: 872 files and 2,374,321,149 bytes restored on Actions,
followed by exact local verification of its returned 128-observation packet.
Finite workflow-artifact retention is not permanent archiving.

The [protected numerical replay contract](protected_numerical_replay_contract.json)
pins an existing **unpublished draft release**, source asset and all numerical
manifest identities. `protected_artifacts cloud-replay` restores the source on
Actions, derives the fixed 128-observation packet without fetching any images,
and uploads and downloads the resulting packet through that same protected
store. Anonymous release/asset access must return 404, and authenticated metadata
must confirm draft status before and after transfer. **Never publish this draft.**
Only an aggregate receipt goes into the ordinary workflow artifact. This route
requires no new external drive or permanent original-image archive; draft access
control is not encryption or a permanent publication archive. A successful
transport test does not authorize ecological fitting or full image processing.

The [subsequent output-preservation receipt](../../reproducibility/v3_executed_outputs_preservation_20260908.json)
adds 337 numerical/provenance files (225,714,529 restored bytes), including the
completed image pilot, raw environment checkpoints, rejected preliminary
selection, source-QC view and corrected candidate selection. The 70,163,053-byte
protected bundle was uploaded, downloaded afresh and restored with every file's
hash verified. Original images and private numerical rows are absent from this
repository and from the public workflow receipt.

The study's measurement scope remains all 27 registered endpoints. The number
14 describes the current original-image operating routes, not a replacement
scientific denominator. The other 13 architecture/surface endpoints are also
measured and retained, including finite QC-failed values and explicit nulls.
Their technical holds remain unchanged until separately specified image-only
evidence resolves them; availability in a new pilot does not itself admit them
to ecological inference. No all-27 complete-case requirement is imposed.

V2 significance, coefficient direction and headline survival are not selection
inputs. Reusing deterministic measurement functions is distinct from selecting
v3 questions by their v2 results. Because the redesign follows v2 and shares
source records, it is not a preregistration or an independent replication.

The bounded `stream_original_traits` command now has an offline output checker:
`python -m analysis.v3.verify_original_stream --input-dir PRIVATE_COMPLETED_PILOT
--out PRIVATE_NEW_VERIFICATION.json`. It reopens every numerical file, checks
all 27 slots for each detected head, reconciles source links and transfer states,
and recomputes the bbox-shift summaries and measurement eligibility. Reported
coverage distinguishes heads, photos and linked observations. It changes no
measurement threshold or ecological route and performs no environmental join.
A missing completion report, changed hash or incomplete denominator is an error,
not a partial pilot promoted to completion. Download/no-detection/QC failures
remain valid recorded outcomes, distinct from integrity failures.

The [completed native pilot](../../reproducibility/v3_native_original_pilot_20260908.json)
covered 128 observations and 196 photo links: 28 lacked an eligible licence;
167 of 168 requests succeeded, one failed, and 30 downloaded images had no
detected head. The 315 detected heads have 8,505 raw endpoint slots and 7,875
bbox-condition slots, all reconciled by the independent output checker. Paired
non-head context was available for 276 heads and green context for 228 heads.
No source image was persisted and no environment values were joined. These are
pilot throughput/support results, not full-cohort coverage or detector accuracy.

### Full-native environmental acquisition

`python -m analysis.v3.production_environment --source PRIVATE_ENRICHED_CSV
--out-dir PRIVATE_NEW_DIRECTORY` extracts all 15 environmental candidates for
the exact 319,244-observation source. It neither samples a smaller cohort nor
reads trait columns. Occupied raster blocks are read once per variable/month;
each numerical checkpoint has an exact hash and the remote object's identity.
`--resume` requires identical source, contracts, implementation and runtime, and
replays committed checkpoints offline. A failed transfer is incomplete acquisition,
not environmental missingness. Coordinates and numerical matrices stay private.

The [candidate process/VIF-10 rule](environment_production_contract.json) records
the hypotheses, literature, proxy limits and one-model proposal. Its four process
groups overlap physically; a conditional block coefficient is not an isolated
causal effect. All 15 columns are retained, including broader context. The
acquisition command deliberately does **not** execute variable selection or ecology.
The [executed acquisition/QC receipt](../../reproducibility/v3_full_native_environment_qc_selection_20260908.json)
records 103 completed checkpoints. A separate source audit found 1,252 unmasked
uint32-maximum GSP values and two BIO12 storage-ceiling values. The original
matrix and its preliminary VIF result are retained; that preliminary selection
is **not admitted**. The [explicit source-QC rule](environment_source_qc_contract_20260908.json)
creates a separate view, retains the original values and reason flags, and masks
only the affected working values. No source observations are deleted or imputed.

Run `python -m analysis.v3.environment_source_qc --matrix PRIVATE_RAW_MATRIX
--out-dir PRIVATE_NEW_QC_DIRECTORY`, then `python -m
analysis.v3.select_production_environment --matrix PRIVATE_QC_MATRIX
--expected-matrix-sha256 EXACT_SHA --out PRIVATE_NEW_SELECTION.json`.
The same weighted VIF-10 rule then retained nine variables in four process
blocks on 317,986 complete observations from 354 taxa; maximum VIF was 6.585
and condition number 6.891. CMI, Tmax and PET were removed from the candidate
representation, not deleted from the data. The integrated design now adopts this
source-QC selection as **one conditional model**, using the unchanged saved
predictor centers/SDs across modules and scales. Its four process blocks retain
4 wetting/moisture, 1 radiation, 3 heat/drying and 1 wind variables. The executable
definition in `environment_model.py` verifies the source receipt and reserves all
36 module-by-process-by-scale test slots; unavailable probabilities stay missing.
This supersedes the original drying/thermal alternatives, without rewriting their
historical receipts. Actual module-specific within/among and nuisance-adjusted
diagnostics, joint test calibration and the covariance-aware ecological runner
remain unfinished. No ecological model was fitted.

The [subsequent source-only identification audit](../../reproducibility/v3_source_environment_identification_20260908.json)
uses all 317,986 complete native source records, without traits. After shared
calendar/year and spherical-basis projection, pooled within/among matrices both
retain rank nine; their maximum VIFs are 5.530 and 9.792. Unadjusted among-taxon
BIO18 VIF is 11.105, so the earlier pooled source VIF is not a guarantee for every
scale. A conservative local nuisance-projection check retains rank nine for only
161 of 354 taxa. This is not the rank of the joint shared-nuisance model and does
not justify deleting the other 193 taxa. Sparse/aliased support must remain
explicit in the hierarchy. Imaging covariates and realized module membership are
absent from this source audit; actual model diagnostics remain necessary.

The following sections document the source inventory and earlier implementation
checkpoints. Their historical completion labels do not override the active
five-stage contract or the execution boundary above.

Current execution boundary (8 September 2026): the source-first design remains,
but full original-image processing and ecological fitting are **held**. The
[implementation correction](ecological_model_review_addendum_20260908.json)
records a synthetic counterexample to splitting common spatial residuals into
taxon-specific regressions. A joint block solver now matches the full interaction
design, including its cross-taxon HC3 covariance. This is verified numerical
machinery, not a completed ecological hierarchy or spatially robust inference.

The completed local pilot retains all 27 raw endpoints separately from the 14
candidate operational routes, individual bbox-shift values, paired background
colour, image hashes, detector geometry and link-level scheduling states. Its
real-image numerical output is verified; synthetic tests additionally exercise
failure cases, not physical accuracy. The worker accepts only an explicitly pinned reconciled
schedule or an explicitly labelled legacy metadata pilot; it rejects mixed inputs.
Ecological fitting still needs calendar/imaging nuisance design and dependence-aware
pooling. Source numerical transfer is verified; each new output bundle must be
preserved and restored through that route. The historical cloud pilot still stops
before retrieval. Earlier aggregate receipts are not retroactively promoted.

The source-cohort cloud workflow is also held and no longer runs automatically
on pushes: its historical recipe published only a report and then removed the
exact private source, authority join and cohort files. It must verify private
preservation and restoration before another cloud source freeze. The local
builders remain available. The subsequent local recovery below restores the exact
cohort and verifies local replay; it does not satisfy the separate off-device
private-archive gate or authorize cloud cleanup.

### Exact native-source recovery and local replay

The [executed recovery receipt](../../reproducibility/v3_source_cohort_recovery_20260908.json)
records the following checks, made without reading traits or fitting ecology:

- The six original archives reproduce all 665,139 source observations. A new LF
  serialization view matches the historical source CSV hash `5632f532...`; the
  Windows-native output is also retained, not overwritten.
- Reacquired WCVP responses and the pinned TDWG geometry reproduce the exact
  319,244-observation, 355-taxon cohort hash `b52503cd...`. The previous cohort
  report also agrees in full. These are source-support counts, not measured/QC
  endpoint counts.
- All 855 authority responses are saved in a private, hash-bound cache. Offline
  replay reproduces the native join, name-resolution and distribution tables
  byte for byte. The old HTTP responses/full join had no original pins: the new
  reconstruction is verified, not retroactively described as historical bytes.
- Enrichment preserves every observation and all 12 inherited fields, attaches
  exact source dates and full-source known dependence, and retains every photo
  link. There are 319,230 known components among the native observations; this
  does not prove all other observations independent.

All native-cohort records are northern-hemisphere records, so southern interaction
columns are structurally constant in this view. Source dates range from 1899 to
2026, including 29 records before 1980. They remain in the source ledger: syntactic
date validity does not verify photography dates. Calendar/year and exposure-support
rules must be fixed before fitting, without silent repair or outcome-led exclusion.
Date terms do not replace observed developmental-stage labels.

The first-page name-search limitation was also checked. Only `Cirsium` (genus)
and `Altissima` (complex) indicated further results; both ranks are outside the
unchanged ecological scope. No eligible source-rank query was truncated in this
reacquired response set. This does not independently verify image identification.

`recover_native_source_authority` injects cached acquisition into an isolated copy
of the unchanged classification helpers. It refuses changed source/authority
identities, altered cached responses, missing offline requests and existing output
directories. Outputs inside this repository must be under ignored `local_data/`.

### Reconciled native-photo schedule and private numerical restoration

The [executed schedule/replay receipt](../../reproducibility/v3_reconciled_schedule_private_replay_20260908.json)
records 319,244 native observations, 548,139 native observation-photo links and
548,123 unique photos. The ledger also keeps six external observation-photo links
sharing these photos and 1,080,801 archived photo versions. All six source archives
were checked against their pins and the reconciliation's original row/line hashes.
The schedule does not substitute the merged photo table or select a preferred
version when source rights or URLs disagree.

Of these photos, 443,811 meet the recorded-license/URL scheduling rules; 104,304
have unavailable or unsupported license codes and eight have missing/invalid URL
support. All states and links remain present. These are source-admission counts,
not successful downloads, detected heads, measurement QC or biological absence.
Only 256,162 of the 319,244 native observations currently have a request candidate;
this restriction must remain visible in the later source-to-measurement coverage
audit. No original images were fetched by this schedule/replay execution.

The 319,230 known native dependence components stay intact across 64 operational
partitions. Each unique photo is scheduled once. These components capture known
shared-photo dependence, not all duplicate individuals or spatial dependence, and
are not independent validation folds. The worker input omits taxon, coordinate,
date and source-user fields while retaining photo attribution and source locators.

The new input adapter selects a hash-ordered prefix of complete components for a
bounded pilot (at most 128 observations), including blocked photos. It neither
splits a component to fill the budget nor substitutes available photos for blocked
ones. The actual 128-observation input contains 196 photo links, of which 168 are
request candidates; the original and restored inputs have identical content hashes.
Offline synthetic tests also execute the worker and verify all endpoint slots and
blocked links. That is not an executed real-image pilot.

`reconciled_photo_schedule` requires `--enriched`, `--reconciliation`, `--archives`,
`--out`, and both `--expected-enriched-sha256` and
`--expected-reconciliation-sha256`. Inputs come from the pinned source chain;
output must be a fresh private directory. To use that schedule in a bounded worker
pilot, supply `stream_original_traits --reconciled-schedule <private-sqlite>` and
`--expected-schedule-sha256 <receipt-hash>`, plus the existing pinned detector,
measurement decision and private output arguments. Do not also supply legacy
`--metadata`, `--cohort` or observation-hash shard arguments. Full production
(`--pilot-observations 0`) remains rejected before network access.

`private_replay snapshot --selection <private-selection.json> --out <fresh-private-directory>`
copies an explicit list of files, each with a relative output `name`, absolute
source `path` and exact `sha256`, under `schema_version: 1` and `files`. It never
sweeps a directory, overwrites a prior snapshot or deletes source files. Restore
with `private_replay restore --snapshot-dir <snapshot> --out <fresh-private-directory>`
and `--expected-manifest-sha256 <snapshot-receipt-hash>`. Prefix both commands with
`python -m analysis.v3.`. Unsupported file types, unsafe paths, conflicting names,
changed hashes and incomplete snapshots fail closed.

The real local test copied and restored 872 files (2,374,321,149 bytes), including
the schedule, source-link database, source/native/enriched tables and 855 cached
authority responses. Every file was byte-verified. Raw acquisition ZIPs and the
large annotation database remain separately retained, not members of this bundle.
All copies were on one device: the tool provides neither encryption nor cloud
transport, and this test does not establish an off-device backup. Private raw
records, coordinates and identifiers must not be uploaded as public artifacts.
The cloud STOP and both production/ecological execution holds remain unchanged.
For Windows source rebuilds, `recover_lf_source(source, destination, expected_sha256)`
creates a separate LF byte view only if it exactly matches the historical pin.

```bash
python -m analysis.v3.recover_native_source_authority \
  --source-csv /private/source-historical-lf.csv \
  --expected-source-sha256 5632f532a63c8babdc20b023bd0d3c47424d69461df5de9793028e93e819f959 \
  --cache-dir /private/authority-cache --out-dir /private/new-native-join

# Replay the same requests with no network; use the saved manifest's exact hash.
python -m analysis.v3.recover_native_source_authority \
  --source-csv /private/source-historical-lf.csv \
  --expected-source-sha256 5632f532a63c8babdc20b023bd0d3c47424d69461df5de9793028e93e819f959 \
  --cache-dir /private/authority-cache --out-dir /private/new-offline-replay \
  --offline --expected-cache-manifest-sha256 SAVED_MANIFEST_SHA256
```

Next connect the pinned enriched cohort and reconciled links to actual image
scheduling, and verify private off-device restoration. Local retention is not a
permanent archive. Full image production and ecological fitting remain held.

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
| Recover | Six original metadata archives, the unthinned processing archive, verified local image bytes and full-source known-dependence groups | Full-source image availability and a rights-aware acquisition/retention schedule; processing records are not image files |
| Measure | Historical five-field recovery; cached detection, all-27 baseline extraction, completed 14-condition perturbations and component-weighted summary; independent scalar arithmetic verification | Full-source image execution and propagation of measurement uncertainty |
| Analyse | Executed full-source observation annotations, kept separate from image measurements; no v3 ecological fitting | Save exact estimands, cohort memberships, formulas, uncertainty, spatial/nesting design and multiplicity before fitting |
| Report | Aggregate receipts, full technical-sensitivity table and reproducible coverage/loss figure; full-source-linked observation measurement view | Model-specific ecological tables/figures from saved outputs |

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

The [scientific design rationale](design_rationale.md) connects the measurement
contribution with the ecological questions and their interpretation limits.
Original-image cloud streaming currently stops before source retrieval: no
authorized durable private numerical destination and verification path is
implemented. Discarding streamed image bytes must not discard numerical results,
source links or processing provenance. Aggregate-only receipts are insufficient.

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

### Full-source observation preparation

`prepare_observation_annotations.py` reads every original metadata row and every
available archived API observation. It retains their exact source locators and
separate metadata versions, then creates one annotation row per reconciled
observation. It does not load image measurements, ecological predictors or model
results. Raw coordinates and source identities remain in the local output only.

Where archived API fields exist, they take precedence over the collector's
flattened metadata. This preserves unknown boolean states that the old collector
could turn into false. API-missing fields are not silently filled from those
defaults. Where raw API is absent, flattened metadata remains an explicit fallback
with its limits. If preferred-source versions disagree, only the affected fields
become unavailable; every version remains saved. No latest-row or most-complete-row
selection is used.

Calendar preparation accepts exact observation dates, preserves their year and
day of year, and uses the actual 365/366-day year in sine/cosine terms. Public
latitude supplies a southern-hemisphere indicator and its interactions with both
calendar terms, allowing the fitted calendar relationship to differ across
hemispheres. Equatorial latitude has its own label. Missing or restricted location
does not receive an inferred hemisphere. These terms are not flowering-stage
measurements or proof that phenological confounding has been removed.

The observation row separately records public-location availability, unknown or
restricted privacy, positional accuracy, captive status and source taxonomy.
Having a public location with reported accuracy is not yet acceptance at any
environmental-grid resolution. Native status is explicitly unassessed in this
preparation; a later versioned range-status join must document its own coverage.
These are reversible annotations, not an ecological inclusion filter.

The [executed preparation receipt](../../reproducibility/v3_observation_annotations_20260907.json)
retains all 665,139 observations, 1,122,855 original photo-metadata rows and 640,141
archived API records. All source locators and record hashes matched the reconciled
ledger. API fields support 640,141 observation rows; 24,998 use the explicit
first-chunk metadata fallback. No preferred-source field conflicts were found.

Exact observation dates are available for 663,255 observations; 1,884 remain
date-missing. Public locations are present for 632,927 observations, while 26,468
are restricted, 5,440 lack valid coordinates and 304 retain a non-usable source
flag. Positional accuracy is absent for 155,718 observations, positive for 509,343
and reported as zero for 78. The hemisphere annotation identifies 611,733 northern,
21,190 southern and four exactly equatorial observations; 32,212 remain unknown.
None of these records was deleted. These counts describe source support, not
environmental alignment, phenology correction or a final ecological sample.

```bash
python -m analysis.v3.prepare_observation_annotations \
  --archives /path/to/upstream-recovery --reconciliation /path/to/source-reconciliation \
  --dependence /path/to/dependence-groups --out-dir /path/to/new-observation-annotations
```

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

The [executed dependence ledger](../../reproducibility/v3_dependence_groups_20260907.json)
retains every observation, photo and link. Shared photo IDs, exact bytes and exact
EXIF-oriented pixels connect observations into 665,103 known-dependence components;
33 components contain multiple observation IDs. The 2,092 cached image objects
belong to 2,082 components. This is known-identity grouping, not proof that every
duplicated scene has been found or that a component shows one biological individual.

The historical 270-image training manifest (211 training and 59 validation images)
was joined through its exact source filenames. No component crossed those recorded
training/validation subsets under the available identity evidence. The broader
historical development pool was conservatively propagated through components:
1,006 cached objects have recorded development exposure. Sixteen components cross
historical usage pools. Five deterministic operational folds keep each component
together, but do not create a fresh independent test set or repair past model use.

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

The [saved-product verification](../../reproducibility/v3_cached_pipeline_verification_20260907.json)
reopened every image object, crop and mask; checked exact hashes, dimensions,
binary masks and raw-to-endpoint table agreement; and confirmed all per-head
endpoint sets. This checks file/relationship integrity, not biological correctness.
Raw measurements, failed eligibility and diagnostic availability remain inspectable
without changing the source or rerunning a significance-selected subset.

Paired colour diagnostics use the same Lab/hue statistics on union floral pixels,
all non-head context and green non-head context. Context masks exclude the union
of all detected head boxes, not only the focal head. The new uniform floral
statistic is separate from legacy chroma, whose engine selects redmagenta pixels
when redmagenta is dominant and the floral union otherwise. Green context is not
verified leaf tissue, and undetected flowers can remain in it. These diagnostics
do not, by themselves, establish illumination correction or flower specificity.

### Technical perturbation evaluation

The saved [14-condition specification](perturbation_contract.json) covers identity
replay, four crop shifts, two resolution reductions, two gamma changes, two intensity
changes, two colour-balance changes and blur. The execution grid contains every
detected head and all 27 endpoint slots in every condition; it is not restricted
to heads with favourable baseline QC. Baseline endpoint values and statuses must
exactly replay the saved measurement before the other conditions proceed. A replay
mismatch stops new scheduling. Failed and unavailable values retain their slots.

Whole-image transforms precede crop extraction. The saved recipes, source hashes,
transformed-pixel hashes, masks' pixel hashes and linked original/mirror/QC results
make each comparison traceable without saving another copy of every transformed
image. Encoded-sRGB changes are specified digital stress tests, not physical camera
calibration or a validated illumination correction. The runner supports bounded
execution and resume of pending jobs, with an operating-system lock and unchanged
input/code/environment identities. A partial receipt cannot become a completed
summary. The [completed cached pass](../../reproducibility/v3_cached_perturbation_20260907.json)
contains 39,942 head-condition rows and 1,078,434 endpoint rows across all 2,853
heads, without baseline-replay, engine or worker errors. This completes the saved
digital probes on the historical local cache, not the full-source image pass.

The summary reports absolute and signed changes, rank agreement and all four QC
transitions. It averages heads within exact image, then images within known
dependence component, and gives components equal weight. Measurement-loss rates
include components with no surviving usable pair. Changes among surviving pairs
must always be read beside that loss. Separate development-exposure strata describe
recorded provenance, not independent validation. Joint hue uses circular angular
distance; the four colour fractions additionally use a joint composition distance.
These diagnostics do not add independent ecological endpoints. No favourable rank
threshold, accuracy claim or ecological conclusion is built into the summary.

The [executed summary](../../reproducibility/v3_technical_sensitivity_20260907.json)
contains all 1,218 rows: 14 conditions, 27 scalar endpoints plus two joint metrics,
and three recorded-exposure strata. The full [aggregate table](../../analysis_outputs/v3/technical_sensitivity_summary_20260907.csv)
is public. An [independent SQLite recomputation](../../reproducibility/v3_technical_summary_verification_20260907.json)
matched the counts, means and loss fractions of all 1,134 scalar rows. This checks
arithmetic, not accuracy; rank correlations, quantiles and joint metrics were not
independently recomputed in that check.

In the all-cached stratum, halving image dimensions removed baseline eligibility
for 55.2% of orientation, 49.2% of colour/display, 56.1% of outline, 80.9% of
architecture and 100% of surface measurements under component weighting. The
surface endpoints had only 53 eligible heads at baseline. This loss includes
crossing the saved minimum-pixel requirements; it is not an engine failure, a
biological absence or proof of inaccurate values. Finite QC-failed values remain
saved. Remaining paired orientation measurements had a rank correlation of 0.984
under halving, so high surviving-pair rank agreement must not hide substantial
eligibility loss.

Absolute shifts also matter: specified warm-channel gains changed paired chroma
by 6.54 Lab-chroma units on average, and blur with sigma 1 changed paired
orientation by 9.85 degrees, with the same image/component weighting. These are
native-unit digital-probe changes, not calibrated real-camera error estimates.
They cannot be compared directly with standardized ecological coefficients or
used alone to invalidate a v2 association. No favourable threshold was fitted to
these outcomes.

`plot_technical_coverage.py` draws all 27 baseline counts beside their eligibility
loss under every non-baseline condition. It requires the exact verified aggregate
and exports a local PNG and vector PDF. Its [figure specification](technical_coverage_figure_contract.json)
separates denominators and missing values explicitly. This is a diagnostic figure,
not a claim of journal-format compliance or a replacement manuscript.

### Observation measurement views without source deletion

The [executed aggregation](../../reproducibility/v3_observation_measurements_20260907.json)
retains all 665,139 observations and references all 1,122,901 source photo links.
Its logical inventory has 17,958,753 observation-endpoint slots, including missing
measurements. It reads image measurements and source identities, but no coordinates,
taxon assignments, environmental predictors or perturbation outcomes.

Eligible heads are averaged within each image and distinct eligible images receive
equal weight within each observation. Exact decoded-pixel aliases count once;
multiple encodings must agree in processing summaries. Unresolved distinct-pixel
versions of one photo ID remain explicit and are not selected by favourable QC.
Shared-image observations retain their dependence links. Hue components and the
four composition parts use joint eligible sets. Quality covariates and paired
flower/context contrasts use the same supporting heads and images as their values.
Raw finite means including QC failures are retained separately, never substituted
for eligible means.

The current cache supports 2,085 observations with selected images and 1,099 with
at least one eligible endpoint: 795 for orientation, 977 for colour/display, 898
for outline, 148 for architecture and 50 for surface endpoints. Paired green-context
colour is available for 745 observations. These are not representative full-source
coverage estimates, a final ecological cohort or a completed negative-control
regression. No v3 ecological model has yet been fitted.

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
python -m analysis.v3.build_dependence_groups \
  --workspace /path/to/image-workspace \
  --training-manifest /path/to/dataset/bootstrap_dataset_manifest.csv \
  --out-dir /path/to/new-dependence-groups
python -m analysis.v3.perturb_cached_heads \
  --workspace /path/to/image-workspace --detection /path/to/cached-detection \
  --measurement /path/to/cached-measurement --dependence /path/to/dependence-groups \
  --out-dir /path/to/new-perturbations --workers 2
python -m analysis.v3.summarize_perturbations \
  --perturbation /path/to/completed-perturbations --measurement /path/to/cached-measurement \
  --dependence /path/to/dependence-groups --out-dir /path/to/new-technical-summary
python -m analysis.v3.verify_perturbation_summary \
  --perturbation /path/to/completed-perturbations --measurement /path/to/cached-measurement \
  --dependence /path/to/dependence-groups --summary /path/to/technical-summary \
  --out-dir /path/to/new-summary-verification
python -m analysis.v3.build_observation_measurements \
  --workspace /path/to/image-workspace --detection /path/to/cached-detection \
  --measurement /path/to/cached-measurement --dependence /path/to/dependence-groups \
  --out-dir /path/to/new-observation-measurements
python -m analysis.v3.plot_technical_coverage \
  --summary analysis_outputs/v3/technical_sensitivity_summary_20260907.csv \
  --verification reproducibility/v3_technical_summary_verification_20260907.json \
  --out-dir outputs/new-technical-coverage-figure
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

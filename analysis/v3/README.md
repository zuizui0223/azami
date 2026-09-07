# Chapter 1 v3: retain the source information, then define analysis views

V3 starts with the recovered acquisition snapshot, not the 46,276 observations
already selected for v2. The current entry point is
`python -m analysis.v3.workflow inventory`.

The archived July 2026 photo metadata contain **665,115 observations and
1,122,854 photos**. These are records of photos, not proof that every image file
has been downloaded. The snapshot is not a count of all iNaturalist records today.

Executed evidence: [source inventory receipt](../../reproducibility/v3_source_inventory_20260907.json).
The full snapshot was inventoried twice with byte-identical copies of all four
outputs; no source row was removed. The 46,276 legacy observations all matched.

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
| Inventory | Exact source hash; streaming ledger preserving every photo row; observation aggregation; duplicate/conflict accounting; optional v2 overlap | Reconcile six pre-merge chunks and original API collection coverage |
| Recover | Source metadata are locally recovered | Reconcile all cached/downloaded photos, exact image versions and rights; raw observations without photos |
| Measure | Existing v2 code and results preserved as references | Full-inventory detector and measurement execution, linked failures, endpoint-specific automated technical evaluation |
| Analyse | Workflow specification; no v3 ecological fitting | Save exact estimands, cohort memberships, formulas, uncertainty, spatial/nesting design and multiplicity before fitting |
| Report | Aggregate source-inventory receipt | Measurement/association tables and figures from the same saved outputs |

The specification is [`workflow_contract.json`](workflow_contract.json).
It fixes procedures and evidential limits, **not the result direction or its
explanation**. V2 results have already been seen: this is retrospective redesign,
not preregistration. Unexpected results and alternative explanations are welcome;
new analyses prompted by them must be labelled as post-inspection exploration.

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

Keep all 27 original endpoints in the registry: 22 previously measured and five
with unfinished functions, not five QC rejections. Future functions or endpoint
extensions need new versions. Joint hue has two columns but one inferential unit.

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

The input was already merged upstream: the saved provenance reports **one
duplicate photo row removed from six chunks** (1,122,855 → 1,122,854).
V3's zero-row-loss inventory does not undo or certify that earlier merge.
Original API/query completeness and records without photos remain explicit
acquisition checks.

## Preserved checkpoints, not the v3 universe

The [offline numerical audit](OFFLINE_AUDIT.md) retains the existing PR #92
arithmetic and original subset verification. It does not define full v3 coverage,
its final ecological cohort or its model plan. Interim native-only preparation
files are retained locally, not published as a competing v3 entry point.

Frozen v2 reproduction remains documented in
[the public runbook](../../reproducibility/README.md).
Prepublication manuscripts and individual comment responses remain outside GitHub.

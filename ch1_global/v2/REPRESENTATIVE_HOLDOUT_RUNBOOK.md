# Representative trait-holdout human audit — runbook

This is the manual, human-labelling step that certifies AI trait measurement
validity (priority 1 in `../SUBMISSION_READINESS.md`). It cannot be automated:
people must look at photos and score traits blind. The chain is
`47 → 31 → 32 → 48`. Everything the humans see is blinded (no AI guess, taxon,
or locality), so the resulting agreement is an honest population estimate.

## Who and how much

- **Annotators:** ≥ 2 people who can score capitulum traits from photos (you +
  one or more lab members/collaborators). Two are required so a double-labelled
  subset yields inter-annotator agreement; more is better.
- **How many tasks:** default `47` draws 120 tasks per trait with 25 %
  double-labelled. For ~6 traits that is ~720 single-view judgements per
  annotator — a few focused sessions. Acceptance in `48` needs ≥ 50 assessable
  tasks and ≥ 10 species with evidence per trait, so 120/trait leaves headroom
  for `unassessable` photos.
- **Time:** budget a few seconds per task; ~1–2 hours per annotator per full
  packet is typical.

## Step 1 — build the blinded, representative, app-ready packet

Draw a representative sample stratified by taxon × 10° spatial block and copy the
images so the audit app can open it directly. `--candidate-units` is a table with
one analysable row per `annotation_unit_id`/`trait_id`
(`annotation_unit_id, trait_id, taxon_name, spatial_block_10deg, ai_candidate_state`);
`--packet-manifest` maps each unit to its `source_image, crop_path,
context_crop_path` (the same head/context crops used elsewhere in v2).

```bash
python ch1_global/v2/47_build_representative_trait_holdout.py \
  --candidate-units   analysis/candidate_units.csv \
  --packet-manifest   analysis/packet_manifest.csv \
  --packet-root       analysis/packet_images \
  --out-dir           analysis/representative_holdout \
  --batch-name        holdout_v1 \
  --n-tasks-per-trait 120 --double-label-fraction 0.25
```

Outputs: `representative_trait_holdout_tasks/` (blinded packet + copied images,
`app_ready`) and `representative_trait_holdout_key/` (the private key — **do not
show annotators**).

## Step 2 — each annotator scores the packet (blinded)

Run the existing Streamlit app once per annotator, each writing to their own
out-dir:

```bash
streamlit run ch1_global/v2/31_trait_audit_app.py -- \
  --tasks   analysis/representative_holdout/representative_trait_holdout_tasks/representative_trait_holdout_tasks.csv \
  --out-dir analysis/responses/rachel \
  --annotator rachel
```

The app shows head crop + peduncle context + source image, hides everything
model/taxon/locality, and lets the annotator pick states or `unassessable`. It
writes `trait_audit_responses.csv`.

## Step 3 — compile responses into canonical annotations

```bash
python ch1_global/v2/32_compile_trait_audit_responses.py \
  --tasks    analysis/representative_holdout/representative_trait_holdout_tasks/representative_trait_holdout_tasks.csv \
  --ontology ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --response rachel=analysis/responses/rachel/trait_audit_responses.csv \
  --response annotator2=analysis/responses/annotator2/trait_audit_responses.csv \
  --out-dir  analysis/canonical
```

## Step 4 — certify and gate (the acceptance decision)

```bash
python ch1_global/v2/48_evaluate_representative_trait_holdout.py \
  --canonical-annotations analysis/canonical/<canonical_annotations>.csv \
  --holdout-key           analysis/representative_holdout/representative_trait_holdout_key/representative_trait_holdout_key.csv \
  --analysis-units        analysis/ai_trait_observation_level_long.csv \
  --ontology              ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --out-dir               analysis/certification \
  --measurement-mode      conservative
```

`representative_holdout_trait_gate.csv` gives each trait a
`accept_for_analysis` / `manual_review_or_exclude` decision. Default thresholds:
population-accuracy Wilson-lower ≥ 0.80; ≥ 10 species with evidence;
worst-species accuracy ≥ 0.60; human–human agreement Wilson-lower ≥ 0.80; ≥ 50
assessable tasks. **Only `accept_for_analysis` traits may headline the
trait–environment atlas and the integration analysis (`49`).** Failing traits
go back to prompt/ontology/detector work or are excluded.

## Also run the detector audit (same idea, boxes not traits)

`17_evaluate_detector_audit.py` needs a small set of human-adjudicated
visible-capitulum boxes to report detector precision/recall (audit finding C1);
build its queue with `16` and adjudicate the same way.

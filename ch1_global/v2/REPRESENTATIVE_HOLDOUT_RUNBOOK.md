# Representative trait-holdout human audit — runbook

This is the independent human-validation step that decides which automated
image traits may headline Chapter 1. It does **not** replace the executed
216-species association atlas. It certifies the measurements that feed that
atlas.

The production chain is:

```text
51 build public/private holdout packets
→ 31 blinded annotation app
→ 32 compile human responses
→ 48 certify each trait
→ retain accepted traits in the existing atlas and integration analysis
```

## Scope

The production holdout currently targets the three whole-capitulum traits:

- `capitulum_orientation`
- `corolla_colour_class`
- `capitulum_shape`

Fine involucral traits remain supplementary until separately validated with
suitable close-up images.

## Why the packet contains two sampling strata

Two validation questions require different samples and must not be mixed.

1. **Population accuracy** — estimated from 120 tasks per trait. Tasks are
   allocated to taxon × 10-degree spatial cells in proportion to the number of
   eligible analysis records in each cell, then selected by a stable hash within
   cells. These rows have `selection_stratum=representative_population` in the
   private key.
2. **Species transferability** — evaluated with a separately marked,
   model-margin-free top-up. Twelve deterministic taxa per trait are brought to
   seven records each when necessary. These extra rows have
   `selection_stratum=species_diagnostic_topup`.

Script `48` uses only the first stratum for population accuracy. It may use both
strata for the per-species minimum and worst-species diagnostic. Therefore the
species top-up cannot inflate the population-accuracy estimate.

## Step 1 — generate the packets from the executed global analysis

The GitHub Actions workflow
`.github/workflows/ch1-build-representative-primary-holdout.yml` uses the merged
216-species analysis artifact and runs:

```bash
python ch1_global/v2/51_build_representative_holdout_packet.py \
  --head-traits global_ai_analysis_units/ai_trait_head_level_with_metadata.csv \
  --yolo-crops merged_shards/global_ai_merged_yolo_crop_metadata.csv \
  --ontology ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --out-dir representative_holdout \
  --traits capitulum_orientation corolla_colour_class capitulum_shape \
  --measurement-mode conservative \
  --n-tasks-per-trait 120 \
  --double-label-fraction 0.25 \
  --species-diagnostic-taxa 12 \
  --species-diagnostic-records 7
```

The builder re-downloads the same iNaturalist **medium** image size used during
AI inference and recreates the head and context crops from the stored YOLO box.
If any selected image cannot be reconstructed, the run fails rather than
silently dropping it after selection.

The workflow publishes three separate artifacts:

- `ch1-representative-primary-public-*` — full blinded packet for the primary annotator;
- `ch1-representative-secondary-public-*` — 25% double-labelled subset;
- `ch1-representative-private-key-*` — taxon, AI state and sampling stratum.

Never give the private-key artifact to an annotator.

## Step 2 — annotate in Streamlit

Upload the primary public ZIP to the Cloud app and use one annotator ID, for
example `rachel`. A second person uploads only the secondary public ZIP under a
different annotator ID.

For local use:

```bash
streamlit run ch1_global/v2/31_trait_audit_app.py -- \
  --tasks representative_trait_holdout_primary/blinded_trait_audit_tasks.csv \
  --out-dir responses/rachel \
  --annotator rachel
```

The app shows the head crop, head-plus-peduncle context and source image. It
hides the AI state, species, locality and sampling stratum. Use `unassessable`
instead of guessing.

Download and retain `trait_audit_responses.csv` regularly because Cloud storage
is temporary.

## Step 3 — compile both response files

```bash
python ch1_global/v2/32_compile_trait_audit_responses.py \
  --tasks representative_trait_holdout_primary/blinded_trait_audit_tasks.csv \
  --ontology ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --response rachel=responses/rachel/trait_audit_responses.csv \
  --response annotator2=responses/annotator2/trait_audit_responses.csv \
  --out-dir canonical
```

The secondary response contains only the designated subset. The compiler joins
responses by `task_id`; it does not require both people to score every task.

## Step 4 — certify the head-level AI states

Use the **head-level** table, not the observation-level table:

```bash
python ch1_global/v2/48_evaluate_representative_trait_holdout.py \
  --canonical-annotations canonical/trait_audit_annotations_canonical.csv \
  --holdout-key representative_trait_holdout_key/representative_trait_holdout_key.csv \
  --analysis-units global_ai_analysis_units/ai_trait_head_level_with_metadata.csv \
  --ontology ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --out-dir certification \
  --measurement-mode conservative \
  --primary-annotator rachel
```

`48` reads `analysis_state_ai_conservative` from the head-level table. Legacy
test fixtures using `observation_ai_conservative_state` remain supported, but
the real certification is head-level because each audit task is one detected
head.

Default gates per trait are:

- at least 50 assessable representative-population tasks;
- population-accuracy Wilson lower 95% bound ≥ 0.80;
- at least 10 taxa with at least 5 scoreable records;
- worst eligible-taxon accuracy ≥ 0.60;
- human–human agreement Wilson lower 95% bound ≥ 0.80.

The key result is `representative_holdout_trait_gate.csv`:

- `accept_for_analysis` — the trait may headline the existing global atlas;
- `manual_review_or_exclude` — retain only as exploratory, revise the
  measurement system, or exclude from confirmatory claims.

## Step 5 — preserve the original analysis

Do not discard or rerun the entire 216-species atlas merely because one trait
fails. Re-headline the existing results using only traits that pass both:

1. the measurement-validity gate above; and
2. the existing spatial-block and held-out-species generalisation checks.

The all-ensemble versus conservative measurement modes remain an
errors-in-variables sensitivity analysis.

## Separate detector gate

Trait validity does not certify the YOLO boxes. Report detector precision and
recall separately with `16_build_detector_audit_manifest.py` and
`17_evaluate_detector_audit.py` on human-adjudicated visible-capitulum boxes.

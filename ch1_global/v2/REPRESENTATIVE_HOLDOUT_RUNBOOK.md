# Chapter 1 representative human check

The design is intentionally simple.

## What is checked

Three whole-capitulum traits:

- `capitulum_orientation`
- `corolla_colour_class`
- `capitulum_shape`

## Packets

- **Primary annotator:** 360 tasks = 120 tasks × 3 traits
- **Second annotator:** 90 tasks = 30 tasks × 3 traits
- **Private key:** used only after both people finish

There is no separate species-diagnostic packet and no extra annotation stratum.
The 120 tasks for each trait are sampled in proportion to the available records
within taxon × 10-degree spatial cells.

## Why this is enough

The human check answers only two questions:

1. Does the AI label agree with a human label?
2. Do two humans agree with each other?

Transfer to unseen species is already tested by the existing held-out-species
analysis in the 216-species atlas, so annotators do not need another species
specific sampling design.

## Annotation

Upload the primary ZIP to the Streamlit app with annotator ID `rachel`.
A second person uploads the smaller secondary ZIP with a different ID.
Use `unassessable` instead of guessing.

Never upload or inspect the private key while annotation is in progress.

## Compile responses

```bash
python ch1_global/v2/32_compile_trait_audit_responses.py \
  --tasks representative_trait_holdout_primary/blinded_trait_audit_tasks.csv \
  --ontology ch1_global/v2/ontology/ch1_trait_ontology.csv \
  --response rachel=responses/rachel/trait_audit_responses.csv \
  --response annotator2=responses/annotator2/trait_audit_responses.csv \
  --out-dir canonical
```

## Evaluate

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

A trait passes when:

- at least 50 tasks are assessable;
- the lower 95% Wilson bound for AI–human agreement is at least 0.80;
- the lower 95% Wilson bound for human–human agreement is at least 0.80.

Per-species accuracy is still exported as a descriptive table, but it is not an
extra production gate. Only traits with `accept_for_analysis` should headline
the existing global atlas.

# Chapter 1 — submission-readiness review and gate

Review date: 2026-07-10. Scope: can the current repository (v2 pipeline +
executed 216-species results) be submitted as Chapter 1? Verdict and the
prioritized path to "yes" are below. This review reads the executed Actions
artifacts' provenance (association-atlas run summary, workflow logs) and the
committed v2 code; it does not re-run the pipeline.

## What already exists (and is genuinely strong)

The v2 pipeline is built and has been executed at research scale. Confirmed
from the successful `Ch1 run trait environment association atlas` run
(2026-07-09) and the merged inference artifacts:

- **216-species global trait dataset** (`ch1-global-ai-trait-analysis-research_grade_216_species_v1`),
  produced by a multi-stage pipeline: multi-class capitulum detector →
  head/context crops → zero-shot CLIP + AI-ensemble trait scoring →
  observation-level analysis units with two measurement modes
  (`observation_ai_all_state`, `observation_ai_conservative_state`).
- **Trait–environment association atlas** over CHELSA climate: 56 trait-state
  outcomes screened, **42 eligible**, 336 coefficient rows, **630 grouped
  validation rows, all successful**.
- **Design rigor already in place** (matches most of `METHODS_AUDIT.md`):
  - Grouped cross-validation by **10-degree spatial block** *and* by
    **held-out species** (climate-only models); species-adjusted (taxon-dummy)
    models as a phylogeny-free taxonomic control.
  - Pre-declared **trait tiers** (primary/secondary/exploratory) and eligibility
    gates (≥100 pos/neg obs, ≥10 positive species, ≥10 spatial blocks, ≥5
    species with within-state variation).
  - **Coordinate-precision sensitivity** (full cohort vs ≤10 km positional
    accuracy) and **measurement-stability** (all-ensemble vs conservative).
  - Trait ontology with explicit assessability rules and `do_not_infer_from:
    species name`; mucilage excluded from automated scoring.
  - Non-causal framing stated in the plan and every report.

This is not a plan — it is an executed, rigorously structured macroecology
result. The earlier assessment that "there are no results" was wrong.

## The decisive blocker (one, and it is the project's own stated rule)

**The automated trait labels are not yet certified as valid measurements, and
the pipeline's own governance says they cannot be with what has been run.**

`33_evaluate_trait_audit_against_clip.py` — the trait-validation step — audits a
deliberately **margin-stratified, non-representative** sample (low-margin
uncertainty + high-margin calibration). Its own report states:

> "reports agreement only within the deliberately selected audit tasks … never
> interprets this non-random audit sample as global population accuracy"
> and: "A separate independent, representative holdout is required before any
> automatic acceptance policy."

Its best possible verdict is `candidate_assist_only_pending_independent_holdout`.
So by design it can never certify the labels as analysis-ready. Two facts follow:

1. **The required representative + species-held-out validation did not exist**
   in the pipeline (only the margin-stratified assist audit did).
2. **No human audit labels are committed at all** (`33` needs
   `--canonical-annotations`; none exist; the only `validation_*` files are the
   frozen v1 orientation artifacts). So even the assist audit has not been run
   with real human data.

Consequently the 216-species atlas — impressive as it is — sits on trait labels
whose **measurement validity is unproven**. The atlas's own grouped CV validates
that *given the labels* the climate associations generalize across space and
species; it does **not** validate that the labels are correct (a different
pillar). A reviewer will reject "we measured colour/shape/armature from
citizen-science photos with CLIP" until measurement validity — including that
the model reads the trait, not species identity — is demonstrated on a
representative, species-aware sample.

## Priority-ordered path to submittable

1. **Certify measurement validity (the linchpin, priority 1).**
   - Build a **representative, taxon × spatial-block stratified, blinded**
     holdout → *added this review:* `ch1_global/v2/47_build_representative_trait_holdout.py`.
   - Collect **human labels** on it via the existing audit app (31), with a
     double-labelled subset for inter-annotator agreement. *(Manual step — only
     the author can do this; it is now unblocked because the sampler + evaluator
     exist.)*
   - Evaluate + gate → *added this review:*
     `ch1_global/v2/48_evaluate_representative_trait_holdout.py`, which reports
     per-trait **population accuracy** (Wilson lower bound), a **per-species
     breakdown with a worst-species floor** (the leave-one-species-out
     diagnostic), and human–human agreement, then emits a per-trait
     `accept_for_analysis` / `manual_review_or_exclude` decision. Only accepted
     traits may headline the atlas.
   - Report **detector precision/recall** via the existing
     `17_evaluate_detector_audit.py` on adjudicated boxes.
   - Acceptance criteria (pre-declared, tunable in `48`): population accuracy
     Wilson-lower ≥ 0.80; ≥ 10 species with evidence; worst-species accuracy
     ≥ 0.60; human–human agreement Wilson-lower ≥ 0.80; ≥ 50 assessable tasks
     per trait.
2. **Re-headline the atlas on certified traits only.** Re-read the atlas
   artifact's per-trait CV AUCs; keep as headline results only traits that
   (a) passed step 1 and (b) generalize under spatial-block CV (the plan's own
   rule: "do not claim generalization if spatial-block validation is poor").
   Carry predicted-state uncertainty (all vs conservative modes) as the
   errors-in-variables sensitivity.
3. **Write the manuscript** (intro / methods / results / discussion + figures).
   No manuscript exists yet; results currently live only as CI artifacts. Frame
   novelty per `README.md` §1.4 (validated multi-trait image phenomics; global
   syndrome/integration structure; a phylogeny-free baseline for an intractable
   clade). Name the pollinator hypothesis as a candidate mechanism for Ch.4;
   analyse no pollinator variable.
4. **Durable data deposit.** The Actions artifacts have a 90-day retention
   (the association atlas expires 2026-10-07) and are not a citable archive.
   Commit the key result tables and deposit the dataset (e.g. Zenodo), and
   record the environmental-extraction provenance (CHELSA/topography/soil
   version, resolution, extraction date) that `METHODS_AUDIT.md` B3 flags.

## Secondary items (do before submission, not blockers to start)

- **Trait integration / constraint analysis** (`README.md` §3.5) — *implemented
  this review:* `ch1_global/v2/49_run_trait_integration_and_constraint.py`
  computes, at species level, pairwise trait association (bias-corrected
  Cramér's V), an integration-magnitude permutation test, and a constraint /
  forbidden-combination test, from the observation-level long table. It is a
  core novelty leg (§1.4-2). Run it on the certified traits once step 1 passes
  (it consumes AI states, so headline only the accepted traits).
- **Aggregation-level sensitivity** (species-mean and species×grid as headline)
  — confirm the atlas reports these, not only observation-level.
- **Multiple-testing discipline** across 42 eligible outcomes: report effect
  sizes + CIs and cross-tier consistency rather than p-value counts.

## Bottom line

The analytical core exists and is strong; the gap is narrow but decisive and is
the project's own declared gate: **certify the automated trait measurement on a
representative, species-aware human holdout, then write it up.** This review
adds the two missing scripts (`47`, `48`, with tests) that turn that gate from
"not implemented" into "run the human audit and evaluate." Everything else is
staged behind that certification.

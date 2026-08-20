# Azami → EAzami handoff and inference boundary

Status: 2026-08-20

## One rule

**Azami measures the global observational pattern. EAzami tests whether candidate biological mechanisms can reproduce that pattern together with independent interaction/fitness evidence.**

The downstream simulation does not make the Azami correlations causal.

## Azami owns the observation layer

The submission-facing Azami chapter owns results generated from the public-image phenomics pipeline, including:

- repeated continuous head-level trait distributions;
- hierarchical visible variance below assigned-taxon means;
- within-taxon environment–trait associations;
- grouped spatial-model results;
- among-taxon environmental sorting;
- image-derived involucre/spine-like proxy associations;
- sensitivity/robustness analyses attached to those observational results.

Canonical frozen result source:

`manuscript/final_claims.json`

The principal handoff patterns currently include:

- below-taxon visible variance across primary endpoints: **0.5886–0.9307**;
- cross-module visible-dispersion × environmental-association rho **−0.0569**, bootstrap interval **−0.265 to 0.155**;
- orientation × BIO1 beta **+0.017059**;
- chroma × BIO12 beta **+0.039324**;
- globally weak/model-dependent gross-outline support;
- involucre-proxy × BIO4 betas: roughness **+0.0975**, spread **+0.0937**, maximum spine-like projection **+0.0911**.

These remain **observational patterns**. They do not by themselves establish plasticity, adaptation, pollinator causation, herbivore defence, evolutionary rate or adaptive radiation.

## EAzami owns the mechanism layer

The separate `zuizui0223/EAzami` repository may import frozen Azami pattern summaries and combine them with independent evidence such as:

- quantitative pollinator responses;
- florivory / predispersal seed-predation responses;
- reproductive-fitness effects;
- orientation/protection experiments;
- population ancestry, plastid and cytotype data;
- focal field manipulations and repeated interaction measurements.

EAzami may then compare generative/mechanistic models asking whether combinations of environment, pollinator benefit, antagonist cost, trait-module independence and population-specific interaction regimes can reproduce the joint pattern vector.

Current downstream files include:

- `data/evidence/capitulum_pattern_reduction_targets_v1.csv`;
- `analysis/simulate_capitulum_pattern_reduction_v1.py`;
- `docs/CROSS_LAYER_CAPITULUM_PATTERN_REDUCTION_2026-08-20.md`;
- subsequent uncertainty-weighted pattern-reduction versions maintained only in EAzami.

## What may cross the repository boundary

### Azami → EAzami

Allowed:

- frozen point estimates and uncertainty intervals;
- endpoint/module identities;
- observation hierarchy and scale;
- cohort/sample sizes;
- explicit noncausal claim boundaries.

Not allowed:

- re-labelling image variance as genetic variance;
- treating environment association as selection coefficient;
- converting image-derived involucre/spine proxies into direct botanical defence traits;
- treating a species mean as a resolved evolutionary state when within-taxon variation is documented.

### EAzami → Azami

Allowed only as downstream context:

- a statement that independent evolutionary/mechanistic work tests explanations for the observed Azami patterns.

Not allowed in the Azami results or causal interpretation:

- simulation-derived claims that a mechanism is true;
- EAzami model ranking presented as evidence that the Azami climate association is adaptive;
- pollinator/herbivore mechanisms inferred back onto Azami observations without direct causal data;
- adaptive-radiation claims.

## Why this separation matters

The two repositories answer different questions:

| layer | question | valid output |
|---|---|---|
| Azami | What global capitulum patterns are visible, and how are they environmentally structured? | observational pattern and uncertainty |
| EAzami | What minimum biological mechanism family can reproduce those patterns plus independent interaction/fitness evidence? | mechanism adequacy / failure diagnosis |
| focal doctoral field work | Do those mechanisms operate in ancestry-resolved populations and affect fitness? | causal biological test |

A successful EAzami simulation therefore **does not upgrade Azami from correlation to causation**. Its value is to identify which additional biological measurements discriminate among plausible explanations.

## Current practical handoff

Azami Chapter 1 is not expanded with simulation, meta-analysis or field-mechanism code. Once its frozen observation patterns are exported, downstream development occurs in EAzami.

If EAzami finds that a proposed mechanism family cannot reproduce the Azami pattern vector, revise the downstream mechanism model or collect new biological data. Do **not** revise the frozen Azami empirical result merely to make a mechanism fit.

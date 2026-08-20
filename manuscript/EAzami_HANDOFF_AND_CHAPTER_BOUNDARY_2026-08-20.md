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

## EAzami owns the mechanism-reduction layer

The separate `zuizui0223/EAzami` repository imports only frozen Azami pattern summaries and combines them with independent interaction/fitness evidence.

Current canonical downstream files are:

- `data/evidence/macro_interaction_pattern_targets_v2.csv` — Azami observation targets plus literature/meta targets;
- `data/evidence/interaction_quantitative_pattern_ledger_v1.csv` — concrete pollination, antagonist, orientation/protection and fitness values/nulls;
- `analysis/simulate_macro_interaction_pattern_reduction_v2.py` — multi-seed ABC-like structural sufficiency screen;
- `data/evidence/macro_interaction_pattern_reduction_result_v2.json` — frozen reproducible simulation result;
- `docs/MACRO_INTERACTION_PATTERN_REDUCTION_RESULT_V2.md` — interpretation and failure diagnosis;
- `docs/INTERACTION_QUANTITATIVE_PATTERN_SYNTHESIS_V1.md` — biological-interaction pattern synthesis.

The current EAzami quantitative interaction bundle includes examples such as:

- *C. purpuratum* display → bumblebee visitation `R² = 0.637` and heads probed `R² = 0.533`;
- *C. purpuratum* display → seed predation `R² = 0.44` (Nikko) and `0.26` (Kawamata);
- direct *Cirsium* reduced-herbivory / ambient-herbivory seed-output meta RR **2.674** (95% CI **2.388–2.993**), equivalent to about **62.6%** loss of potential seed output under ambient herbivory in the harmonized set;
- *C. canescens* developing-seed insect damage **21–54%** across three sites;
- *C. pitcheri* infested capitula with **60% fewer mature seeds**;
- *C. arvense* scent compounds attracting both pollinators (>10 species) and florivores (16 species);
- *Cremanthodium campanulatum* natural nodding versus artificial erect achene set **56.3% vs 15.7%**, while pollinator orientation preference was not detected.

These values constrain downstream mechanisms. They are not results of Azami Chapter 1.

## Current pattern-reduction result

EAzami v2 compares five mechanism families across **720 prior draws per family**, retaining the top 5% core-fitting draws and asking whether they also reproduce held-out interaction patterns.

The robust ordering is:

1. full environment + pollinator + antagonist trade-off with modular evolvability;
2. full trade-off with common lability;
3. antagonist only;
4. pollinator only;
5. environment only.

The important result is the separation between **full multi-driver models and single-driver models**, not the small difference between the two full models.

Accepted-set median pattern distance / held-out reproduction are approximately:

- modular full: **0.287 / 0.706**;
- common-lability full: **0.297 / 0.689**;
- antagonist only: **0.419 / 0.417**;
- pollinator only: **0.429 / 0.317**;
- environment only: **0.586 / 0.133**.

Current inference:

> **Azami's environmental/trait patterns plus independent biological-interaction effects are structurally difficult to reproduce with environment, pollinators or antagonists alone. Joint ecological trade-offs are strongly required; modular evolvability remains plausible but is not yet decisively distinguished from common lability.**

This is a downstream structural-sufficiency result, not a causal interpretation of Azami.

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

- a statement that independent evolutionary/mechanistic work tests explanations for the observed Azami patterns;
- a pointer to the current EAzami target registry/result.

Not allowed in the Azami results or causal interpretation:

- simulation-derived claims that a mechanism is true;
- EAzami model ranking presented as evidence that the Azami climate association is adaptive;
- pollinator/herbivore mechanisms inferred back onto Azami observations without direct causal data;
- adaptive-radiation claims.

## Why this separation matters

| layer | question | valid output |
|---|---|---|
| Azami | What global capitulum patterns are visible, and how are they environmentally structured? | observational pattern and uncertainty |
| EAzami | What minimum biological mechanism family can reproduce those patterns plus independent interaction/fitness evidence? | mechanism adequacy / failure diagnosis |
| focal doctoral field work | Do those mechanisms operate in ancestry-resolved populations and affect fitness? | causal biological test |

A successful EAzami simulation therefore **does not upgrade Azami from correlation to causation**. Its value is to identify which additional biological measurements discriminate among plausible explanations.

## Current unresolved downstream mechanisms

EAzami v2 specifically fails to represent several held-out patterns cleanly:

- morph-specific colour choice;
- trait-specific defence nulls such as *C. discolor* stickiness;
- orientation-mediated rain/UV protection;
- early-morning versus all-day orientation effects;
- year-dependent tolerance;
- pollinator dependence / reproductive assurance.

The next model expansion should be targeted to these failure modes rather than increasing generic model complexity.

## Current practical handoff

Azami Chapter 1 is **not** expanded with simulation, meta-analysis or field-mechanism code. Once its frozen observation patterns are exported, downstream development occurs in EAzami.

If EAzami finds that a proposed mechanism family cannot reproduce the Azami pattern vector, revise the downstream mechanism model or collect new biological data. Do **not** revise the frozen Azami empirical result merely to make a mechanism fit.

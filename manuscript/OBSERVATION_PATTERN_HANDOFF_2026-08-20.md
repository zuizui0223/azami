# Observation-pattern handoff from Azami to EAzami

Status: 2026-08-20

## Boundary

Azami Chapter 1 remains the **observational global macro layer**. Its job is to estimate reproducible patterns in continuous capitulum phenotypes and their environmental/spatial structure.

It does **not** fit pollinator–antagonist mechanism models or use simulation to turn correlation into causal adaptation.

EAzami may consume selected, frozen Chapter 1 summaries as **patterns to be reproduced**, not as causal parameters.

## Canonical observational source

The machine-readable source is:

`manuscript/final_claims.json`

The initial EAzami reduction programme currently uses the following Azami pattern classes:

1. below-taxon visible variance range across primary endpoints;
2. near-zero common relationship between visible dispersion and environmental-association strength;
3. module-specific environmental structure, especially orientation and visible colour;
4. weak/model-dependent global gross-outline environment structure;
5. temperature-seasonality associations in the auxiliary involucre/spine-like image proxies.

These summaries preserve their original limitations.

## What can be handed off

A valid downstream target records:

- source cohort and endpoint family;
- numerical estimate or interval;
- observation level / scale;
- whether it is primary, auxiliary or sensitivity evidence;
- original claim boundary.

Example:

```text
Azami observation:
  involucre_projection_roughness × BIO4 beta = +0.0975

Downstream use:
  pattern requiring a mechanism model to permit a positive seasonality association

Not allowed:
  beta = defence-selection coefficient
  BIO4 caused stronger spines
  roughness is a botanically validated defence trait
```

## What stays out of Azami

The following belong downstream in EAzami or focal studies:

- pollinator-benefit functions;
- florivore / seed-predator cost functions;
- seed-output effect-size meta-analysis;
- trait-to-fitness manipulations;
- population ancestry;
- regulatory / molecular mechanisms;
- generative or ABC simulations comparing causal mechanism families.

## Cross-layer reduction logic

The intended programme is:

```text
Azami
observed global phenotype/environment patterns
        ↓ frozen target summaries

independent interaction literature
pollination + antagonist + fitness effect sizes
        ↓

EAzami
minimal mechanism models / simulation
        ↓
ask whether the same mechanism family can reproduce both evidence layers
        ↓

focal field data
replace external priors and discriminate remaining mechanisms
```

A simulation success does not retroactively upgrade Azami Chapter 1 to a causal study. A simulation failure likewise does not invalidate the measured Azami association; it rejects or weakens the tested mechanism family.

## Stop rule

Do not modify the Chapter 1 headline or submission claims because a downstream simulation reproduces them. Chapter 1 changes only if its own measurement, cohort, robustness or external-validation evidence changes.

# Floral trait architecture field protocol

## Purpose

This protocol supports the `azami` programme by measuring head orientation
alongside floral display, floral barriers, abiotic exposure, interactions, and
reproductive outcomes on shared biological units. Its purpose is to distinguish
plausible functional routes, not to relabel trait correlations as adaptation.

## Study-unit hierarchy

Keep the following identifiers at every stage:

```text
species -> population -> site -> plant -> capitulum -> observation period
```

The focal unit for trait-to-interaction analysis is usually a capitulum nested
within plant and population. Record repeated observation periods separately.

## Minimum panel

For every focal capitulum, record:

1. **Trait state:** orientation class, head angle, anthesis stage, head diameter,
   display height, and the predeclared defensive or exposure traits relevant to
   the system.
2. **Pollination route:** observation effort, visitor identity to the finest
   feasible taxon, visit count, and legitimate anther/stigma contact when it can
   be observed.
3. **Antagonism route:** floral damage score, antagonist identity when observed,
   and whether damage preceded or followed the interaction observation.
4. **Abiotic route:** wetness or rainfall state, wind class, and a standardised
   photograph documenting head exposure.
5. **Outcome:** fruit set and/or filled seed production. Add mating-system or
   paternity data only under an independent genetic protocol.

## Timing and standardisation

- Measure morphological traits at a declared anthesis stage.
- Photograph every focal head with a scale and an orientation reference.
- Record visitors as effort-standardised head-hours; do not compare raw visit
  counts across unequal observation durations.
- Score floral damage before terminal fruit collection whenever possible.
- Record weather at each observation, because head orientation can change the
  relationship between exposure and interaction outcomes.
- Keep manipulated and unmanipulated heads explicitly labelled; never merge
  them during exploratory summaries.

## Candidate functional routes

The following routes are tested separately before any integrated conclusion:

```text
orientation/display -> visitation or legitimate handling -> reproductive outcome
orientation/barrier -> floral damage -> reproductive outcome
orientation/barrier -> pollinator obstruction -> reproductive outcome
orientation/cover -> rain or wind exposure -> reproductive outcome
```

A positive association between a trait and seed set does not identify any one
route. A visit count alone does not establish pollination efficiency, and a
trait–damage association does not establish defensive efficacy without a
predeclared comparison or experiment.

## Analysis boundary

1. Use mixed or hierarchical models that retain population and plant structure.
2. Treat floral traits as potentially correlated; report collinearity and avoid
   making causal claims from a single selected predictor.
3. For comparative analyses, incorporate phylogeny or make the non-phylogenetic
   scope explicit.
4. Analyse floral and leaf modules separately unless a predeclared cross-organ
   bridge is measured.
5. Report null and ambiguous pathway results, not only supporting patterns.

## Relationship to general theory

The panel can be mapped to the general attraction–barrier framework in
`biotic-interaction-trait-architecture`, but `azami` remains the biological home
for Cirsium-specific assumptions, traits, protocols, and data. This mapping does
not turn the field panel into evidence of a universal trade-off.

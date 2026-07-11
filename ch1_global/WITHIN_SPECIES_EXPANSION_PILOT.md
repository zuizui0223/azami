# Within-species expansion pilot

## Why

The current global atlas caps each species at 40 observations and each 5-degree
block at three observations. That design protects between-species comparisons
from domination by common species, but it is underpowered for within-species
trait–climate inference.

The expansion pilot therefore creates a second, independent sampling layer. The
40-observation atlas remains frozen for between-species analysis.

## Pilot design

Default eligibility:

- research-grade species-level iNaturalist observations;
- non-captive;
- public coordinates usable for environment extraction;
- positional accuracy <=10 km;
- at least 60 eligible observations per species;
- at least eight 2-degree spatial blocks;
- up to 250 selected observations per species;
- up to 15 observations per species per block;
- one primary photo per observation;
- no measured trait or flower annotation used for selection.

Selection proceeds round-robin across spatial blocks before taking additional
observations from the same block. This prioritizes geographic and climatic spread
rather than simply retaining the most photographed localities.

## Purpose

The pilot asks whether increasing within-species density changes:

1. usable observations after detector and trait QC;
2. within-species climatic range;
3. precision of demeaned climate coefficients;
4. consistency between full and <=10 km coordinate cohorts;
5. power to detect small and moderate standardized effects.

It does not replace the 216-taxon between-species dataset and does not make a
trait–climate association causal.

## Recommended first run

Start with the default eligibility criteria and inspect
`within_species_expansion_species_audit.csv`. Do not immediately process every
eligible species. Select approximately 10–20 taxa spanning continents, head
orientation and existing sample coverage, then run detector and continuous-trait
measurement using the existing scripts 54–57.

The pilot should continue only when most selected taxa retain:

- at least 50 QC-usable observations for the relevant endpoint;
- at least eight independent spatial blocks;
- non-trivial BIO1/BIO4/BIO12/BIO15 ranges;
- no single block contributing a large majority of observations.

## Interpretation

A stronger or newly detected coefficient after expansion would show that the
40-observation atlas limited detectability. A result that remains weak after
high-density, broad-gradient sampling would provide more informative evidence
than the current null result, but still would not demonstrate absence of a
biological effect.

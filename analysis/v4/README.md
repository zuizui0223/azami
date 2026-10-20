# Chapter 1 v4: one workflow from source to functional-trait inference

PR142 is exploratory. The frozen v3 results are a reference, not the conclusion
of the redesigned v4 analysis. Execution and verification belong in GitHub Actions.

## Inputs and reuse

- Reuse the published [Zenodo v3 numerical inputs](https://doi.org/10.5281/zenodo.22703215)
  through the existing `Reproduce current Chapter 1 analysis` Action. Its eight
  stages and 16 reference outputs verify the inherited numerical baseline.
- Zenodo does not contain the full 665,115-observation upstream metadata or
  original photographs. Do not silently substitute the old 46,276-observation
  cohort for the full source.
- The existing `V4 all 665115 upstream native climate replay` Action restores the
  exact original source from the repository secret `AZAMI_ORIGINAL_METADATA_URL`,
  an owner-authorized HTTPS download URL. The SHA-256 remains fixed. The former
  Actions artifact expired; a normal private Drive sharing link is not a direct
  authenticated download for GitHub runners.
- Reuse existing WCVP evidence, CHELSA raster caches and valid archived image
  measurements. Missing or unusable measurements remain explicit. Reacquisition
  is for missing coverage or measurements whose definition must change.

## Primary analysis sequence

1. Start with all source observations and apply native-range eligibility before
   selecting on any measured trait. Report source, native, measured and analysed
   counts separately. See [full-source audit](NATIVE_665115_FULLSOURCE_20261020.md).
2. Analyse colour, orientation and phyllary opening separately. Each uses all
   available observations for that response, without requiring the other trait
   groups or old involucre shape descriptors. Do not average responses to species
   or spatial-cell means as the primary analysis.
3. For each response, jointly include temperature, radiation, VPD and rainfall;
   inspect correlations and support before interpretation. Separate within-species
   environmental deviations from between-species exposure and account for repeated
   observations by species. A species random intercept alone is not this separation
   and is not a phylogenetic correction.
4. Evaluate species influence. Broad spatial structure and phylogenetic placement
   are subsequent checks; soil is an additional comparison on matched records.
5. Restrict only the supplementary whole-head integration analysis to records with
   all three trait groups. Evaluate within- and between-species covariation separately.
6. Report effect estimates, uncertainty, measurement coverage and failed or
   non-evaluable hypotheses together. Environmental association tests predictions
   of a functional hypothesis; it does not establish adaptation or a mechanism.

## Measurement and execution limits

Whole colour uses the observed lightness, chroma and circular hue information;
it is not saturation alone. The archived hue definition must be preserved rather
than relabelled as CIELAB a*/b* or anthocyanin concentration.

Archived bract projection/contour descriptors are not individual phyllary opening
angles. A defence-related opening measure needs identified anatomical reference
points and validation. The existing original-byte anatomical Action replays five
annotated examples; its success does not establish full-source measurement coverage.
See [upstream measurement boundary](UPSTREAM.md).

The processed-CSV `functional_pilot.py`, earlier species/cell models and local
experimental scripts are diagnostic history, not the primary workflow above.
The full-source classification Action currently ends at native climate extraction;
full-source trait measurement, the redesigned models and their ecological conclusions
are not yet completed in Actions. Do not report a green diagnostic job as completion.

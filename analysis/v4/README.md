# Chapter 1 v4: one workflow from source to functional-trait inference

PR142 is exploratory. The frozen v3 results are a reference, not the conclusion
of the redesigned v4 analysis. Execution and verification belong in GitHub Actions.

## Goal and completion criteria

The active goal is to explain what quantitative photographic traits make
comparable beyond species categories, then evaluate environmental predictions
about capitulum function. Reuse existing data and processing instead of adding
parallel analysis routes. Completion requires one Actions-backed report with
trait-specific sample sizes, jointly adjusted environmental estimates and
uncertainty, within/between-species covariation, species influence, additional
soil/spatial/phylogenetic checks, and explicit unsupported or non-evaluable
functional claims. A successful data-preparation job alone does not meet the goal.

Verified starting point: [run 39407663676](https://github.com/zuizui0223/azami/actions/runs/39407663676)
reproduced 665,115 source observations, 313,350 native-eligible observations and
313,240 complete nine-predictor records. [Zenodo replay 39405950824](https://github.com/zuizui0223/azami/actions/runs/39405950824)
matched all 16 inherited numerical reference outputs. These establish reusable
inputs and a historical baseline, not the redesigned v4 biological conclusions.

The next input gate passed in [run 39408603488](https://github.com/zuizui0223/azami/actions/runs/39408603488):
all 313,350 native rows remain, including 96,246 without a full-source archived
measurement row. All 27 endpoint columns and their separate eligibility flags
are retained. No extra trait averaging or cross-trait complete-case filter is
applied. With the four common climate predictors, orientation has 170,066 usable
observations, lightness/chroma 186,439, circular hue 186,390, and archived
involucre/projection descriptors 1,303. These last descriptors are not opening
angles. The colour cohort has temperature–VPD correlation 0.91045 and VPD VIF
10.21468; joint effects require uncertainty and collinearity diagnostics.

The 1,303 phyllary-proxy records are **inherited coverage, not the result of
measuring all native photographs**. Published cohort and measurement receipts
trace 406,582 detector-positive observations → 203,020 with qualifying
high-resolution heads → 22,186 after intersecting the old strict spatial cohort
and retaining one head per observation → 22,177 reconstructed crops → 2,326
usable contour measurements, of which 1,303 enter the current native analysis.
These historical stages use different denominators from the original 665,115
source records. The existing manual Action now audits per-endpoint failure
dispositions from the checksum-verified Zenodo archive before fitting models.
Failure codes describe the old processing method, not the inherent measurability
of the original photograph. Full native anatomical measurement remains unfinished.

Downstream input: encrypted artifact `12257794094`,
`v4-native-observation-input-encrypted`; checksum/provenance and endpoint counts
are in aggregate artifact `12257893990`. The prepared table contains every native
observation and the original per-observation values, including ineligible values
with false eligibility flags. It must be filtered by each response's own flags at
model fitting, not by completeness across all 27 endpoints. Its success is input
readiness, not a fitted v4 result.

## Inputs and reuse

- Reuse the published [Zenodo v3 numerical inputs](https://doi.org/10.5281/zenodo.22703215)
  through the existing `Reproduce current Chapter 1 analysis` Action. Its eight
  stages and 16 reference outputs verify the inherited numerical baseline.
- Zenodo does not contain the full 665,115-observation upstream metadata or
  original photographs. Do not silently substitute the old 46,276-observation
  cohort for the full source.
- The existing `V4 all 665115 upstream native climate replay` Action restores the
  exact original source from an AES-256-GCM encrypted asset in the existing draft
  Release `private-v3-numerical-20260908`. The separate `AZAMI_METADATA_KEY` secret
  decrypts it only inside the runner. Both encrypted and original byte hashes are
  checked; only aggregate receipts are uploaded. The former Actions artifact
  expired. The original Drive backup remains private and unchanged.
- Reuse existing WCVP evidence, CHELSA raster caches and valid archived image
  measurements. Missing or unusable measurements remain explicit. Reacquisition
  is for missing coverage or measurements whose definition must change.

## Primary analysis sequence

Primary model specification: [observation-level model contract](OBSERVATION_MODEL_CONTRACT.md).
The existing `V4 functional trait analysis` Action runs these models on manual
dispatch, reusing the verified encrypted observation input. PR execution retains
the historical pilot as a regression check. These are distinct execution modes.
Manual execution first calls the existing original-photo method workflow, then
audits inherited coverage and fits the archived-endpoint reference models.
The new measurement coordinates are not silently substituted into inherited
columns. A green demonstration replay does not certify botanical validity or
native-wide completion; the numerical report explicitly retains incomplete
new-measurement status until validated annotations and full-source execution exist.
GEB is the intended journal; the scientific contribution must concern general
macroecological patterns and what continuous within-species photographic
variation reveals beyond species categories.

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

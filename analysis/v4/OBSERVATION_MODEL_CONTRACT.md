# V4 environmental response workflow

Latest model instruction: use raw observation environmental values simultaneously
with a species random intercept, without within/between decomposition in the main
model. `environment_response.py` now uses `fit_module(..., decompose=False)`.
The eight selected predictors remain unchanged; colour and orientation are main,
bract projection is SI. One family contains 3 x 8 = 24 joint tests. Coefficients
are common slopes informed by both within- and between-species variation, under
the random-intercept model assumptions, and cannot be labelled within-only or
between-only effects. The successful 39416125522 decomposed analysis remains
archived as a different estimand; its 48-test correction is not rewritten.
The decomposition description below documents that archived fit.

Question: which environmental gradients accompany changes in capitulum traits,
conditional on the other included environmental dimensions?

## Inputs and scope

Reuse the checksum-verified observation table from Actions run 39408603488
(encrypted artifact 12257794094; aggregate receipt 12257893990). The source
contains 665,115 observations; native-range classification yields 313,350.
Do not remeasure photographs or average trait responses. Use each endpoint's
own eligible observations, including species with fewer than five observations.
Do not impose a complete-case filter across trait modules.

## Select environmental variables before fitting traits

Use all nine stored variables from the six groups defined in
`analysis/ch1/capitulum_environment_blocks_contract.json`:

| Group | Candidates |
|---|---|
| Temperature | BIO1, BIO4 |
| Precipitation | BIO12, BIO15 |
| Radiation and atmospheric drying | radiation, VPD |
| Mechanical exposure | wind speed |
| Growing-season water input | growing-season precipitation |
| Climatic productivity | NPP |

`prepare_vifstep.py` verifies the input and exports only these environmental
columns. Selection uses all native rows with all nine values finite, without
trait values, measurement eligibility, response associations or taxon averaging.
`select_environment.R` calls actual `usdm::vifstep` with threshold 10, no `keep`
variables and `size=nrow(x)+1` (no subsampling). Do not force each group to survive.
Save candidate coverage, input identity, exclusion order, per-step and final VIF,
correlations, package versions and code/run identifiers.

Run 39415814454 selected eight variables on 313,240 rows, excluding BIO1 only.
This is environmental redundancy selection, not proof of causal identification.
Do not change its threshold or protected variables after inspecting trait results.

## Fit the selected set unchanged

`environment_response.py` reads the selection artifact and verifies that it refers
to the same input. All three modules use exactly the selected variables:

- Whole visible colour: archived L*, C*, hue sine and cosine jointly tested.
- Presentation: image-vertical angle; higher means more downward in the image
  reference, not an independently measured gravity angle.
- Bract contour projection: relative 95th-percentile projection and spread
  fraction jointly tested; not anatomical opening angles or measured defence.

For every selected variable, compute species means from all native observations
with that variable measured, and deviations from those means. Retain individual
response observations. Center/scale exposures on the native population, before
trait filtering. Fit response-specific Gaussian random species intercept REML
models with all selected within/between terms simultaneously. Joint tests use
stacked species-cluster scores, not an unrestricted multivariate likelihood.
Species clustering allows within-species dependence but is not phylogenetic or
spatial correction. It does not establish a common slope in every species.

Use one Holm family of 3 modules x 2 scales x number of selected variables
(48 tests for eight variables). Failed tests remain unidentifiable and count as
p=1 for multiplicity; do not interpret them as evidence of no association.
Report effect sizes per native environmental SD, pointwise 95% intervals,
endpoint-specific N/species counts and exact within/between design diagnostics.
Selection on raw environmental rows does not guarantee low VIF after decomposing
within/between or restricting to measured endpoints; report this transparently,
without an outcome-driven second selection step.

## Execution and evidence

Run `.github/workflows/v4-functional-exploration.yml` on the PR branch. Default
manual dispatch stops after environmental selection. `fit_after_selection=true`
executes selection followed by the three trait models and verifies that the
selected set is unchanged. Tests and scientific computation run in Actions.
The public artifact contains aggregate JSON, effect CSV, inventory and report;
observation identifiers, coordinates and private input tables are not uploaded.

The question is conditional association, not adaptive mechanism, causal effect,
plasticity, pigment identity, defence or genetic modularity. Earlier results have
been seen, so this is retrospective analysis, not preregistered confirmation.
No omission-model search, additional spatial/soil/phylogenetic models or whole-head
integration are part of this narrowed execution path. Historical code/results
remain available in Git and prior artifacts and are not selected for significance.

# Environmental scope and selection for v4

The estimand is an observational environmental gradient conditional on the other
retained exposures, separately within and between species. It is not the causal
effect of changing one climate variable. Existing results have already been seen;
this is a transparent revised analysis, not a preregistered independent test.

## Candidate universe

The nine available candidates come from the six groups in
`analysis/ch1/capitulum_environment_blocks_contract.json`. That inherited
contract defines the candidate universe, not the v4 statistical test.

| Candidates | Environmental description | Interpretation boundary |
|---|---|---|
| BIO1, BIO4 | Mean temperature and thermal seasonality | Seasonality cannot substitute for mean temperature |
| BIO12, BIO15 | Annual precipitation and its seasonality | Annual rain is not rain during flowering or direct floral wetting |
| rsds, VPD | Shortwave radiation and atmospheric drying | Shortwave radiation is not UV; VPD is not soil moisture |
| sfcwind | Mean near-surface wind exposure | Candidate mechanical/drying context; not head-level loading or a demonstrated confounder |
| GSP | Precipitation in the climatic growing season | Shares water-input information with annual precipitation |
| NPP | Modelled potential productivity | Composite climatic index, not measured local productivity or an independent mechanism |

These groups are practical descriptors of the available data, not six independent
biological mechanisms. Wind was inherited as the mechanical-exposure candidate;
its inclusion was not motivated by its observed trait P value. NPP and GSP may
share or derive information from other climate variables. Interpreting coefficients
as independent causal drivers is therefore inappropriate.

## Executed rule

1. Use all native records with complete candidate environments, without response
   access or measurement eligibility filters.
2. Run actual `usdm::vifstep`, threshold 10, no protected variable, and no random
   subsampling. Save all exclusion steps, correlations and package versions.
3. Freeze the retained predictor list for every trait module.
4. Split every retained exposure into its native-observation species mean and
   each observation's deviation. Fit both sets simultaneously, preserving each
   endpoint's eligible response rows.
5. Report correlation/VIF/rank diagnostics for each endpoint's between, within
   and combined designs. Report uncertainty rather than pruning further to
   recover significance.

The existing rule excludes BIO1 only. Consequently, the current analysis does
not estimate a mean-temperature coefficient. VIF 10 is a declared redundancy
threshold, not proof of adequate causal adjustment. Selection on raw environments
does not ensure that decomposed or phenotype-specific designs have low VIF.

## Required reporting and remaining audit

- Quantify conditional slope contrasts and pointwise intervals in original
  phenotype units; retain the whole-colour and whole-outline joint tests.
- Show separate between/within exposure support. Pooled within-species exposure
  ranges are not a range supported within every species.
- Record source raster temporal coverage, units, resolution and derivation.
  These provenance details must be reconciled before interpreting physical-unit
  gradients; figures currently identify the stored environmental units explicitly.
- If diagnostics make individual gradients poorly resolved, report that result.
  Any later candidate change needs a documented scientific rationale and a
  separately labelled specification, not replacement of inconvenient results.
- Species random intercepts do not account for spatial dependence or phylogeny.

Method reference: https://babaknaimi.r-universe.dev/usdm/doc/manual.html

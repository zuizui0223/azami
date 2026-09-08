# Chapter 1 v3 scientific design rationale

V3 has two linked aims: make the image-to-trait measurements traceable and test
ecological associations within their measured support. Source definition,
measurement evaluation and the model specification come before interpretation.
V2 results have already been inspected, so this is retrospective redesign, not
preregistration. This document explains the design; it is not an execution receipt.

The active order is source definition, measurement, assessability, primary
within/among inference, then optional breadth synthesis. See the
[integrated contract](integrated_workflow_contract.json) for explicit revisions to
earlier model mappings and the [evidence index](integrated_evidence_index.json) for
what is actually implemented. A successful preflight is not production permission.

## 1. Preserve the source; define the ecological target separately

The master ledger retains every recovered observation-photo link, including
native, introduced and unresolved records. It preserves source taxon names and
the evidence used to resolve them to accepted taxon concepts. Name resolution is
not independent validation of the identification in a photograph. Uncertain,
hybrid and unresolved assignments remain visible in coverage accounting.

Primary ecology concerns native-range records after a versioned WCVP/TDWG join
and question-specific taxonomic, wild-status, date and location requirements.
Introduced and unresolved records remain useful for methodological coverage and
separately labelled range-scope comparisons. Native restriction is an analysis
target, not deletion from the acquisition source or a significance-triggered check.

The original 46,276-observation view and the capped environment-diagnostic sample
do not set the final ecological denominator. Preserve unthinned source records;
report the realized source, image-availability and endpoint-specific denominators.

## 2. Start from biological questions, then select exposure representations

The [capitulum hypotheses](capitulum_abiotic_hypotheses_v3.md) distinguish:

- orientation and presentation in relation to precipitation and wetting exposure;
- joint visible colour in relation to radiation, atmospheric drying, water
  availability and thermal exposure;
- gross shape and architecture in relation to physical exposure, only where the
  image construct has sufficient measurement support;
- wind as secondary presentation or architecture exposure, since photographs do
  not measure peduncle mechanics or three-dimensional drag.

Observation-month precipitation climatology is closer to photographed seasonal
exposure than annual precipitation alone, but it is not the weather experienced
by the individual head. BIO12, BIO18 and growing-season precipitation are broader
alternatives, not rival predictors selected by their trait P values. Temperature
can be relevant without assuming a direct vertical-angle temperature mechanism.
Ecosystem productivity is not interchangeable with a direct physical exposure.

Use the environment-only cohort to assess coverage, Pearson/Spearman correlation,
rank, condition number, VIF and redundancy. Select representations using biological
meaning and environmental support before joining traits. The frozen drying and
thermal alternatives are reported together, not chosen by the more favourable
outcome. Marginal slopes cannot establish independent environmental effects.

The current revision uses the same four-exposure drying and thermal formulations
for each primary module. Precipitation is the focal orientation coefficient;
that focal hypothesis is not a reason to leave correlated exposures unadjusted.
This supersedes the earlier precipitation-only primary orientation model, without
rewriting its historical contract. The choice changes the conditional question;
it does not demonstrate that all confounders have been measured. Removing annual
or growing-season representations also narrows the temporal question, not just
redundant information.

Before fitting, diagnose the realized endpoint/module cohort at both scales and
after nuisance adjustment. A low VIF in a capped source-diagnostic sample cannot
certify identification in the final cohort or in every taxon. Keep singular and
weakly supported contrasts visible as non-estimable or uncertain.

Chapter 1 is abiotic-only. Globally comparable pollinator exposure is not included;
pollinator filtering remains an untested alternative explanation, not a proxy to
add after inspecting an abiotic result.

## 3. Evaluate the measurements before ecological interpretation

YOLO locates visible heads; deterministic functions define continuous image
features. Each endpoint needs a versioned definition, units, required image
evidence, QC states and an interpretation limit. Technical consistency is not
detector truth, physical botanical accuracy or a requirement for new human labels.

Retain crop-shift, resolution, sharpness, replay and context diagnostics with
their measurement links. Apply the declared measurement-support rules without
using environmental coefficients to tune thresholds. Orientation and gross shape
need the specified bbox-stable formulation and its bounded alternatives. Colour
needs paired flower/non-head/green-context values, not only availability counts.
Background associations limit floral specificity; their absence does not establish
calibrated reflectance.

The current legacy chroma function changes its pixel subset with the dominant
colour class. Preserve its raw values and version history. A consistently defined
floral-union chroma is a candidate new measurement, not an automatically qualified
replacement; its image-only evaluation must precede primary colour fitting.

Lightness, chroma and circular hue are interpreted jointly. Low chroma alone cannot
distinguish pale from dark colour or quantify anthocyanin. The four colour fractions
form one closed composition, not four independent biological traits. Absolute head
size is not inferred without a scale reference.

Source images may be streamed without retaining a full image archive, but this
does not permit losing the numerical measurements, source identities or processing
provenance. Cloud streaming is stopped until authorized durable private numerical
storage and its verification are implemented. Public outputs remain scientific
code and aggregate receipts; raw links and private records are not uploaded here.

## 4. Respect observation and taxon dependence

Keep head-to-photo-to-observation-to-taxon links. Multiple heads and photographs
do not become independent observations or automatically identify the same biological
individual. Known shared-photo and exact-content components remain linked.

Within-taxon inference concerns the distribution of conditional taxon slopes and
its partially pooled summary, including between-taxon variation and uncertainty.
Wide-ranging taxa can therefore have their own supported estimates. Sparse taxa
remain in coverage accounting when their slopes are not estimable or individually
interpretable. Source-support rules are fixed before coefficients are inspected.

Among-taxon inference uses one uncertainty-aware phenotype summary and matched,
equal-observation exposure summary per taxon. Raw photograph abundance is not a
biological taxon weight. Opposite supported slopes describe heterogeneity rather
than automatically representing failed replication.

Spatial structure belongs in the primary specification. Residual diagnostics and
bounded geographic uncertainty comparisons describe remaining dependence, not
sequential gates that promote surviving associations. The implemented estimator
must be checked against the intended hierarchy before use on trait outcomes.

Exact dates and hemisphere-specific calendar harmonics belong in the source-to-
model path, alongside endpoint-matched imaging terms. They do not identify
developmental stage. The primary estimand describes the photographed mixture of
stages; conditioning on stage or restricting to anthesis asks a different question.
Partial pooling likewise does not recover missing recording effort. Assessability
compares the eligible native source with attempted, downloaded, detected and
QC-usable records, using fixed source-taxon weights and endpoint/module cohorts.
Unreported processing states must not be called failures or detector negatives.

The central ecological comparison is direct: estimate the difference between
within- and among-taxon associations on matched support and a common measurement
scale, with joint uncertainty. Comparing a significant slope with a nonsignificant
slope is not such a test. Within-group centering separates the two association
scales, but does not by itself establish a causal interpretation
([van de Pol & Wright, 2009](https://doi.org/10.1016/j.anbehav.2008.11.006)).

## 5. Report contributions without promoting association to mechanism

The methodological contribution includes recoverable source coverage, endpoint
definitions, technical sensitivity, comparability and reproducible saved products.
The ecological contribution concerns associations and variation within the
observed sampling and exposure support, not required survival of v2 candidates.

Head/photo/observation/taxon variance describes visible image-observation variation.
It cannot separate biological and photographic variance without suitable matching
and additional evidence. The hierarchical variance analysis remains a planned
analysis until an executed receipt and numerical results exist.

Joint colour and shape modules organize the main questions; individual endpoint
slopes decompose those questions. Whole-capitulum synthesis is secondary and cannot
rescue an unsupported module by searching for a favourable combination.

The reporting structure therefore remains: image-to-trait coverage and technical
evaluation; measured distributions and module covariance; within/among ecological
estimates and their direct contrast; then optional breadth synthesis. The method's
contribution is this traceable distribution-level evidence, not the invention of
YOLO, PCA or hypervolume. The ecological contribution is the scale comparison,
whether its associations are clear, weak, absent or different from v2.

Hypervolume A can visualize environmental coverage at assessability, but does not
replace retention and density diagnostics. Hypervolume B is a secondary comparison
of environmental and phenotype breadth on the same observation sets. Use common
axes and units across taxa, sample-size controls, declared bandwidth and probability
mass, and convergence checks. Different spaces cannot be compared by direct set
overlap. Bandwidth, dimension and sample size affect estimated hypervolumes
([Blonder et al., 2018](https://doi.org/10.1111/2041-210X.12865)); these choices need
outcome-blind qualification. B is optional and is not independent corroboration of
the same-data primary association. Observed environmental breadth is not a complete
fundamental niche, and image-phenotype breadth is not plasticity.

Flora or monograph descriptions can provide a separate taxon-matched external
consistency benchmark, not individual-level physical calibration. Phylogenetic
placement comparisons must retain the same applicable model and report direct-tip
coverage, lambda and coefficient movement. Fifty-two placements are not fifty-two
independent confirmations; negligible leverage should be stated directly.

The claim remains an image-feature association that can motivate functional
validation. Adaptation, mechanism, pigment amount and true gravitational
orientation require evidence not supplied by spatial association alone. Report
weak, absent, non-estimable and unexpected results with their support and limits.

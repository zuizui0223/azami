# Chapter 1 v3 supervisor-comment integration

This ledger maps the substantive supervisor comments in `Azami_Chapter1_Main（チョウ）_ss_v2_0907ABE.docx` to the v3 workflow. It is a redesign ledger, not a rebuttal letter. The goal is to move concerns upstream into source definition, measurement, model specification and claim limits rather than defend selected v2 results with a long post-hoc sensitivity ladder.

## Author-labelled comments (`作成者`) — 7/7

| Comment ID | Comment | v3 integration | Status before ecological fitting |
|---|---|---|---|
| 153 | Introduced records are discussed, but misidentification should also be acknowledged. | Preserve source taxon names; require a versioned resolved taxon concept/accepted-name join; uncertain/hybrid/unresolved assignments do not silently enter primary taxon-level ecology. Name resolution is explicitly not image-ID validation. | Design integrated; taxonomic join still to execute. |
| 283 | Environmental variables are mostly tested individually; show correlations and identify which analysis actually tests added information beyond the climate core. | Primary abiotic inference is block-structured: core thermal/hydric climate, then explicit nested tests of radiation/VPD, wind, growing-season water and productivity beyond the same core/cohort. Correlation matrix, matrix rank, condition number and VIF are required before fitting. Univariate atlas slopes are descriptive only. | Design integrated; v3 ecological models not yet fitted. |
| 352 | `adaptive-pattern candidates` overstates adaptation; use functional-validation language. | `adaptive-pattern candidate` is disallowed primary language. Allowed label is `candidate association for functional validation`; adaptation/mechanism remain outside the claim ceiling. | Integrated. |
| 360 | Orientation is sensitive to bounding-box shifts; how does that instability affect the precipitation result? Can the relationship be shown using stable measurements? | Measurement robustness is upstream. Orientation requires bbox/resolution/sharpness/replay diagnostics and a predeclared stable-measurement subset or explicit uncertainty propagation before environmental outcomes are interpreted. Stability thresholds must be chosen from image-only perturbations, not from the precipitation coefficient. | Design integrated; full-source propagation/stable-cohort ecological fit remains pending. |
| 403 | If the 52 phylogenetic placements barely change the result, are they really informative? | Placement trees are no longer independent evidence votes or a candidate-promotion gate. Report direct-tip coverage, Pagel lambda and coefficient/uncertainty movement; if placement has negligible leverage, state that directly. | Integrated. |
| 493 | The variation partition can be mistaken for biological intraspecific variation; image-observed variation should be distinguished from biological variance. | Replace the single `within-taxon variance` story with a required head→photo→observation→taxon decomposition and use `visible image-phenotype/image-observation variation` terminology. Ecological importance of below-taxon image variation is a downstream hypothesis. | Design integrated; hierarchical variance model remains pending. |
| 562 | Lower chroma is ambiguous because low chroma can occur in both pale/high-lightness and dark/low-lightness colours. | Colour is interpreted jointly: lightness + chroma + circular hue first at the colour-module level; four colour fractions remain one closed composition. Chroma alone cannot be translated into pale/dark or anthocyanin amount. | Integrated; joint v3 colour-environment fit remains pending. |

## Sakaguchi comments that change the scientific workflow

Many of the 72 Sakaguchi comments are wording, citation, figure or terminology edits. The comments below are the ones that alter the inferential workflow.

### 1. Native range and introduced records — comments 17 and 305

V2 treated introduced records in the primary cohort and moved native-only to a later sensitivity. V3 reverses that logic.

- The master source ledger retains native, introduced and unresolved records.
- A versioned native-range join is required before ecology.
- **Primary ecological inference is native-range only.**
- Introduced/unresolved records are retained for source/methodological coverage and a separately labelled range-scope/transportability comparison, not as evidence defining the primary abiotic relationship.

This is an analysis-target decision, not a significance-triggered sensitivity.

### 2. Dominance by very few taxa and very small taxa — comments 444 and 446

V2 was dominated by `Cirsium vulgare` and `C. arvense` at the row level and contained many sparsely sampled taxa. V3 does not solve this by repeatedly deleting the top taxa after significance.

- Multiple photos/heads do not gain independent observation weight.
- Within-taxon inference uses taxon-specific **partial-pooled random slopes** rather than one raw pooled row-level slope as the sole primary estimand.
- The primary result is the distribution of taxon slopes and its hyperdistribution, with sample/environment-range support shown per taxon.
- Among-taxon inference uses one uncertainty-aware taxon summary per taxon; raw photograph abundance is not a biological weight.
- Taxa lacking sufficient independent observations/environmental range are retained in coverage accounting but are not individually interpreted. Exact support thresholds must be frozen before v3 ecological outcomes are inspected.

### 3. Wide-ranging taxa: do individual taxa show environment–trait relationships? — comment 463

This becomes part of the primary within-taxon model rather than a follow-up plot. The hierarchical random-slope model returns taxon-specific slope posteriors/estimates, between-taxon slope variation, sample support and environmental-range support.

### 4. Environmental variables should be coherent, not nine separate stories — comments 254, 262 and the `作成者` comment 283

The primary sequence is:

1. diagnose covariance/collinearity on the exact cohort;
2. fit the thermal+hydric climate core;
3. test each predeclared abiotic block for information beyond that core;
4. only then decompose a supported module/block pattern to endpoint coefficients.

The individual predictor atlas can remain descriptive, but it is not the primary independent-effect test.

### 5. BIO18 / seasonal precipitation — comment 266

BIO18 is now a predeclared alternative representation of warm-season precipitation. It must be compared against the existing growing-season-water formulation using the same cohort/response and cannot be retained merely because it gives a smaller P value.

### 6. Pollinator filtering — comment 274

The project decision is **abiotic-only Chapter 1**. No directly comparable global pollinator exposure layer is available. Pollinator filtering is therefore outside the estimand and can appear only as an untested alternative mechanism in Discussion. Its absence is not a submission/analysis-completion gate and no coarse post-hoc pollinator proxy will be added to defend an orientation result.

### 7. Clustered SE / non-independence — comment 277

The v2 clustered-SE explanation is superseded as the main design. V3 keeps head/photo/observation/taxon nesting explicitly and uses hierarchical models; known exact-image/shared-observation dependence is preserved in component groups.

### 8. Within/between variance meaning and figure logic — comments 221, 241, 243, 244, 496, 504, 720

The primary question is no longer a two-bar AMOVA-like claim that the residual is biological within-species variation. V3 decomposes variance across head/photo/observation/taxon levels and reports image-observation variation separately from biological interpretation.

### 9. Trait extraction accuracy/comparability and YOLO role — comments 195, 198, 202, 353

- YOLO performs localization only.
- Deterministic functions produce traits.
- Every endpoint has a versioned definition, required image evidence, QC, units and interpretation ceiling.
- Perturbation/replay/resolution/sharpness/context diagnostics are retained per endpoint.
- Technical repeatability does not equal detector truth or physical botanical accuracy.
- Absolute head size is not inferred without a scale reference (comment 615).

### 10. Trait–trait correlation / whole-capitulum structure — comments 340 and 477

V2 treated the complete-18 synthesis as secondary. V3 promotes trait structure upstream: orientation, visible colour, gross shape and architecture/fine geometry are coherent modules. Module/multivariate structure is tested before turning every endpoint×environment pair into a separate biological story; endpoint coefficients are decomposition of supported structure.

### 11. Literature/field-guide benchmark — comment 586

Authoritative flora/monograph descriptions remain a separate external consistency benchmark where taxon concepts can be matched. They are not substituted for physical individual-level calibration and cannot be used to manufacture significance by changing taxon coverage.

### 12. Phylogenetic procedure/reproducibility — comments 337 and the `作成者` comment 403

Phylogeny is no longer a sequential `52/52 survived` evidence ladder. The same ecological model may be refit across defensible placements, but the report emphasizes direct-tip coverage, lambda and coefficient movement. Placement sensitivity is not 52 independent confirmations.

### 13. Scope, data source, geography and naming — comments 157, 160, 162

The v3 master ledger records exact source snapshots, source taxon assignment, observation-photo links, geographic support and taxonomic-resolution provenance before any ecology. The final manuscript still needs a clean prose statement of the realized genus/taxon and geographic scope after the v3 taxonomic/range joins are executed.

## Workflow shape after integration

```text
full source ledger
  ↓
taxon-concept + native-range + date/location/support annotations
  ↓
image localization → endpoint-specific measurement + technical uncertainty
  ↓
head → photo → observation → taxon hierarchy
  ↓
native-range ecological analysis view
  ↓
trait modules × abiotic environmental blocks
  ├─ within taxon: partial-pooled taxon slopes
  └─ among taxa: uncertainty-aware taxon summaries
  ↓
spatial structure included in the primary specification
  ↓
limited bounded sensitivities (range transportability, measurement formulation,
phylogenetic placement where applicable)
  ↓
association for functional validation, not adaptation/mechanism
```

The v3 design intentionally avoids the v2 pattern of: find FDR-positive rows → delete dominant taxa → native-only check → spatial gate → 52-tree gate → call survivors stronger candidates. Controls that define the estimand or error structure are moved into the main design. Remaining sensitivities quantify bounded assumption changes; they do not accumulate as evidence votes.

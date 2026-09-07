# Chapter 1 v3 supervisor-comment integration

This ledger maps the substantive supervisor comments in `Azami_Chapter1_Main（チョウ）_ss_v2_0907ABE.docx` to the v3 workflow. It is a redesign ledger, not a rebuttal letter. The goal is to move concerns upstream into source definition, biological hypotheses, measurement, environmental representation, model specification and claim limits rather than defend selected v2 results with a long post-hoc sensitivity ladder.

## Author-labelled comments (`作成者`) — 7/7

| Comment ID | Comment | v3 integration | Status before ecological fitting |
|---|---|---|---|
| 153 | Introduced records are discussed, but misidentification should also be acknowledged. | Preserve source taxon names; require a versioned resolved taxon concept/accepted-name join; uncertain/hybrid/unresolved assignments do not silently enter primary taxon-level ecology. Name resolution is explicitly not image-ID validation. | Design integrated; full-source taxonomic/native-range join code added, execution still required. |
| 283 | Environmental variables are mostly tested individually; show correlations and identify which analysis actually tests added information beyond the climate core. | The former privileged `climate core` has been removed from v3 design. Biological capitulum×abiotic hypotheses are specified first. Candidate exposures then enter an **environment-only, phenotype-blind matrix** on equal footing. Coverage, Pearson/Spearman correlation, matrix rank, condition number, VIF and redundancy components are diagnosed before any trait join. The final minimal environmental representation is frozen from biological meaning + environment-only redundancy, not from endpoint P values. Individual predictor slopes are decomposition/descriptive outputs, not independent-effect discovery tests. | Design and code integrated; CHELSA source pilot queued/running, final native-range environment matrix still to execute. |
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
- A versioned WCVP/TDWG native-range join is required before ecology.
- **Primary ecological inference is native-range only.**
- Introduced/unresolved records remain available for source/methodological coverage and a separately labelled range-scope/transportability analysis, not as evidence defining the primary abiotic relationship.

This is an analysis-target decision, not a significance-triggered sensitivity. The full-source join has its own v3 contract and builder and does not reuse the old 46,276-row native-only result as the new denominator.

### 2. Dominance by very few taxa and very small taxa — comments 444 and 446

V2 was dominated by `Cirsium vulgare` and `C. arvense` at the row level and contained many sparsely sampled taxa. V3 does not solve this by repeatedly deleting the top taxa after significance.

- Multiple photos/heads do not gain independent observation weight.
- Within-taxon inference uses taxon-specific **partial-pooled random slopes** rather than one raw pooled row-level slope as the sole primary estimand.
- The primary result is the distribution of taxon slopes and its hyperdistribution, with sample/environment-range support shown per taxon.
- Among-taxon inference uses one uncertainty-aware taxon summary per taxon; raw photograph abundance is not a biological weight.
- Taxa lacking sufficient independent observations/environmental range remain in coverage accounting but are not individually interpreted. Exact support thresholds must be frozen from source support and measurement precision before v3 ecological outcomes are inspected.

### 3. Wide-ranging taxa: do individual taxa show environment–trait relationships? — comment 463

This becomes part of the primary within-taxon model rather than a follow-up plot. The hierarchical random-slope model returns taxon-specific slope estimates/posteriors, between-taxon slope variation, sample support and environmental-range support.

### 4. Environmental variables should be coherent, not nine separate stories — comments 254, 262 and the `作成者` comment 283

V3 no longer starts from a fixed thermal/hydric `core`. The sequence is now:

1. define biological hypotheses from capitulum function/exposure **before choosing variables**;
2. align direct climatological exposures to the observation month where possible (`pr`, `tas`/`tasmax`, `rsds`, `vpd`, `sfcWind`; PET/CMI only as candidate representations);
3. construct the intended native-range **environment-only** cohort without reading trait values;
4. report coverage and diagnose Pearson/Spearman correlation, rank, condition number, VIF and redundancy components;
5. compare alternative representations within the same biological exposure family;
6. freeze the smallest biologically interpretable, non-redundant environmental representation;
7. only then join capitulum traits and fit within/among ecological models.

Annual BIO variables remain broad comparability descriptors where useful, not an a priori privileged core. The individual predictor atlas can remain descriptive, but it is not the primary independent-effect test.

### 5. BIO18 / seasonal precipitation — comment 266

BIO18 and GSP are no longer rival trait predictors chosen by which produces a stronger association. The primary hydrometeor hypothesis uses **long-term precipitation climatology matched to the photographed observation month**. BIO18, GSP and BIO12 are broader seasonal/annual representations. Before trait fitting they are compared against observation-month precipitation and each other for coverage, redundancy and biological correspondence. If redundant, keep the smallest representation that best matches the hypothesized exposure; do not retain whichever yields the smaller trait P value.

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

V2 treated the complete-18 synthesis as secondary. V3 promotes trait structure upstream: orientation, visible colour, gross shape and architecture/fine geometry are coherent measurement/phenotype modules. Trait-specific biological hypotheses are defined first; whole-capitulum multivariate synthesis then asks whether responses collapse to one exposure syndrome or remain module-specific. Endpoint coefficients decompose a supported module/environment relationship rather than replacing it with endpoint×raster significance hunting.

### 11. Literature/field-guide benchmark — comment 586

Authoritative flora/monograph descriptions remain a separate external consistency benchmark where taxon concepts can be matched. They are not substituted for physical individual-level calibration and cannot be used to manufacture significance by changing taxon coverage.

### 12. Phylogenetic procedure/reproducibility — comments 337 and the `作成者` comment 403

Phylogeny is no longer a sequential `52/52 survived` evidence ladder. The same ecological model may be refit across defensible placements, but the report emphasizes direct-tip coverage, lambda and coefficient movement. Placement sensitivity is not 52 independent confirmations.

### 13. Scope, data source, geography and naming — comments 157, 160, 162

The v3 master ledger records exact source snapshots, source taxon assignment, observation-photo links, geographic support and taxonomic-resolution provenance before any ecology. The final manuscript still needs a clean prose statement of the realized genus/taxon and geographic scope after the v3 taxonomic/range joins are executed.

## Biological hypothesis layer added before environment selection

The environmental design is now downstream of the capitulum hypotheses in `capitulum_abiotic_hypotheses_v3.md`:

- **H1 orientation / presentation:** flowering-period precipitation and wetting exposure are the strongest directional abiotic hypothesis. Monthly precipitation aligned to the observation month is the primary representation candidate; BIO12/BIO18/GSP are broader alternatives.
- **H2 visible colour:** joint lightness/chroma/hue state is tested against radiation, atmospheric drying/water availability and thermal exposure; chroma alone has no privileged directional prediction.
- **H3 head/involucre architecture:** desiccation/radiation exposure is tested only for image constructs whose biological meaning and technical robustness are adequate.
- **H4 wind:** secondary presentation/architecture exposure because current images do not measure peduncle mechanics or 3D drag directly.
- **H5 thermal context:** biologically relevant to the whole capitulum and colour; a specific vertical-angle temperature mechanism is not assumed from azimuthal orientation literature.

NPP is no longer given equal mechanistic status with direct physical exposures merely because it existed in v2.

## Workflow shape after integration

```text
FULL SOURCE LEDGER
  ├─ taxon concept + native/introduced/unresolved + date/location support
  │      ↓
  │   NATIVE-RANGE CANDIDATE OBSERVATIONS
  │      ↓
  │   PHENOTYPE-BLIND ENVIRONMENT MATRIX
  │   month-aligned pr / temperature / radiation / VPD / wind
  │   + broader alternatives (BIO12 / BIO18 / GSP etc.)
  │      ↓
  │   coverage + correlation + rank + condition number + VIF + redundancy
  │      ↓
  │   FREEZE ENVIRONMENT REPRESENTATION
  │
  └─ image localization → endpoint-specific measurements + technical uncertainty
         ↓
     head → photo → observation → taxon hierarchy

ONLY AFTER BOTH BRANCHES ARE FROZEN:
  native-range nested phenotype × environment join
      ↓
  hypothesis-linked trait modules × selected abiotic exposures
      ├─ within taxon: partial-pooled taxon slopes
      └─ among taxa: uncertainty-aware taxon summaries
      ↓
  spatial structure in the primary specification
      ↓
  bounded alternatives only
  (range transportability / predeclared measurement formulation /
   phylogenetic placement where applicable)
      ↓
  association for functional validation, not adaptation/mechanism
```

The v3 design intentionally avoids the v2 pattern of: find FDR-positive rows → delete dominant taxa → native-only check → spatial gate → 52-tree gate → call survivors stronger candidates. Controls that define the estimand, exposure representation or error structure are moved into the main design. Remaining sensitivities quantify bounded assumption changes; they do not accumulate as evidence votes.

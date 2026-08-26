# Bias-control reanalysis protocol

## Current decision

Chapter 1 is on **submission hold**. The frozen analysis and its provenance remain intact, but the labels `computational robustness complete` and `submission ready` are not used while the controls below remain unresolved.

This protocol was fixed before executing the raw-calendar seasonal, dominant-taxon and image-resolution outcomes. Separate hemisphere-season and native-range contracts were subsequently frozen before their respective outcomes were inspected. None permits post-outcome threshold changes or rescue analyses. Exact machine-readable rules are in `analysis/ch1/bias_control_contract.json`, `analysis/ch1/hemisphere_season_sensitivity_contract.json` and `analysis/ch1/native_range_sensitivity_contract.json`.

## 1. Seasonal composition

All 46,276 observations in the strict spatial cohort contain a parseable `observed_on` date. For each of the nine endpoints and four CHELSA predictors, the frozen taxon-mean model is reproduced and compared with a taxon-specific cyclic-season model. The latter residualizes both outcome and predictor on a taxon-specific intercept plus sine and cosine of day-of-year before fitting the one-predictor within-taxon model with taxon-clustered uncertainty.

This adjustment addresses cyclic collection timing. It does not identify bud, anthesis or post-anthesis heads and therefore cannot by itself close the developmental-stage gate.

Because the raw calendar assigns opposite seasons to the same date across hemispheres, the separately locked sensitivity uses two definitions: adding one half-cycle to Southern Hemisphere dates, and fitting separate taxon-specific intercept plus sine/cosine curves by hemisphere where both occur. Both definitions reuse the 36-row BH family and outcome-blind omission set. All four non-circular rows retained their frozen sign, BH support and omission-sign stability under both definitions; the developmental-stage gate remains open.

## 2. Dominant-taxon concentration

Omission taxa are selected only from observation counts, without reading trait outcomes. The ten most represented taxa are removed one at a time, and the two most represented taxa are also removed jointly. Each omission reruns all 36 component models under both the baseline and seasonal definitions, with BH correction kept separate by model and omission scenario.

## 3. Visible-colour negative controls

The production colour pipeline uses uncalibrated JPEG pixels and does not retain a non-flower climate control. A new image pass must retain crop-border and nonfloral-green Lab statistics. Specificity must be tested as a paired `region × climate` contrast; a significant corolla slope combined with a nonsignificant background slope is not, by itself, evidence of a difference.

The chroma–BIO12 row was subsequently withdrawn from the current headline after failing the native-only BH rule. The paired colour control remains necessary before broader image-colour precipitation patterns can be interpreted as corolla-specific.

## 4. Involucre resolution and sharpness

The high-resolution auxiliary models are rerun with log minimum head dimension and log1p sharpness as within-taxon covariates. BIO4 must remain positive and BH-supported in the full adjusted family and positive within every predeclared resolution stratum with sufficient data: 150–199, 200–299 and at least 300 pixels.

## 5. Repeated photographs

The source metadata are used to recover every available photograph for observations in the strict cohort. The intended hierarchy is taxon → observation → photograph → head. The between-photograph component is described as repeat-photograph variance, not pure measurement error, unless manual review verifies that the photographs show the same individual.

## 6. Native-range sensitivity

Native and introduced regions require an authority-backed, versioned regional source. The contract pinned WCVP checklist DOI 10.15468/6h8ucr, its exact dataset timestamp and licence, plus TDWG level-3 commit `52da7828aba9d461dd133c27b3bd7a4407161f54` and the geometry SHA-256. No native status was inferred from prevalence or taxon identity. Ambiguous accepted-name joins, unlisted units and unmapped coordinates were excluded fail-closed.

The completed join resolved 245 of 259 source taxa and classified 27,066 of 46,276 observations as native. Native-only refits retained orientation–BIO1 and aspect ratio–BIO4 under the locked same-sign plus BH rule. Chroma–BIO12 and aspect ratio–BIO12 retained their signs but failed native-only BH correction and were withdrawn from the current headline. These outcomes may not be rescued by altering the mask or multiplicity family.

## 7. Claim and figure consequences

- A failed or sign-reversed control removes the corresponding row from the current headline in `final_claims.json`; its frozen primary estimate remains visible as provenance and is not rescued by changing the model family.
- Main Figure 6 is demoted to the Supporting Information because 54 of 216 atlas taxa are direct backbone tips and the retained Pagel-λ fits estimate λ = 0.
- The vacated display position is reserved for bias-control or repeat-photograph evidence only after that evidence is complete.
- Independent detector boxes and continuous-trait reference measurements remain submission-blocking regardless of these computational results.

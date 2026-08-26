# Supporting Information

## Appendix S1. Supplementary methods, tables and figures

### S1.1 Cohorts and operational taxonomic units

The balanced atlas, exhaustive post-detection stream, coordinate and positional-accuracy subsets, primary taxon × 0.25° cohort, and high-resolution subset were retained as distinct analytical cohorts (Table S1.1). Source-assigned names were treated as operational taxonomic units because public observations do not constitute a resolved genus-wide taxonomy. The primary analyses preserve those identifiers; simultaneous collapse of World Checklist of Vascular Plants (WCVP) synonym conflicts provides a sensitivity analysis (Section S1.8; Table S1.10).

Geographic density and filtering are shown in Figure S1.6. Trait assessability varies among regions because public images differ in viewpoint, resolution, lighting and occlusion (Figure S1.7). An unassessable image is measurement missingness, not biological absence.

### S1.2 Detection and trait measurement

A single-class You Only Look Once version 11 nano (YOLO11n) detector localized visible capitula at a confidence threshold of 0.25. Detection supplied regions of interest only. Deterministic image functions then measured orientation relative to Exchangeable Image File Format (EXIF)-oriented image vertical, Commission Internationale de l'Éclairage L*a*b* (CIELAB) corolla lightness and chroma, sine and cosine components of circular hue, and four two-dimensional outline statistics. The definitions, assessability rules and interpretation limits for all endpoints are listed in Table S1.2.

Detector evaluation was separated from detector development by excluding every development photograph and observation before sampling 1,000 audit images; 250 were assigned for double annotation. Human adjudication is incomplete, so precision and recall are not reported. Horizontal mirroring assesses technical repeatability but does not establish biological accuracy. Independent reference measurements are likewise needed for orientation, colour and outline.

### S1.3 Nested visible-variance decomposition

Visible image sums of squares were partitioned among source-assigned taxa, among photographs within taxa and among heads within photographs. The exact components are given in Table S1.4. Two sensitivities retained one deterministic head per photograph or repeatedly sampled 10 photographs per eligible taxon. Both preserved substantial below-taxon variation (Figure S1.1). These components describe visible image phenotypes and are not estimates of genetic variance.

### S1.4 Within-taxon environmental models

Primary ordinary least-squares models used within-taxon centred and standardized traits and CHELSA annual mean temperature, temperature seasonality, annual precipitation and precipitation seasonality. Taxon-clustered standard errors were used, with Benjamini–Hochberg (BH) correction across 36 endpoint-component tests. Supported rows are listed in Table S1.3 and the complete coefficient pattern is shown in Figure S1.2.

Grouped stochastic partial differential equation models fitted by integrated nested Laplace approximation (SPDE-INLA) compared climate, climate plus topography, climate plus soil and full-environment predictor sets. Each endpoint used the same complete cases across groups, a taxon random intercept and a Matérn spatial field. Stable supported patterns are reported in Table S1.5, and within-endpoint differences in the widely applicable information criterion (WAIC) are shown in Figure S1.3.

The high-resolution layer used three two-dimensional involucral contour proxies: projection roughness, outward-positive radial-projection fraction and maximum spine-like projection relative to head radius. The original ≤10-km BIO4 rows were refitted with log minimum head dimension and log1p sharpness and stratified at 150–199, 200–299 and ≥300 pixels. None met the predeclared retention rule (Table S1.6). These proxies do not distinguish longer bract tips from spreading or recurvature and are not direct measurements of defence.

All 36 primary models were also refitted after taxon-specific sine/cosine day-of-year residualization. The ten most represented taxa were omitted separately and the two most represented taxa jointly, with omission membership fixed from observation counts alone. All four non-circular frozen rows retained BH support in the full raw-calendar cyclic-season model and retained sign across omissions, although omission magnitudes varied, especially for chroma–BIO12 (Table S1.11).

A separately locked hemisphere-aware sensitivity then shifted Southern Hemisphere dates by one half-cycle or fitted separate taxon-specific cyclic curves by hemisphere. The cohort contained 2,356 Southern Hemisphere observations, and only *Cirsium vulgare*, *C. arvense* and *C. palustre* were sampled in both hemispheres. All four non-circular rows retained their frozen sign, BH support and omission-sign stability under both definitions (Table S1.13). Day-of-year does not classify developmental stage. A metadata-only preflight identified 20,073 strict-cohort observations with at least two public photographs (Table S1.12), but repeat-photo traits have not been remeasured.

The native-range sensitivity used a separately frozen contract, the WCVP checklist dataset release identified by DOI 10.15468/6h8ucr and a pinned TDWG level-3 geometry commit. Taxon matches with more than one accepted WCVP key and coordinates without an unambiguous level-3 mapping were excluded. Of 46,276 observations, 27,066 were classified native, 10,554 introduced, 5,491 unresolved at the taxon join, 2,100 geographically unmapped and 1,065 unlisted for the accepted taxon. Native-only models retained orientation–BIO1 and aspect ratio–BIO4 under the frozen same-sign plus BH rule. Chroma–BIO12 and aspect ratio–BIO12 retained their primary signs but failed native-only BH correction and were withdrawn from the current headline (Table S1.14).

### S1.5 Among-taxon environmental and phylogenetic analyses

For taxa represented by at least five observations, lower and upper trait quartiles were compared in three-axis environmental principal-component space using centroid distance and Gaussian Bhattacharyya overlap. Each trait was permuted 10,000 times while environmental coordinates, availability and trait distributions were held fixed. Results for all complete taxa and the direct-dated-backbone subset are listed in Table S1.7; Figure S1.8 shows the observed separation geometry.

Taxon-level climate associations were fitted across 50 randomized within-genus placements using Pagel-λ phylogenetic generalized least squares. Table S1.8 and Figure S1.9 report the rows supported in every randomized tree. Only 54 of 216 atlas taxa were direct dated-backbone tips, retained Pagel-λ estimates were zero and direct-backbone non-circular signal was unsupported; these results therefore quantify sensitivity to tested placements rather than resolve thistle phylogeny.

### S1.6 Spatial robustness

For each primary endpoint, the selected SPDE specification was refitted with the same centring, predictor set, mesh and priors to obtain residuals. Global Moran's I used eight nearest neighbours on unit-sphere coordinates and 999 permutations. Residual autocorrelation results are listed in Table S1.9 and shown in Figure S1.4. Taxon trait ranks were also recomputed after omitting each broad geographic region; minimum rank correlations are shown in Figure S1.5.

### S1.7 Multiplicity and interpretation

False-discovery-rate families remained separate by cohort and analysis tier. Circular hue components were interpreted jointly rather than as independent colour states. Cross-sectional trait–environment associations were not interpreted as plasticity, local adaptation, temporal response, pollinator-mediated selection, evolutionary rate or adaptive radiation.

### S1.8 Taxonomic sensitivity

Eight WCVP synonym candidates would merge an active source unit with another active source unit. All eight were collapsed simultaneously before recomputing nested variance and 36 primary climate models. Candidate mappings and outcomes are listed in Table S1.10. This sensitivity tests taxonomic uncertainty while leaving source observation identities transparent.

**Figure S1.1. Nested visible-variance sensitivities in the global thistle atlas.** For nine primary endpoints, below-taxon visible variance in the full decomposition is compared with one-head-per-photograph and equal-replication sensitivities. Substantial below-taxon variation remains under both alternatives.

**Figure S1.2. Primary within-taxon climate coefficients for global thistle photographs.** Standardized coefficients are shown for nine endpoint components × four CHELSA predictors. Asterisks identify the eight rows with BH q < 0.05. Hue sine and cosine are components of one circular colour response.

**Figure S1.3. Grouped spatial-model comparison for global thistle photographs.** Differences in WAIC within each endpoint compare climate, climate plus topography, climate plus soil and full-environment SPDE-INLA models; zero denotes the best-fitting group. The figure compares model fit, not coefficient significance.

**Figure S1.4. Residual spatial autocorrelation after spatial modelling of global thistle photographs.** Moran's I is shown for residuals from the selected SPDE specification for each endpoint, using eight nearest neighbours and 999 permutations. No endpoint had permutation P < 0.05.

**Figure S1.5. Taxon-rank stability under broad-region omission in the global thistle data.** Each value is the minimum Spearman correlation between full-data taxon rankings and rankings recomputed after omitting one geographic region. All endpoint minima remain high.

**Figure S1.6. Sampling geography across filtering stages for global *Cirsium* and allied Cardueae photographs.** One-degree densities are shown on a common log-count scale for (a) detector-positive observations, (b) coordinate-usable observations, (c) observations with positional accuracy ≤10 km and (d) the primary taxon × 0.25° cohort. Maps use a Mollweide equal-area projection and include a 5,000-km equatorial scale bar. They diagnose sampling, not global representativeness.

**Figure S1.7. Geographic variation in trait assessability across global thistle photographs.** Two-degree cells with at least 20 coordinate-usable observations show the usable fraction for (a) orientation, (b) complete visible-colour components, (c) complete gross-outline components and (d) the mean across all nine primary endpoints. Maps use a Mollweide equal-area projection and include a 5,000-km equatorial scale bar. Unassessable images are treated as missing.

**Figure S1.8. Environmental sorting of low- and high-trait thistle taxa.** Observed environmental centroid distance is plotted against Gaussian Bhattacharyya overlap in three-axis environmental principal-component space. Higher distance and lower overlap indicate stronger separation. Point size marks endpoints with at least one all-taxa metric passing BH-corrected permutation tests; rings mark support in the direct-backbone sensitivity. The figure does not imply within-taxon response or causal adaptation.

**Figure S1.9. Historical-placement sensitivity of global thistle trait–climate associations.** (a) Frozen primary within-taxon and randomized-tree among-taxon climate associations are juxtaposed without combining inferential families. (b) Six standardized PGLS coefficients passed BH correction across all 50 randomized trees. This is a placement sensitivity, not a resolved phylogenetic correction: only 54 of 216 atlas taxa are direct dated-backbone tips, retained Pagel-λ estimates are zero and no non-circular direct-backbone trait has FDR-supported phylogenetic signal.

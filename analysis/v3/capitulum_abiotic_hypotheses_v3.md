# Chapter 1 v3 — capitulum × abiotic-environment hypotheses

Status: biological hypothesis redesign before any v3 ecological fitting. This is a retrospective redesign after v2 results were seen, not preregistration. The purpose is to define which abiotic processes are biologically defensible for each measured capitulum module before choosing the final environmental representation.

## 1. General rule

Do not start from a preferred list of BIOCLIM variables and then ask which traits correlate with them. Start from the reproductive structure and ask which abiotic exposures can plausibly act on the measured image phenotype.

The capitulum is treated as a reproductive exposure unit with partially distinct modules:

- orientation / presentation;
- visible corolla colour;
- gross head shape;
- involucre / fine architecture.

Not every module is expected to respond to every environmental variable. Trait-specific mechanistic hypotheses are primary; a whole-capitulum multivariate analysis asks whether those responses form one exposure syndrome or remain multidimensional.

## 2. Environmental time scale: flowering-period exposure before annual summaries

The observation ledger retains exact dates for nearly all source observations. Therefore the primary abiotic representation should be aligned to the photographed reproductive period rather than defaulting to annual BIO variables.

For each observation with a usable date and location, derive long-term monthly climatological exposure for the observation month from CHELSA v2.1 monthly/climatology products where available:

- precipitation (`pr`);
- mean / maximum temperature (`tas`, `tasmax`);
- shortwave radiation (`rsds`);
- vapour-pressure deficit (`vpd`);
- near-surface wind (`sfcWind`).

This is a climatological exposure aligned to the observation month, not actual weather on the photography date. It is also not proof of flowering stage. A bounded alternative may use a centred multi-month window to reduce date/stage noise, but this choice must be fixed without looking at trait coefficients.

Annual BIOCLIM variables, BIO18, growing-season precipitation and NPP are not automatically privileged. They can be compared as broader climatic representations only after the hypothesis-linked monthly exposure matrix is constructed and its coverage/redundancy is quantified.

## 3. H1 — hydrometeor shielding and capitulum orientation [strongest directional hypothesis]

### Biological premise

A more downward-facing capitulum can reduce direct wetting of exposed reproductive tissues. Experimental work in nodding Asteraceae shows that water exposure can reduce pollen viability and that nodding presentation can protect reproduction from wet conditions. This provides a direct functional precedent for an abiotic rain-exposure hypothesis.

### Primary phenotype

- continuous image-referenced orientation angle;
- only after the v3 orientation-stability rule is frozen from image-only perturbation data.

### Primary environmental exposure

- precipitation climatology aligned to the observation month.

### Directional prediction

Higher flowering-period precipitation is associated with a larger image-referenced orientation angle (more downward presentation).

### Secondary representations

- centred multi-month precipitation window;
- annual precipitation / BIO12;
- BIO18 or growing-season precipitation only as broader seasonal alternatives, not competing variables selected by significance.

### Important exclusions

- image vertical is not gravity;
- precipitation association does not prove rain interception, pollen protection or adaptation;
- east/west thermal orientation literature cannot be mapped directly onto this vertical-angle metric, so temperature/radiation are not co-equal directional orientation hypotheses.

## 4. H2 — visible floral colour and abiotic stress environment [strong multivariate hypothesis, weak single-coordinate direction]

### Biological premise

Floral pigmentation and visible colour can covary with solar radiation, water availability and temperature. Anthocyanin-related responses to radiation and drought provide a mechanistic precedent, but JPEG-derived colour is not pigment concentration or spectral reflectance.

### Primary phenotype

Treat colour jointly before interpreting any coordinate:

- lightness;
- chroma;
- circular hue;
- closed colour-composition fractions where available.

### Primary environmental exposures

- shortwave radiation aligned to observation month;
- atmospheric drying (`vpd`) aligned to observation month;
- precipitation / water availability aligned to observation month;
- temperature / maximum temperature aligned to observation month.

### Prediction

The joint visible-colour state changes along radiation / water-stress / thermal exposure gradients.

A fixed directional prediction for chroma alone is intentionally not primary because low chroma can represent both pale and dark colours and the mapping from camera-space chroma to pigment concentration is not generally monotonic.

### Interpretation ceiling

A supported colour-environment association is an optical-phenotype geography result. It is not evidence for anthocyanin concentration, UV pattern, pollinator perception or adaptive pigmentation.

## 5. H3 — protective enclosure / architecture and atmospheric water stress [moderate, measurement-limited hypothesis]

### Biological premise

The Asteraceae involucre forms a protective envelope; bracts and associated structures can contribute to protection from desiccation, radiation and other threats. In a nodding Asteraceae system, elongated involucral bracts substantially reduced UV-B exposure to pollen. These precedents motivate an exposure-buffering hypothesis for head architecture.

### Candidate phenotypes

Use only endpoints whose image construct can defensibly represent enclosure / compactness / architecture after measurement QC. Fine projection and surface proxies must not be renamed as spine length, stiffness, gland density or secretion.

Potential module-level information includes:

- head compactness / gross shape;
- involucre length:width / taper geometry;
- validated cover/display geometry if available;
- projection architecture only as image geometry.

### Primary environmental exposures

- VPD / climatic water demand;
- precipitation / moisture availability;
- shortwave radiation;
- maximum temperature.

### Prediction

Protective-architecture phenotype space is structured along desiccation / radiation exposure gradients.

No universal signed endpoint prediction is frozen until the image metrics are botanically calibrated enough to define what an increase in each metric means biologically.

## 6. H4 — mechanical exposure and head presentation / architecture [secondary hypothesis]

Wind can alter mechanical loading and boundary-layer conditions around reproductive structures, but current image traits do not directly measure stem stiffness, peduncle mechanics or three-dimensional drag area.

Therefore near-surface wind is retained as a secondary abiotic exposure for orientation / architecture modules, not a mandatory primary driver and not a free-standing story if only one proxy endpoint is significant.

## 7. H5 — thermal exposure is a whole-capitulum context, not automatically an orientation hypothesis

Recent Cardueae work shows that living capitula can maintain temperatures below ambient air under hot conditions, indicating that floral thermal balance is biologically important. Other Asteraceae experiments show that orientation can affect flower temperature, but many such results concern azimuthal east/west orientation rather than the vertical nodding angle measured here.

Therefore:

- flowering-period temperature / maximum temperature belongs in the abiotic environment matrix;
- thermal association with colour and whole-capitulum phenotype is biologically motivated;
- a specific vertical-orientation × temperature prediction is secondary unless a stronger mechanistic bridge is established.

## 8. Variables that should not be primary by default

### NPP

Potential NPP is a broad resource/productivity context rather than a direct physical exposure to the capitulum. Retain only if the literature question explicitly concerns resource limitation or if it is needed as a broad habitat-context covariate. It should not have equal mechanistic status with precipitation, radiation, VPD, temperature or wind.

### BIO18 vs growing-season precipitation

Do not choose between these by which gives a better trait result. Once observation-month precipitation is available, both become broader seasonal representations. Compare their redundancy and biological correspondence to photographed reproductive timing before deciding whether either adds information.

### Annual BIO1/BIO4/BIO12/BIO15

These remain useful broad climate descriptors and comparability variables, but should not be called the a priori core merely because they were used in v2.

## 9. Environment-selection rule after hypothesis definition

Build a phenotype-blind candidate environment matrix on the intended native-range ecological cohort. Then, before looking at trait coefficients:

1. report coverage and missingness;
2. standardize units;
3. calculate pairwise correlation and nonlinear dependence diagnostics where useful;
4. calculate matrix rank, condition number and VIF;
5. cluster highly redundant environmental representations;
6. within each biological exposure family, choose the smallest representation that preserves the process meaning;
7. record the choice and rejected redundant alternatives;
8. only then join traits and fit ecological models.

The environmental variables are therefore equal at entry, but not forced to remain equal when they are redundant or biologically distal.

## 10. Proposed hierarchy of ecological hypotheses

### Tier 1 — primary, direct literature bridge

- H1: orientation ↔ flowering-period precipitation / wetting exposure;
- H2: joint visible colour ↔ flowering-period radiation / water stress / thermal exposure.

### Tier 2 — primary only if measurement construct is sufficiently defensible

- H3: head/involucre architecture ↔ desiccation / radiation exposure.

### Tier 3 — contextual / secondary

- H4: wind ↔ presentation / architecture;
- H5: thermal environment ↔ whole-capitulum phenotype, without assuming a vertical-orientation mechanism.

## 11. Whole-capitulum test

After the trait-specific hypotheses are defined, test whether the environment-associated responses collapse into a common exposure-buffering syndrome or remain module-specific. The expected outcome is not fixed. Endpoint-specific results are decompositions of a supported module/environment relationship, not a fishing screen across all endpoint × raster combinations.

## 12. Claim boundary

These hypotheses concern spatial alignment between image-defined phenotype and climatological abiotic exposure. They do not establish plasticity, adaptation, selection, mechanism, pigment chemistry, true three-dimensional orientation, or direct protection from rain / radiation / drought.

Key literature anchors for this redesign include work on nodding Asteraceae and rain/UV protection, floral thermal ecology and orientation, global flower-colour macroecology, and Asteraceae capitulum/involucre functional morphology. Full bibliographic metadata should be synchronized with the submission bibliography after the hypothesis design is finalized.

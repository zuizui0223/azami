# Introduction literature review

## Purpose

This file records the literature base used to draft `manuscript/01_introduction.md`. It separates established findings from claims that still require source-level verification. The review is organized by the logical work each literature group performs in the Introduction rather than by chronology.

## Search scope

### Concept groups

1. intraspecific trait variation and limits of species means;
2. global plant-trait databases and within-species coverage;
3. citizen-science biodiversity images and image-based phenotyping;
4. measurement error, assessability and observation bias in uncontrolled images;
5. Cardueae and Carduus–Cirsium systematics, hybridization and polyploidy;
6. ecological interpretation of capitulum orientation, colour and form.

### Search strings used or scheduled

- `intraspecific trait variation species means community ecology review`
- `global meta-analysis intraspecific trait variation plant communities`
- `TRY database intraspecific variation species mean`
- `BIEN plant trait database multiple records within species`
- `citizen science photographs plant traits image phenomics`
- `iNaturalist continuous plant traits image analysis`
- `image assessability measurement error biodiversity photographs`
- `Carduus Cirsium phylogeny hybridization polyploidy`
- `Cirsium cytotype polyploid species differentiation`
- `Cardueae Hyb-Seq nuclear plastid incongruence`
- `capitulum orientation rain pollinator floral temperature`
- `anthocyanin flower colour climate soil pollinator`

### Inclusion priorities

- peer-reviewed syntheses for broad conceptual claims;
- database description papers for TRY and BIEN;
- primary methodological papers for image recognition and image-derived traits;
- primary phylogenetic, cytogenetic and systematic studies for Cardueae and Cirsium;
- papers that distinguish within-species from among-species inference;
- papers documenting limitations or bias, not only successful applications.

### Exclusions

- studies using images only for species identification when they make no claim about traits, except as background for the field's current emphasis;
- unverified claims that Cirsium underwent adaptive radiation;
- studies of unrelated taxa used to assert a Cirsium-specific mechanism;
- web summaries as final manuscript citations.

## 1. Intraspecific trait variation and species means

### Core references

- **Bolnick et al. 2011, Trends in Ecology & Evolution — Why intraspecific trait variation matters in community ecology.** Foundational synthesis arguing that ecological models based on homogeneous species can miss interaction and community consequences of individual variation.
- **Albert et al. 2011, Perspectives in Plant Ecology, Evolution and Systematics — When and how should intraspecific variability be considered in trait-based plant ecology?** Plant-focused framework for when ITV changes trait-based inference.
- **Violle et al. 2012, Trends in Ecology & Evolution — The return of the variance: intraspecific variability in community ecology.** Conceptual argument for retaining trait distributions rather than means alone.
- **Siefert et al. 2015, Ecology Letters — A global meta-analysis of the relative extent of intraspecific trait variation in plant communities.** Quantifies the contribution of ITV across plant communities and provides direct support for scale-aware trait ecology.

### Claims supported

- ITV can alter ecological inference and is not necessarily negligible relative to interspecific differences.
- Species means are useful summaries but erase distributional information.
- The importance of ITV depends on trait, spatial scale and ecological question.

### Claims not supported without additional evidence

- that all image-derived within-species dispersion is biological;
- that ITV necessarily represents plasticity or local adaptation;
- that one universal threshold determines when species means are invalid.

## 2. Global trait databases

### Core references

- **Kattge et al. 2011, Global Change Biology — TRY: a global database of plant traits.** Establishes the original global trait-data integration framework.
- **Kattge et al. 2020, Global Change Biology — TRY plant trait database: enhanced coverage and open access.** Documents expanded coverage and database architecture.
- **Maitner et al. 2018, Methods in Ecology and Evolution — The BIEN R package.** Describes access to integrated plant occurrence, taxonomy and trait data.

### Balanced interpretation

The Introduction should not portray TRY or BIEN as deficient. They were built to integrate heterogeneous observations and enable broad synthesis. The gap is that downstream analyses often reduce retained records to species-level summaries, while geographic replication and environmental slopes remain sparse or uneven for many species. Public images provide complementary replication, not a replacement for standardized trait measurements.

### Review tasks still open

- identify studies quantifying geographic and taxonomic gaps in TRY;
- identify papers comparing species means with population-level trait records;
- verify how TRY currently exposes multiple observations, location metadata and measurement context;
- distinguish database structure from common downstream use.

## 3. Citizen-science images and image phenomics

### Core references

- **Wäldchen and Mäder 2018, Methods in Ecology and Evolution — Machine learning for image-based species identification.** Reviews automated identification and clarifies the field's dominant taxonomic-classification objective.
- **Van Horn et al. 2018 — The iNaturalist Species Classification and Detection Dataset.** Demonstrates the scale and heterogeneity of citizen-science image data for computer vision.
- **Rzanny et al. 2019, Plant Methods — Flowers, leaves or both? How to obtain suitable images for automated plant identification.** Shows that organ choice and imaging perspective materially influence model performance.
- **Katal et al. 2022, Frontiers in Plant Science — Deep learning in plant phenological research: a systematic literature review.** Documents the rapid growth of image-based phenology and its methodological limitations.

### Current frontier

The literature is rich in identification and phenological-state classification but thinner in globally replicated, continuous, observation-level plant traits linked to local environmental gradients. The strongest novelty claim is therefore not that photographs have never been used to measure traits, but that this study combines:

- organ detection;
- continuous multi-module trait extraction;
- trait-specific assessability;
- observation-level distributions;
- within-species environmental slopes;
- among-species niche contrasts;
- spatial and phylogenetic sensitivity.

### Bias literature required

Add primary studies or reviews on:

- photographer and accessibility bias;
- taxonomic and geographic unevenness in iNaturalist/GBIF;
- illumination and colour-calibration limitations;
- viewpoint dependence in morphology;
- false negatives and detector recall;
- uncertainty propagation from model output to ecological inference.

## 4. Why Cirsium and the Carduus–Cirsium group?

### Core references

- **Ackerfield et al. 2020, TAXON — A prickly puzzle: Generic delimitations in the Carduus–Cirsium group.** Demonstrates conflict between traditional generic boundaries and molecular phylogenetic structure.
- **Herrando-Moraira et al., Hyb-Seq phylogeny of Cardueae.** Provides nuclear and plastid evidence and a temporal framework at tribe/subtribe scales. Exact bibliographic details must be verified before final citation.
- **Bureš et al. 2023, Preslia — Evolution of genome size and GC content in Carduinae.** Documents cytogenomic diversity, including polyploidy and strong phylogenetic structure at broader taxonomic scales.
- **Moreyra et al. 2023, Plants — African mountain thistles: three new genera in the Carduus–Cirsium group.** Shows that mountain lineages formerly placed in broad traditional genera represent distinct clades.
- regional Cirsium studies documenting cytotype differences, natural hybrids and geographically restricted lineages; each should be added with primary-source details.

### Safe claims

- the group is species-rich and taxonomically difficult;
- molecular data challenge some morphology-based generic or sectional circumscriptions;
- hybridization and polyploid taxa occur within Cirsium and close relatives;
- nuclear, plastid, morphological and taxonomic histories may not coincide;
- these features make a single mean plus a single grafted tree an incomplete representation.

### Claims to avoid

- Cirsium is proven to be an adaptive radiation;
- weak phylogenetic signal in this study was caused by polyploidy or hybridization;
- all Cirsium lineages diversified rapidly;
- chloroplast phylogeny alone resolves species history.

## 5. Ecological interpretation of capitulum modules

### Orientation

Literature is needed on floral orientation in relation to rain protection, floral temperature, solar radiation, pollinator access and mechanical support. The Introduction should state these as plausible interfaces, not as mechanisms demonstrated here.

### Colour

Review anthocyanin physiology, temperature and light responses, edaphic effects, pollinator colour vision and repeated white/pink/purple transitions. The Taiwan–Ryukyu white-flower system is a future mechanistic extension, not evidence for the present global analysis.

### Form

Review developmental and biomechanical determinants of capitulum shape, but acknowledge that two-dimensional image outline is viewpoint dependent. Shape is valuable precisely because it can reveal high visible dispersion with weak climate tracking, but biological and technical components must remain separated.

## 6. Gap established by the combined literature

Existing research has established, separately, that:

1. ITV matters;
2. trait databases enable global comparison;
3. citizen-science images support taxonomic and phenological inference;
4. image measurements are sensitive to perspective and quality;
5. thistle evolutionary history is taxonomically and cytogenetically complex.

What remains uncommon is a single global analysis that retains observation-level floral-trait distributions and explicitly separates:

- assessability;
- visible within-species variation;
- within-species environmental responsiveness;
- among-species environmental sorting;
- uncertain historical placement.

This combined gap is the Introduction's central justification.

## Verification status

Before journal submission, every citation in `01_introduction.md` must be checked against the full paper for:

- exact author list and year;
- title, journal, volume and pages/article number;
- DOI;
- whether the cited claim appears in the source;
- whether a newer synthesis changes the interpretation.

The bibliography should label records as `verified_full_text`, `verified_abstract`, or `candidate` during drafting.
# Introduction literature review

## Purpose

This file records the literature base used to draft `manuscript/01_introduction.md`. It separates established findings from claims that still require source-level verification. The review is organized by the logical work each literature group performs in the Introduction rather than by chronology.

## Search scope and procedure

The review targeted six concept groups:

1. intraspecific trait variation and limits of species means;
2. global plant-trait databases and within-species coverage;
3. citizen-science biodiversity images and image-based phenotyping;
4. measurement error, assessability and observation bias in uncontrolled images;
5. Cardueae and Carduus–Cirsium systematics, hybridization and polyploidy;
6. ecological interpretation of capitulum orientation, colour and form.

Searches combined exact article titles, topical terms and taxon-specific terms. Priority was given to peer-reviewed syntheses for broad claims, database description papers for TRY and BIEN, primary methodological studies for image analysis, and primary systematic, cytogenetic or hybridization studies for Cardueae and Cirsium. Web summaries were used only to locate publications and are not manuscript references.

Excluded claims include: a demonstrated adaptive radiation of Cirsium; causal attribution of the present trait patterns to polyploidy or hybridization; treating all image-derived dispersion as biological; and treating missing image measurements as trait absence.

## 1. Intraspecific trait variation and ecological scale

The conceptual base is now strong. Bolnick et al. (2011) establish that individual-level variation can alter species interactions and community processes. Albert et al. (2011) provide a plant-specific framework for deciding when ITV must be incorporated. Messier et al. (2010) explicitly organize trait variance across hierarchical ecological scales. Violle et al. (2012) argue for retaining trait distributions rather than only means, and Siefert et al. (2015) quantify the variable contribution of ITV across plant communities.

Together these papers support three statements used in the Introduction:

- species means are useful but incomplete summaries;
- ITV is often non-negligible and varies among traits and scales;
- phenotypic dispersion should not automatically be interpreted as plasticity or adaptation.

They do not justify interpreting the present 82–99% visible within-species fractions as purely biological variance. The manuscript therefore uses `visible variation` and explicitly includes viewpoint and measurement components.

## 2. Global trait databases

TRY and BIEN are presented as enabling infrastructures rather than as deficient resources. Kattge et al. (2011, 2020) document TRY's integration of heterogeneous plant-trait measurements and its expanded coverage. Maitner et al. (2018) documents programmatic access to BIEN's integrated taxonomy, occurrence and trait information.

The defensible gap is not that these databases contain only means. It is that conventional field trait records are sparse and uneven for many species and traits, and downstream global analyses commonly aggregate them to species-level summaries. The present study is positioned as a complement that provides dense replicated measurements of visible reproductive traits, not as a replacement for calibrated functional-trait measurements.

## 3. Citizen-science images and image phenomics

Wäldchen and Mäder (2018), Van Horn et al. (2018), Rzanny et al. (2019) and Katal et al. (2022) establish the development of image-based identification, detection and phenological classification. They also show why uncontrolled photographs are difficult: classes are long-tailed, images vary in quality and context, and organ choice and viewpoint affect performance.

A targeted novelty search identified Sharma et al. (2025), PlantTraitNet, as a close methodological predecessor that infers continuous plant traits from citizen-science photographs using multimodal learning and uncertainty. It is currently recorded as a preprint and its publication status must be checked immediately before submission. Its existence means the manuscript must not claim to be the first study to infer plant traits from citizen-science images.

The narrower, stronger novelty claim is the integration of:

- detection of a repeated reproductive observation unit;
- continuous orientation, colour and outline measurements;
- trait-specific assessability rather than forced prediction;
- observation-level within-species distributions;
- species-specific environmental slopes;
- among-species environmental-niche contrasts;
- spatial and historical-placement sensitivity.

No located predecessor combined all of these levels in one global floral-trait study.

## 4. Observation and image bias

Troudet et al. (2017) establishes broad taxonomic and societal bias in biodiversity data. Di Cecco et al. (2021) shows that iNaturalist data production is structured by participant behaviour and contribution patterns. These papers support the claim that public observations are not probability samples of organisms or environments.

The manuscript adds image-specific concerns that follow directly from the measurement problem: illumination and camera processing affect colour; distance and occlusion affect detectability; viewpoint affects orientation and outline. These are framed as measurement limitations and motivations for assessability/QC, not as effects quantified by the cited platform-bias studies.

Further literature on calibrated colour extraction and formal propagation of detector false negatives would strengthen Methods and Discussion, but is not needed to establish the Introduction's central gap.

## 5. Why the Carduus–Cirsium group?

Ackerfield et al. (2020) provides the main evidence that traditional generic delimitations do not map cleanly onto molecular phylogenetic structure. Moreyra et al. (2023) provides a well-resolved example in which African mountain thistles previously placed in broad traditional genera represent distinct clades and required new genera.

Bureš et al. (2023) documents genome-size diversity, dysploidy and polyploidy across Carduinae. Bureš et al. (2018) is a species-level example in which ploidy, morphology and genetic evidence supported separation of lineages formerly treated together. Michálková et al. (2023) provides named primary evidence that contemporary hybridization can affect a rare Cirsium lineage.

These sources support describing the group as taxonomically and cytogenetically complex and as containing reticulate histories. They do not establish that Cirsium as a whole is a rapid adaptive radiation. The Introduction therefore treats the group as a stringent test for species-mean and single-tree simplifications, not as a demonstrated radiation.

Publisher-level metadata for Ackerfield et al. (2020), Bureš et al. (2018, 2023) and Michálková et al. (2023) remains to be completed before submission. Moreyra et al. (2023) is fully identified by DOI.

## 6. Capitulum modules

The Introduction currently uses conservative language: orientation, colour and form have different plausible ecological interfaces. It does not state that the study demonstrates rain protection, pollinator selection, thermal adaptation or pigment physiology.

A dedicated search remains necessary for the Discussion, where specific results will be interpreted:

- floral orientation experiments involving precipitation, solar radiation, flower temperature and visitation;
- anthocyanin regulation by temperature, light and water stress;
- edaphic influences on visible floral colour;
- pollinator colour vision and preference;
- developmental and biomechanical determinants of capitulum shape.

These mechanistic papers should not be allowed to retroactively turn the observational global associations into causal claims.

## Gap established by the combined literature

Prior work has established separately that ITV matters, global trait databases enable synthesis, citizen-science photographs support identification and increasingly trait inference, public observations are non-random, and thistle evolutionary history is taxonomically and cytogenetically complex.

What remains uncommon is a global analysis that retains observation-level floral-trait distributions and explicitly separates:

1. measurement assessability;
2. visible within-species variation;
3. within-species environmental responsiveness;
4. among-species environmental sorting;
5. uncertain historical placement.

This combined gap—not a claim that photographs or image-derived traits are unprecedented—is the Introduction's central justification.

## Reference management

- Machine-readable bibliography: `manuscript/references_introduction.bib`
- Human-readable list: `manuscript/REFERENCES_INTRODUCTION.md`
- Claim-to-reference audit: `manuscript/literature/introduction_reference_matrix.csv`

Before submission, every citation must be checked against the publisher record and, for claim-critical sources, the full text. Drafting status is tracked as `verified_bibliographic`, `partial_verification`, or `candidate_current` in the reference matrix.

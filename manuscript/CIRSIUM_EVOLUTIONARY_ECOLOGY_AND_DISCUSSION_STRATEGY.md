# Cirsium evolutionary ecology and Chapter 1 discussion strategy

## Purpose

This document connects the executed Chapter 1 results to the evolutionary ecology of *Cirsium* without converting observational climate associations into claims of adaptation. It also specifies how to retain the concrete trait–environment results while making the paper more general than a conventional species-mean trait–climate study.

## 1. Biological context of the genus

*Cirsium* is a useful system precisely because morphology, taxonomy and genealogy are not expected to align neatly.

### 1.1 A species-rich and taxonomically difficult lineage

The Carduus–Cirsium group is species rich, and molecular studies have repeatedly shown that traditional generic and sectional boundaries are difficult to reconcile with a single clean bifurcating history. Ackerfield et al. described this as a “prickly puzzle” of generic delimitation. Nuclear and plastid analyses of Cardueae, including later Hyb-Seq work, likewise show that historical inference depends strongly on marker type, taxon sampling and treatment of reticulation.

This makes *Cirsium* biologically suitable for asking whether visible traits are conserved, repeatedly assembled or environmentally sorted. It also makes it unsuitable for claiming that one grafted megaphylogeny represents the true species history.

### 1.2 Hybridization is part of the evolutionary background

Numerous natural interspecific hybrids are documented in Eurasian *Cirsium*. Some species hybridize with several sympatric congeners, and recent work has shown that hybridization can threaten narrow endemics through genetic swamping. The invasive tetraploid *C. vulgare* has also been proposed to have an intergeneric hybrid origin.

Hybridization provides a plausible route by which combinations of capitulum traits can move across taxonomic boundaries or appear repeatedly in different lineages. Chapter 1 does not test introgression, but weak direct-tip phylogenetic signal and strong visible heterogeneity are compatible with a history that is more reticulate than a simple conserved-trait model.

### 1.3 Polyploidy changes the meaning of species-level comparisons

Chromosome-number and genome-size studies show diploid and polyploid members within the broader Carduinae and within species complexes traditionally treated as single taxa. The separation of diploid *C. greimleri* from tetraploid *C. waldsteinii* illustrates that ploidy can reveal evolutionary units obscured by morphology alone. Polyploidy can alter reproductive isolation, effective population size, gene flow and the amount of standing variation.

Therefore, broad visible variation in Chapter 1 should not automatically be interpreted as plasticity. It may combine population differentiation, cryptic taxonomic structure, hybrid ancestry, polyploid history, ontogeny, environment and image-related variance.

### 1.4 Rapid radiation is a hypothesis, not an established result of this dataset

The high regional endemism, repeated descriptions of narrowly distributed species, island and mountain differentiation, hybridization and incomplete phylogenetic resolution are consistent with relatively recent or rapid diversification in parts of the genus. However, Chapter 1 does not estimate diversification rates or demonstrate an adaptive radiation.

Use the wording:

> The genus exhibits the hallmarks of a rapidly and reticulately diversifying lineage in several regions, but the timing, rate and adaptive basis of diversification remain unresolved.

Do not write:

> *Cirsium* underwent a proven recent adaptive radiation.

## 2. What the executed results add to this biological context

### 2.1 Visible trait diversity is not organized by one phylogenetic or ecological axis

Across nine primary endpoints, most observed variance occurs within species rather than among species. Species-level PCA is also multidimensional: PC1 explains only 32.9%, and the first three components explain 69.3%. Together with the failure of non-circular endpoints to retain FDR-supported phylogenetic signal in the 54 direct-backbone taxa, this argues against a simple model in which closely related species inherit one conserved capitulum syndrome.

The defensible interpretation is not “phylogeny does not matter.” Instead:

> The current data do not support a single strongly conserved axis of capitulum morphology, and apparent historical signal is sensitive to how unresolved taxa are placed.

In a genus affected by hybridization, polyploidy and uncertain taxonomic boundaries, that negative result is biologically informative.

### 2.2 Variation, environmental tracking and species sorting are different processes

The central empirical result is that species-level within-variation and environmental responsiveness are negatively associated (Spearman rho = -0.333) and occupy all four quadrants. This separates three processes that species-mean trait studies commonly merge:

1. visible heterogeneity within a species;
2. orderly trait change along environmental gradients within species;
3. replacement of species with different mean phenotypes across environments.

The strict within-species analysis finds no BH-FDR-supported primary effects, whereas the expanded pooled table contains four small effects. At the same time, environmental niche contrasts between extreme species quartiles are appreciable for orientation and colour. Thus climate structure is clearer at the species-sorting scale than as a universal within-species response.

This pattern is especially plausible in a lineage where geographic isolation, hybridization, polyploidy and taxonomic turnover can generate sharp differences among populations or species without requiring strong contemporary slopes within every species.

### 2.3 Trait modules retain concrete ecological meaning

The conceptual story must not erase the original trait–environment results.

#### Orientation

Orientation shows one of the largest environmental niche centroid separations between low- and high-trait species (2.398). In grouped spatial models, climate-only models are preferred for orientation, and repeated environmental effects are stronger than for most shape endpoints.

Possible biological interpretation:

- orientation changes exposure to precipitation, solar radiation and convective cooling;
- it changes the visual presentation and accessibility of the capitulum to visitors;
- stem and peduncle mechanics may impose developmental constraints;
- species replacement along climate gradients may therefore be more important than uniform within-species adjustment.

Do not claim rain protection or pollinator adaptation without field tests. Present these as mechanisms that make the observed environmental sorting biologically credible.

#### Colour

Colour has the lowest median total variation among the three modules but the highest median environmental responsiveness. Hue and chroma also show relatively strong environmental niche separation, and grouped spatial effects are concentrated in colour. Soil improves the preferred model for lightness and chroma relative to climate-only models.

Possible biological interpretation:

- anthocyanin production is sensitive to temperature, radiation, nutrient status and stress physiology;
- floral colour also mediates pollinator perception;
- colour may therefore show modest but orderly environmental tracking despite lower total visible dispersion;
- repeated white, pink and purple states across *Cirsium* may reflect multiple developmental or evolutionary routes rather than one conserved transition.

The Chapter 1 data cannot distinguish plastic pigment regulation, local genetic differentiation, pollinator-mediated selection or shared ancestry. These become explicit hypotheses for field, pigment and transcriptomic work.

#### Shape

Shape has the highest median visible within-species variation but weaker environmental responsiveness, high environmental niche overlap for aspect ratio, circularity and solidity, and no BH-supported grouped spatial effects across the four shape endpoints.

Possible interpretation:

- much shape variation may occur at developmental, individual or viewpoint scales;
- the four CHELSA gradients may not capture the relevant selective environment;
- mechanical support, flowering stage, herbivory or local biotic interactions may matter more than broad climate;
- high visible variation is therefore not evidence of high environmental lability.

This is the clearest empirical demonstration of why dispersion and responsiveness must be separated.

## 3. Why this is stronger than a conventional trait–environment big-data paper

A conventional study would assign one trait value to each species, join those values to climate, fit a phylogenetic regression and report significant coefficients. Chapter 1 improves on that design in six ways.

### 3.1 It retains individual observations

Public photographs provide repeated, georeferenced observations rather than one literature-derived species value. This permits estimation of within-species distributions and environmental slopes.

### 3.2 It distinguishes missingness from biological absence

Each trait has explicit assessability and QC. A head that cannot be measured because of viewpoint, occlusion, colour cast or resolution is not coded as a biological state.

### 3.3 It separates biological scales

Among-species climatic structure, strict within-species association, expanded pooled association, spatial-model results and species-level lability are reported as different inferential layers rather than interchangeable tests.

### 3.4 It treats morphology as modular and multidimensional

Orientation, colour and outline are measured continuously and compared as modules. Circular hue is retained as sine/cosine coordinates rather than forced into an invalid linear angle model.

### 3.5 It propagates historical uncertainty

The study reports direct-backbone coverage and random-grafting sensitivity rather than presenting a poorly resolved megatree as exact.

### 3.6 It is reproducible at image scale

The workflow records detector output, crops, continuous measurements, QC decisions, spatial filtering, climate extraction, software versions, workflow runs and checksums. This is a substantive methodological contribution because processing thousands of uncontrolled photographs is not equivalent to downloading a trait table.

## 4. Relation to general plant trait databases

Global plant trait databases have transformed ecology, but their common analytical unit is a species-level value or a set of measurements with highly uneven replication. Global syntheses have shown that intraspecific trait variation contributes substantially to community-level trait variation, yet many comparative analyses still collapse it into a mean because replicated measurements are sparse.

Chapter 1 uses citizen-science imagery to address that sampling bottleneck for visible reproductive traits. Its strongest general contribution is not that public photographs replace standardized measurements. It is that they can reveal when the species-mean approximation is inadequate and identify which taxa, traits and environments deserve targeted standardized sampling.

Recommended wording:

> Citizen-science images complement rather than replace conventional trait databases: they trade experimental standardization for geographic and within-species replication, while explicit assessability models make that trade-off visible.

## 5. How the large photograph set was processed

The Methods and Discussion should make the computational design understandable, not bury it as implementation detail.

1. Public observations were treated as provenance-bearing biological records, not as an unlabeled image scrape.
2. A detector identified visible capitula and created both tight measurement crops and context-preserving crops.
3. Continuous colour, orientation and outline measurements were extracted deterministically.
4. Trait-specific QC rejected unsuitable views without treating rejection as trait absence.
5. Horizontal-mirror technical replicates were used to identify unstable image-derived measurements.
6. Observation-level values were preserved before species aggregation.
7. Spatial thinning and coordinate-quality rules reduced observer and geographic pseudoreplication.
8. Climate extraction, coefficient tables and all sensitivity outputs were versioned and checksummed.
9. The final analysis retained failed and missing measurements in provenance tables rather than silently filtering them from the record.

This workflow turns heterogeneous photographs into an auditable phenomic dataset. The novelty is the combination of scale, continuous multi-trait measurement, within-species replication, explicit missingness and multi-scale inference.

## 6. Recommended final manuscript story

### Opening problem

Trait-based macroecology usually treats species as fixed points in trait space. This is especially limiting in diverse, hybridizing and polyploid lineages, where species boundaries and visible phenotypes need not map neatly onto a bifurcating tree.

### Study opportunity

Public biodiversity photographs contain replicated visible phenotypes across geography, but their heterogeneity requires detection, trait-specific QC and explicit treatment of assessability.

### Main result sequence

1. Most visible capitulum-trait variation lies within species, so species means discard substantial structure.
2. Within-species dispersion and environmental responsiveness are not equivalent and occupy a two-dimensional lability space.
3. Colour, orientation and shape occupy different regions of that space.
4. Concrete climate structure is strongest for colour and orientation and is clearer in among-species environmental sorting than in strict universal within-species slopes.
5. Apparent phylogenetic conservatism weakens under direct-tip restriction, consistent with the difficulty of imposing one resolved history on a hybridizing and polyploid lineage.

### General conclusion

> Species means, single lability rankings and single-tree explanations are inadequate for global reproductive-trait ecology. Public-image phenomics can recover the missing within-species and modular structure, provided measurement assessability, scale and historical uncertainty are modelled explicitly.

## 7. Discussion section architecture

### Discussion 1 — Species means conceal the dominant visible scale of variation

Lead with the 82–99% within-species fractions, immediately state the image-related caveat, and explain why distributions are more informative than species means.

### Discussion 2 — Lability is at least two-dimensional

Interpret rho = -0.333 and the four occupied quadrants. Use shape as the clearest example of high dispersion without high climate tracking.

### Discussion 3 — Specific ecological signals remain visible

Discuss orientation and colour separately, retaining niche divergence, grouped spatial-model preferences and the small pooled effects. Avoid replacing concrete results with purely abstract language.

### Discussion 4 — Why *Cirsium* may generate weak phylogenetic correspondence

Introduce hybridization, polyploidy, taxonomic complexity and regional differentiation. Present these as plausible background processes, not mechanisms proven by Chapter 1.

### Discussion 5 — Public-image phenomics as a complement to trait databases

Explain scale, processing, QC, provenance and what becomes observable only because thousands of images were retained at observation level.

### Discussion 6 — Boundaries and next tests

Prioritize detector validation, standardized colour photography, common-garden or reciprocal-transplant work, pollinator observation, pigment chemistry, nuclear phylogenomics and explicit ploidy/hybrid sampling.

## 8. Literature anchors for the manuscript bibliography

The exact bibliographic metadata and DOI should be checked before submission.

- Ackerfield J. et al. 2020. A prickly puzzle: Generic delimitations in the Carduus–Cirsium group. *Taxon*.
- Susanna A. et al. The Cardueae revisited: nuclear ITS and plastid trnL–trnF/matK evidence.
- Herrando-Moraira S. et al. Nuclear and plastid Hyb-Seq phylogeny of Cardueae and its temporal framework.
- Bureš P. et al. 2018. *Cirsium greimleri*: a new species endemic to the Eastern Alps and Dinarides. *Preslia*.
- Bureš P. et al. 2023. Evolution of genome size and GC content in Carduinae: dysploidy, polyploidy and phylogenetic structure. *Preslia*.
- Bureš P. et al. 2024. Intergeneric hybrid origin of the invasive tetraploid *Cirsium vulgare*. *Plant Biology*.
- Siefert A. et al. 2015. A global meta-analysis of the relative extent of intraspecific trait variation in plant communities. *Ecology Letters*.
- Kattge J. et al. 2020. TRY plant trait database—enhanced coverage and open access. *Global Change Biology*.
- Van Horn G. et al. 2018. The iNaturalist species classification and detection dataset. CVPR.

## 9. Claims to avoid

- *Cirsium* is demonstrated here to be an adaptive radiation.
- Polyploidy caused the observed capitulum-trait patterns.
- Hybridization caused weak phylogenetic signal.
- Colour differences are pollinator adaptations.
- Nodding heads evolved to protect pollen from rain.
- Public photographs measure true biological variance without error.

These are testable next hypotheses, not Chapter 1 conclusions.

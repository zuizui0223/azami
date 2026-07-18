# Methods

## Study design and analytical separation

We used public biodiversity photographs to quantify continuous capitulum traits across *Cirsium* and related thistles while retaining repeated observations within species. The analysis was designed to distinguish three structures that are often conflated in comparative trait studies: visible within-species variation, within-species environmental responsiveness and among-species environmental sorting. These structures were estimated in separate analysis layers and were not combined into a single inference about plasticity, adaptation or evolutionary rate.

The submission-facing workflow is defined by the canonical Chapter 1 pipeline and frozen claim registry. Any change to accepted taxa, trait definitions, measurement thresholds, environmental data or historical-placement scenarios requires a new analysis version and full validation rather than an undocumented manuscript-layer edit.

## Public-image observations and taxonomic scope

Source observations were obtained from public biodiversity records with image, taxon and spatial provenance retained throughout processing. The continuous image-analysis dataset contained 3,725 public observations, 6,626 detected capitula and 216 accepted image-analysis taxa. The strict within-species layer contained 46,276 spatially thinned observations from 259 taxa. The primary two-axis lability analysis used the 102 taxa meeting all predeclared completeness criteria.

Taxon names were retained exactly as they occurred in frozen analysis inputs and are linked to a dated, human-reviewed accepted-name decision table. The taxonomic audit requires one decision per source name, an external authority record, decision date and rationale. Synonym collapses are recorded explicitly. Any decision that merges or removes active analysis units triggers a new analysis version.

## Capitulum detection and context-preserving crops

Capitula were detected from source photographs using the frozen production detector. Each retained detection preserved the source observation identifier, image provenance, bounding box and contextual crop. Tight crops were used for measurements requiring the capitulum boundary, whereas context crops retained stem and display information needed for orientation assessment.

Detector credibility is evaluated independently using manually annotated source images and one-to-one intersection-over-union matching. Precision, recall, F1, false positives, false negatives, true no-capitulum images and human-unassessable images are reported separately. These metrics and the final licence-safe demonstration panels remain submission gates and are not inferred from the analysis dataset itself.

## Continuous trait measurements

Nine primary endpoints represented three capitulum modules: orientation, visible corolla colour and two-dimensional outline. Measurements were generated deterministically from the frozen production functions rather than by manual scoring.

Orientation was quantified from the capitulum–stem geometry retained in the context crop. Visible colour was summarized from the retained corolla mask using continuous colour coordinates, with hue represented jointly by sine and cosine components to respect circularity. Outline measurements were calculated from the segmented capitulum boundary and convex hull and included aspect-ratio and compactness-related endpoints.

Trait-specific assessability rules were applied before analysis. A failed or unassessable measurement was retained as missing data rather than converted to biological absence or a zero value. The measurement audit displays assessable and unassessable examples, tight and contextual crops, production-generated orientation lines, colour masks, outlines and convex hulls, together with horizontal-mirror technical checks and low-, middle- and high-value examples.

## Quality control and observation-level retention

All trait values were retained at observation level before aggregation. Quality-control rules were specific to each trait module because orientation, colour and outline fail for different photographic reasons. Coordinates, taxon names, source identifiers, image licences, detector outputs and measurement status were propagated into downstream tables.

The continuous comparison layer and strict within-species layer are distinct cohorts. Counts, false-discovery-rate summaries and sensitivity results are therefore always reported with their cohort and analysis scope. Missing endpoints were not interpreted as taxon-level absence.

## Environmental data and spatial cohorts

Environmental predictors were extracted for georeferenced observations from the frozen environmental data sources recorded in the submission manifest. The strict primary within-species cohort was restricted to observations with coordinate uncertainty no greater than 10 km. Spatial thinning was applied before species-specific environmental modelling to reduce clustering of repeated public observations.

The accepted analysis includes separate coordinate-quality and pooled summaries. The strict primary cohort and the expanded pooled coefficient table are not combined into a single significance count. Environmental associations are interpreted as observational climate tracking or environmental responsiveness, not as demonstrated plasticity or local adaptation.

## Partitioning visible variation and species-level trait structure

For each primary endpoint, visible variance was partitioned into within- and among-species components using the frozen continuous-analysis workflow. The within-species fraction describes dispersion among uncontrolled photographs assigned to the same species and can include biological heterogeneity, developmental state, viewpoint and measurement components.

Species-level trait structure was summarized with principal-component analysis on the accepted continuous endpoints. This analysis was used to assess whether orientation, colour and outline could be reduced to one dominant species-level axis.

## Within-species variation and environmental responsiveness

Species-specific within-variation values were calculated for eligible trait–species combinations from the strict spatial observation table. Environmental responsiveness was calculated from standardized species-specific trait–environment coefficients in the exhaustive within-species climate analysis.

The primary species-level lability cohort required at least 10 observations per eligible trait–species combination, at least six measured traits per species and at least three environmental predictors per trait. Species-level indices were built only for the 102 taxa satisfying all completeness criteria. The within-variation and environmental-responsiveness indices were compared using Spearman rank correlation and divided at their cohort medians to describe four high/low combinations. Module summaries used medians rather than means.

Sensitivity analyses repeated the species rankings with minimum sample-size thresholds of five and 20 observations. Rank correlations with the primary threshold quantified whether the central two-axis result depended on the selected minimum sample size.

## Among-species environmental sorting

Among-species environmental structure was analysed separately from within-species responsiveness. Grouped spatial models evaluated orientation, colour and outline endpoints against climate and soil predictor groups while accounting for spatial structure. Model comparison and globally adjusted fixed-effect summaries were used to identify which trait modules showed the clearest environmental sorting.

A complementary trait-extreme analysis compared environmental niches of species at opposite ends of each species-level trait distribution. Environmental niche centroid separation and overlap were interpreted as among-species sorting and not as direct evidence of selection or adaptive divergence.

## Residual spatial and broad-region robustness

The accepted models are not replaced by the diagnostic audit. Observation-level fitted values or residuals are joined back to frozen observations through explicit observation and endpoint identifiers. Residual spatial structure is quantified with endpoint-specific k-nearest-neighbour Moran's I and permutation p-values using a recorded seed.

A manually reviewed broad-region lookup is required; regions are not inferred automatically from coordinates. Regional observation and taxon coverage are summarized, and leave-one-region-out analyses quantify the stability of taxon rankings after each broad region is omitted. Final diagnostic values and any resulting limitation language remain pending until the frozen prediction/residual export and reviewed region lookup are supplied.

## Historical-placement sensitivity and molecular-data audit

Historical sensitivity was evaluated across a dated backbone and alternative deterministic and randomized placements for taxa absent from that backbone. Direct-backbone and grafted-tree results were retained separately. Phylogenetic-signal and historical regression outputs were classified according to whether support survived restriction to directly represented tips and alternative placements.

Hue sine and cosine components were not interpreted as independent biological tests; circular colour inference requires a joint treatment. The molecular-database audit summarized available nucleotide, ITS, plastid, plastome and sequence-read records and was used to evaluate whether a uniformly sampled resolved nuclear species tree could be assembled for the complete analysis set.

## Multiplicity, sensitivity and claim control

False-discovery-rate correction was applied within the predefined reporting families of the frozen workflows. Strict and expanded analyses were not pooled after correction. Effect size, cohort definition and model scope were considered together rather than using adjusted significance alone.

All manuscript-facing numerical claims are validated against `manuscript/final_claims.json`. The manuscript uses the terms visible within-species variation, environmental responsiveness, climate tracking, among-species environmental sorting, module dependence and historical-placement sensitivity. It does not treat the observational analyses as demonstrations of plasticity, local adaptation, selection, rain protection, pollinator causation, evolutionary rate or a resolved *Cirsium* species tree.

## Reproducibility and pending submission gates

The canonical scripts, inputs and outputs are mapped in `analysis/ch1/pipeline.json` and executed through versioned GitHub Actions workflows. The durable submission bundle will include final CSV tables, PNG and PDF figures, software environments, workflow run identifiers, commit SHA, artifact digests, SHA-256 checksums and a file-level output map.

Before this section is frozen for submission, three externally supplied audit products must be inserted without altering the accepted analyses: independent detector and measurement-audit outputs, the completed taxonomic decision table, and the residual spatial and broad-region diagnostic results. Public photographs will not be redistributed beyond their licence terms; the archive will preserve source identifiers and licence metadata.

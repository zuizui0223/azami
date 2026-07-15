# Phylogenetic-signal feasibility for Chapter 1

Status: decision document for the submission analysis

## Question

Can phylogenetic signal be estimated for the two lability axes and the orientation,
colour and shape modules without treating a plastid genealogy or a heavily grafted
megaphylogeny as the true species history of *Cirsium*?

## Current in-repository coverage

The image-derived species table contains 216 accepted Chapter 1 taxa. In the dated
GBOTB/LCVP backbone audit, only 54 taxa are direct backbone tips; 162 taxa are
inserted within *Cirsium* by a grafting rule. Direct dated-tip coverage is therefore
25.0%. The existing historical sensitivity already includes deterministic S1 and
S3 placements and 50 randomized S2 trees.

This is enough for a **tree-uncertainty sensitivity analysis**, but it is not enough
to present one estimated value of Blomberg's K or Pagel's lambda as a property of a
resolved *Cirsium* species tree.

## Sequence and tree resources that exist

### Broad taxonomic backbones

- Open Tree of Life can resolve names and return an induced synthetic topology.
- GBOTB.extended.LCVP provides dated vascular-plant branch lengths but sparse direct
  coverage within *Cirsium*.
- These sources are suitable for auditing and sensitivity, not for resolving recent
  reticulate relationships within the genus.

### Public molecular data

- GenBank/ENA/DDBJ contain many *Cirsium* accessions, including nuclear ribosomal ITS
  and several plastid loci; complete plastomes are available for a subset of taxa.
- SRA contains genomic and transcriptomic projects for a smaller and taxonomically
  uneven subset.
- Cardueae-wide studies have combined nuclear ITS with plastid trnL-trnF/matK, and
  more recent tribe-level work uses Hyb-Seq nuclear and plastid data.
- No currently identified public study supplies a densely sampled nuclear species
  tree matching the 216 Chapter 1 taxa.

Database counts change continuously. Exact accession and taxon coverage must be
frozen by an executable audit at submission rather than copied manually into the
manuscript. The audit should record query date, query string, accession list,
accepted-name mapping and the number of Chapter 1 taxa with each data type.

## Why a chloroplast-only tree is not adequate

A plastid tree represents one cytoplasmic lineage. In *Cirsium*, hybridization,
chloroplast capture and geographic introgression can make plastid history differ
from the nuclear species history. Polyploid and especially allopolyploid taxa may
have more than one parental history, which cannot be represented fully by a single
strictly bifurcating tree.

Consequently:

- do not construct the headline tree from whole plastomes alone;
- do not interpret plastid clades as the unique species relationships;
- do not use one chloroplast topology to claim trait conservatism or repeated
  evolution;
- treat nuclear/plastid disagreement as biological uncertainty, not a nuisance to
  be removed by concatenation.

## Polyploidy and hybrid evidence

Ploidy is relevant both biologically and statistically, but available records are
heterogeneous. Chromosome-count databases, floras and individual cytological papers
must be mapped to accepted names with vouchers and source citations. Absence of a
record is `unknown`, not diploid.

The Chapter 1 policy remains:

- retain all taxa in the main descriptive analysis;
- create reviewed evidence fields for ploidy and hybrid origin;
- run documented-hybrid/polyploid exclusions only as sensitivity analyses;
- do not delete taxa based on automated text inference or missing evidence.

## Recommended phylogenetic-signal analysis

### Traits to test

1. species within-variation index;
2. species environmental-responsiveness index;
3. module-level variation and responsiveness for orientation, colour and shape;
4. individual trait summaries only in the supplement.

### Tree sets

Run the same signal calculation on three nested sets:

1. **Direct-tip analysis** — only taxa that are direct dated-backbone tips. This is
   the most defensible but least powerful test.
2. **All-taxon tree ensemble** — deterministic S1/S3 plus 50 randomized S2 trees.
   Report the distribution, not one preferred estimate.
3. **Evidence-filtered sensitivity** — repeat after excluding only taxa with
   documented hybrid or polyploid complications, while leaving unknown taxa intact.

### Statistics

For every continuous endpoint and tree:

- estimate Pagel's lambda by maximum likelihood with a likelihood-ratio comparison
  to lambda = 0;
- calculate Blomberg's K with tip-label permutation;
- record sample size, tree identifier, direct/grafted status composition and model
  convergence;
- summarize median, 2.5–97.5% range, sign/support frequency and the fraction of
  trees with nominal support.

K and lambda answer related but different questions. Neither should be translated
as an evolutionary rate. A stable non-zero signal means related taxa resemble one
another more than expected under tip permutation or independence on the tested tree
ensemble; it does not prove that the trait evolved by Brownian motion.

## Decision rule for the manuscript

Promote phylogenetic signal to the main text only when:

- direct-tip and all-taxon ensemble results agree qualitatively;
- the central estimate and uncertainty are stable across randomized grafting trees;
- the result is not driven by one region, one module or documented hybrid/polyploid
  taxa;
- at least approximately 40 direct tips have non-missing values for the endpoint.

Otherwise report it as a historical-constraint sensitivity in the supplement.

## Current feasibility judgement

- **Pagel's lambda across the existing 50-tree ensemble: executable now.**
- **Blomberg's K across the same ensemble: executable now.**
- **A definitive genus-wide nuclear species-tree analysis: not executable from the
  current repository alone.**
- **A plastome-only replacement tree: technically executable but scientifically
  unsuitable as the headline tree.**
- **A new target-capture/Hyb-Seq species tree: high value for a later chapter, but a
  separate data-generation project rather than a requirement for Chapter 1.**

## Required next implementation

1. freeze a machine-readable Chapter 1 taxon list;
2. query and save NCBI nucleotide, organelle and SRA coverage by accepted taxon;
3. join reviewed chromosome-count/ploidy evidence;
4. add an R script that calculates K and lambda over the existing tree ensemble;
5. output direct-tip-only, all-tree and ploidy/hybrid sensitivity tables;
6. add the tree IDs, software versions, seed and commit hash to the submission
   bundle.

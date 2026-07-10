# Historical constraints for Chapter 1

## Why this is a sensitivity analysis

The global continuous-trait dataset spans 216 *Cirsium* taxa, but there is no
single, densely sampled nuclear species tree for all of them. A plastid tree
would follow one cytoplasmic lineage and can disagree with the nuclear history
after hybridization or chloroplast capture. Allopolyploids additionally have
more than one parental history and cannot be represented completely by one
strictly bifurcating tree.

For that reason, Chapter 1 does not use one mega-tree as unquestioned truth.
Historical structure is introduced in layers.

## Layers

1. **Within-species environmental models**
   
   Outcomes and predictors are demeaned within species. These models estimate
   within-species associations and remove all time-invariant between-species
   differences without requiring a tree. They are the least tree-dependent
   historical control.

2. **Open Tree audit**
   
   Open Tree of Life is used to audit taxonomic name resolution and to obtain an
   independent synthetic induced topology. It is not assigned ad hoc divergence
   times for the production PGLS.

3. **Dated GBOTB sensitivity**
   
   `V.PhyloMaker2` is used with `GBOTB.extended.LCVP`. The output records which
   taxa are direct backbone tips and which are grafted within *Cirsium*.

4. **Alternative grafting scenarios**
   
   Deterministic scenarios S1 and S3 are compared. Scenario S2 is repeated on 50
   randomized trees for the nine main continuous endpoints. A coefficient is a
   robust historical-sensitivity candidate only when its direction is stable in
   at least 90% of random trees, the 2.5–97.5% tree range excludes zero, and at
   least 80% of fits have nominal `p < 0.05`.

5. **Ploidy and hybrid evidence**
   
   Every taxon begins as `unknown`. A taxon is excluded as a documented or
   suspected hybrid/polyploid only after a reviewed source, page/table/figure and
   supporting quotation are recorded. Unknown is never converted to diploid or
   non-hybrid.

## Models

For each species-level endpoint, the same four standardized CHELSA predictors
are fitted:

- annual mean temperature;
- temperature seasonality;
- annual precipitation;
- precipitation seasonality.

The deterministic tree analysis reports:

- ordinary least-squares reference estimates;
- Brownian-motion PGLS;
- PGLS with Pagel's lambda estimated by maximum likelihood.

The randomized-tree analysis uses Pagel-lambda PGLS and summarizes coefficient
and lambda distributions across trees.

## Interpretation boundary

A stable PGLS result means that a coefficient is not easily removed by the
particular dated backbone and grafting uncertainties tested here. It does not
show that:

- the GBOTB tree is the true *Cirsium* species tree;
- plastid and nuclear histories agree;
- hybridization or allopolyploid parentage has been resolved;
- the environment caused the observed trait difference.

The strongest claims require agreement among measurement-QC sensitivity,
within-species models, positional-accuracy sensitivity, and the historical-tree
sensitivity described here.

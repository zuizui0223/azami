# V4 native phylogenetic placement sensitivity

Successful numerical and report run:
https://github.com/zuizui0223/azami/actions/runs/39435174409

Code: `c1266ec807f83b60e111a5d9125e7b14abc29786`.
Inputs and tree-bundle hashes: `receipt.json`.
Exact source-name matches and excluded taxa: `taxon_mapping.json`, `coverage.json`.
Full coefficients, variances and boundary flags: `coefficients.csv`.

The frozen five environmental variables and their native-population
within/between decomposition are unchanged. Nine primary coordinates use
endpoint-specific observations. This is observation-level phylogenetic
Gaussian REML, not regression of species-averaged responses.

There are 468 phylogenetic fits (nine coordinates by 52 placements) and nine
independent-species fits on identical covered cohorts. No failures were
recorded. 194 taxa are covered; missing taxa are mostly sparsely observed,
but their exclusion still limits taxonomic scope.

For orientation, 169912 observations are retained. The between-species
rainfall coefficient remains positive on all placements: 0.4722–0.5078
degrees per 100 mm; the envelope of conditional 95% intervals is
[0.2051, 0.7748]. Within-species coefficients are about −0.0455 degrees
per 100 mm and every interval includes zero.

The temperature–chroma direction remains positive between species and
negative within species. Between-species radiation coefficients retain
positive lightness and negative chroma. The radiation–width-profile
direction remains negative between species and positive within species.
These are component-level sensitivities, not new whole-colour/outline tests.

Intervals are model-based, conditional on fitted variance components and
pointwise (not multiplicity-adjusted). They are not interchangeable with
the previous species-cluster intervals. Do not promote newly narrower
within-species intervals as fresh discoveries. The 52 trees are grafting
scenarios, not independent confirmations or resolved molecular species
histories. This analysis does not estimate cross-trait evolutionary
covariance or establish adaptation, plasticity or speciation.

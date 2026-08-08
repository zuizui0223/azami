# Completed SPDE residual and broad-region robustness audit

## Frozen source and diagnostic-only refit

The accepted grouped SPDE-INLA coefficient tables remain unchanged. The diagnostic workflow refitted only the lowest-WAIC frozen specification for each of the nine primary endpoints in order to recover observation-level fitted values and residuals that the compact production workflow had not archived.

Workflow run: `31152400475`  
Artifact: `8983877726`  
Artifact digest: `sha256:151161f926dbc92d35832a0bd71622e4f1f3018c9047a7ba8b4ad94926273dca`

The diagnostic refits reproduced the frozen model support closely. Across the nine endpoints, diagnostic-refit WAIC minus frozen WAIC ranged from -2.56 to +2.50, with identical frozen observation and species counts for every endpoint.

## Residual Moran's I

Global k-nearest-neighbour weights were constructed on 3-D unit-sphere coordinates rather than raw longitude/latitude degrees. The neighbour graph was fixed within endpoint and residual values were permuted 999 times. Up to 10,000 rows per endpoint were sampled deterministically for the permutation diagnostic.

| Endpoint | Moran's I | Permutation P |
|---|---:|---:|
| Corolla hue cosine | 0.0024 | 0.604 |
| Corolla hue sine | 0.0028 | 0.574 |
| Corolla chroma | -0.0096 | 0.059 |
| Corolla lightness | -0.0019 | 0.716 |
| Orientation | -0.0068 | 0.126 |
| Aspect ratio | 0.0023 | 0.658 |
| Circularity | 0.0030 | 0.532 |
| Solidity | -0.0031 | 0.508 |
| Width-profile CV | 0.0013 | 0.774 |

No endpoint retained detectable residual spatial autocorrelation at P < 0.05, and all observed Moran statistics were close to zero.

## Broad-region coverage and leave-one-region-out stability

Broad regions were assigned reproducibly from Natural Earth 1:50m admin-0 polygons. Polygon-contained points received the country continent directly; small-island omissions used an explicitly flagged nearest-country fallback. Two of 46,276 source observations remained unmapped and were retained as a separate diagnostic stratum rather than silently assigned.

Source observation coverage was dominated by Europe and North America, with additional representation from Asia, Oceania, South America and Africa. Because the residual diagnostic contains one row per endpoint when assessable, coverage totals in the diagnostic table repeat observations across endpoints.

For each endpoint, species trait means were recomputed after removing one broad region at a time and compared with the full-data species ranking. Minimum leave-one-region-out Spearman correlations were:

| Endpoint | Minimum rho |
|---|---:|
| Solidity | 0.856 |
| Aspect ratio | 0.859 |
| Hue sine | 0.866 |
| Orientation | 0.874 |
| Hue cosine | 0.877 |
| Circularity | 0.879 |
| Width-profile CV | 0.882 |
| Chroma | 0.936 |
| Lightness | 0.972 |

The lowest values occurred when Europe was removed. No endpoint showed a qualitative collapse of species rankings under broad-region omission.

## Submission interpretation

This audit supports the adequacy of the frozen spatial control for the descriptive and within-species results: the best supported frozen SPDE specifications leave negligible residual global spatial structure, and species rankings are robust to broad-region omission. These are diagnostics, not a new inferential model family, and they do not convert observational associations into causal or adaptive effects.

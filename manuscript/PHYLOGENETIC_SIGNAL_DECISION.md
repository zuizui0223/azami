# Phylogenetic signal decision

## Frozen source

This decision is based on the validated workflow run `29391037056` from commit `6bf05b63e89e6bbf256924839045f7de5cf39740`.

The source analysis contained:

- 216 taxa,
- 54 taxa directly represented in the GBOTB/LCVP backbone,
- 50 randomized grafting trees,
- 12 registered continuous endpoints,
- 636 successful fits and zero failed fits.

## Molecular database coverage

The live NCBI audit found the following discovery-level coverage among the 216 taxa:

| Resource | Taxa with at least one record | Coverage |
|---|---:|---:|
| Any nucleotide record | 138 | 63.9% |
| ITS-related nucleotide record | 133 | 61.6% |
| Plastid/chloroplast record | 129 | 59.7% |
| Complete plastome | 11 | 5.1% |
| SRA record | 120 | 55.6% |

These counts do not establish orthology, voucher correctness, ploidy, or suitability for species-tree inference. Complete plastomes are especially sparse, and plastid-only coverage cannot resolve reticulate or allopolyploid histories.

## Signal robustness rule

A phylogenetic signal is treated as directly supported only when Blomberg's K or Pagel's lambda remains significant after Benjamini-Hochberg correction in the 54 direct-backbone taxa.

Signals detected only in the full grafted tree are classified as grafting-sensitive. Hue sine and cosine are components of one circular trait and are not interpreted separately, even when one component is significant.

## Result

After FDR correction and the direct-tip requirement:

- no non-circular endpoint had direct-backbone-supported signal,
- shape circularity and shape solidity were grafting-sensitive,
- both hue components require a joint circular test,
- the remaining eight endpoints had no FDR-robust signal.

## Manuscript decision

Phylogenetic signal remains a supplementary sensitivity analysis rather than a central result.

The defensible conclusion is not that capitulum traits lack evolutionary history. It is that the present broad vascular-plant backbone does not provide robust evidence for trait conservatism once uncertain Cirsium placements are separated from directly represented tips.

This supports the main paper's emphasis on visible within-species variation and climate tracking rather than on a strongly resolved macroevolutionary history.

## Future route

A stronger evolutionary analysis would require a nuclear multi-locus or target-capture species framework that explicitly evaluates:

- hybridization,
- chloroplast capture,
- auto- and allopolyploidy,
- gene-tree discordance,
- uncertain taxonomic identity and vouchers.

Plastid-only reconstruction should not be used as the primary species tree for this dataset.

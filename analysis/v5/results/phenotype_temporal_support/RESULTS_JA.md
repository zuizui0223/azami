# Temporal support audit: rapid change in novel climate (2026-11-26)

[Run 45771204026](https://github.com/zuizui0223/azami/actions/runs/45771204026) successful; artifact 15093604841, downloaded ZIP SHA256 27e199c1ddf1b5ede26ce6f18462c7c14b9ff4572140a9139ecd9c9551d2c895 matches published digest. Actual outputs: 108 date coverage rows and 1728 repeated support rows. Primary window pair 2016–2020 vs 2021–2025, 0.25 degree cell x month, at least 3 observations, 3 observers, 2 years in each window.

| Species | Range | Climate | Chroma cell-month strata | Local cells | Blocks |
|---|---|---|---:|---:|---:|
| C. arvense | introduced | analogue | 536 | 293 | 9 |
| C. arvense | introduced | novel | 0 | 0 | 0 |
| C. vulgare | introduced | analogue | 576 | 293 | 10 |
| C. vulgare | introduced | novel | 24 | 18 | 3 |
| C. arvense | native | all | 261 | 162 | 9 |
| C. vulgare | native | all | 177 | 129 | 7 |

Neither species nor any of the eight endpoints achieves the frozen gate (>=30 local cells and >=8 blocks) for novel climate in any grid, window or record-count specification. Across all introduced novel groups, max local-cell support is 2/24 for arvense/vulgare, max broad blocks 2/7 respectively (including all_source). Chroma specifically: arvense max 1 cell/1 block, vulgare max 18 cells, max 6 blocks (different specifications). No qualifying early-to-late repeated cells in novel climates for chroma. Large recent source cohorts do not fix this: introduced arvense chroma 385/10136/34381 records in early/middle/late, vulgare 760/13400/36338; the counts are not local repeated measures. No trait changes fitted. Stop the novel-climate temporal analysis, without loosening gates. This is lack of identification, not evidence of no rapid change. Observation dates are not introduction dates. Mixed climate groups are not interchangeable with novel climate.

Next: inventory public genomic resources and documented introduction history. The contemporary photo findings remain descriptive/exploratory; the hypothesis of adaptation from standing genetic variation requires variant origin, linked phenotypes, ancestry control and fitness or credible selection evidence. No use of herbarium dried-petal RGB as equivalent to living flower chroma.
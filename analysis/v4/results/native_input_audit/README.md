# Native input audit

Actions: https://github.com/zuizui0223/azami/actions/runs/39426282205
Code: 87f9587. Five focused tests passed; decoding and aggregate audit completed.
This is an input audit, not a fitted ecological result.

- Native population: 313,350 observations / 306 accepted taxa.
- Eligible orientation: 170,066 / 244 taxa.
- Eligible L* and C*: 186,439 each / 244 taxa; complete colour: 186,390.
- Eligible four-component outline: 176,048 / 243 taxa.
- All three main modules: 157,779 / 242 taxa, before environmental completeness
  and minimum per-species integration support.
- Bract SI: 1,303 / 111 taxa. It is not required for main-module completeness.
- Correctly decoded nine-environment complete observations: 312,001.

Raw temperature–VPD r=0.9083; rainfall–GSP r=0.8599. Those high correlations
persist at between and within scales. Candidate reduction must be justified
before fitting, not by which conditional coefficient becomes significant.

Coordinates are present for all observations. The tested observer-ID and region
field names are absent. This does not prove that observer information is absent
from every upstream archive. No photographer adjustment is currently supported
by this prepared table.

Full overlap counts and three-scale correlation matrices are in `audit.json`.
The decoding receipt tracks the raw source and corrected table checksums.
No private observation identifiers or coordinates are included.

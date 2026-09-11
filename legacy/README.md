# Historical implementations

Start current numerical reproduction at [the current runbook](../reproducibility/CURRENT_ANALYSIS.md), not here.

- `ch1_global/v2/`: preserved image acquisition, detector development, categorical experiments, continuous measurements and numbered endpoint scripts. These are historical implementations, not additional current entry points. Their presence does not authorize new image acquisition or model training.
- `v2/analysis/`: frozen endpoint-level analysis and recovery implementations.
- `v2/figures/`: frozen endpoint-level figure builders, not the current construct-level figures.
- `v2/ORIGINAL_REPRODUCTION.md`: original runbook for immutable revision `584af97b050d15701f26ce1facea212d5b648d4d`. Its commands apply to that checkout, before relocation.
- `tests/`: historical tests whose target implementations had already been removed before this cleanup. Run them on their original revision, not as current tests.
- `workflows/`: original per-stage Actions definitions; current execution is unified in `.github/workflows/reproduce-current-analysis.yml`.

The current construct sensitivity chain explicitly imports two retained numerical implementations from `legacy.v2.analysis`: spatial sensitivity and historical-placement sensitivity. These are shared dependencies, not independent current paper entry points. Keeping their numerical implementation avoids changing the analysis during repository cleanup.

The complete old-to-new path mapping is [legacy_path_map.json](../reproducibility/legacy_path_map.json). Frozen results in `analysis_outputs/`, reference figures in `reproducibility/figures/`, artifact hashes, and immutable historical Git references are preserved. Historical receipts can legitimately contain pre-relocation paths; use the map to locate them in this checkout.

No historical statistical result, missing measurement or failed audit has been removed or promoted by this reorganization.

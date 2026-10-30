# Legacy research versions

The current GEB manuscript (Results 3.1 BIO18–angle, 3.2 colour gradients, 3.3 local climate-analogue resemblance) is described in [root README](../README.md) and the [canonical manuscript entry](../manuscripts/README.md). These directories are preserved for provenance and historical reproduction, **not** promoted as present V4 results.

## V3 archival snapshot

- [`v3/README.md`](v3/README.md) lists the exact immutable pre-migration main revision and the old→archive path map.
- `v3/analysis/v3/` contains all prior V3 biological-construct analytical scripts and their original scientific results.
- `v3/reproducibility/` contains the V3 numerical runner, frozen public figures, release manifests and source receipts (including historical V2 reference illustrations).
- `v3/workflows/active_as_of_20260928/` contains GitHub Actions previously active on main, including V3 publication steps; these are **deliberately disabled** on current main.
- `v3/workflows/previously_archived/` contains the older per-stage V3 workflow definitions.
- `v3/tests/` contains V3-specific tests and fixtures from the former active analysis.

Run historical V3 reproduction from pre-migration commit `b463dbc938e454baa54b2e480f7f1d8badd7795d` rather than assuming moved path imports work unchanged. File-tree moves preserve each original blob SHA.

## Earlier V2 files

- `ch1_global/v2/`: early acquisition, detector development and continuous measurement scripts.
- `v2/analysis/`: preserved frozen endpoint-level numerical methods and recovery.
- `v2/figures/`: original v2 figure generation.
- [`v2/ORIGINAL_REPRODUCTION.md`](v2/ORIGINAL_REPRODUCTION.md): replay instructions pinned to `584af97b050d15701f26ce1facea212d5b648d4d`.
- `tests/`: already archived earlier test scripts unrelated to active V4.
- `v3/reproducibility/legacy_path_map.json`: prior reorganization's path mapping.

No archived result should be used in a V4 manuscript as if it belongs to the 24 October 2026 science freeze.

## Deferred V4 research and deletion policy

- [`v4_s2_precipitation_niche_20261025/`](v4_s2_precipitation_niche_20261025/README.md) preserves an exploratory climate-niche strand that is **not** part of the current GEB claim family.
- [Current V4 source scripts/results](../analysis/v4/README.md) and the [article](../manuscripts/v4_s1_s2_frozen_20261025/GEB_RESEARCH_ARTICLE.md) stay outside `legacy/`; the H3c result was incorporated later than the original S1/S2 freeze.
- [Archive and branch triage](../docs/BRANCH_ARCHIVE_TRIAGE_20261030.md) requires link, PR-base, SHA and workflow checks before any branch or historical-file removal.

**No historical source, workflow or artifact is made authoritative merely because a current navigation page links to it.** Old publication workflows remain disabled; the code license does not override external image/data rights.

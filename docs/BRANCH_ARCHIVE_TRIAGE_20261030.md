# Branch and archive triage — non-destructive inventory (2026-10-30)

The repository contains **137 branch references** in the GitHub REST listing at this date; **49 distinct branch names point to the same archived commit `7d28f2340e1d1b7fed5bd7a5611f52853e09a388`**. These are *candidate* redundant references, not evidence that deletion is safe. No branch or scientific file has been deleted by this cleanup.

## Branch groups and decisions

| Group | Examples | Decision |
| --- | --- | --- |
| **Current GEB candidate** | [PR #159](https://github.com/zuizui0223/azami/pull/159) → `manuscript/v4-lego-conditional-interaction-refinement-20261027` | **Keep active.** Its paper, figures and verification assets are authoritative for the current *working* manuscript, not yet for published science. |
| **Separate exploratory V4.1** | [PR #160](https://github.com/zuizui0223/azami/pull/160) → `analysis/v4-bio18-spatial-holdout-20261029` | **Keep separate.** Do not backport its new claims into the frozen GEB paper without explicit evaluation. Its PR base is the PR #159 branch. |
| **Other open draft analysis branches** | [PR #144](https://github.com/zuizui0223/azami/pull/144), [#145](https://github.com/zuizui0223/azami/pull/145), [#146](https://github.com/zuizui0223/azami/pull/146), [#148](https://github.com/zuizui0223/azami/pull/148) | **Do not delete.** Record research question, frozen scientific outputs, PR base and dependency status, then decide whether to merge, close or preserve as a named snapshot. |
| **Explicit archives and snapshots** | `archive/*`, `snapshot/v4-final-bio18-holm-20261024` | Preserve as references until release evidence/DOIs/checksums and open-PR ancestry are audited. |
| **Many names at one old SHA** | Earlier `agent/*`, `feature/*`, `run/*` and others pointing to `7d28f234...` | **Review candidates only.** Check release tags, links, deployed Actions, PRs and consumers before any removal. |
| **Default `main`** | `main` (old three-part README until PR #159 merge) | Do **not** force-push. Merge or promote only after review of open PR dependencies and CI. |

## Archive ownership

- [`legacy/v3/`](../legacy/v3/README.md): frozen prior Chapter 1/V3 manuscripts, code and evidence.
- [`legacy/v2/`](../legacy/v2/ORIGINAL_REPRODUCTION.md): historical V2; its Zenodo numerical DOI must never be relabelled V4.
- [`legacy/v4_s2_precipitation_niche_20261025/`](../legacy/v4_s2_precipitation_niche_20261025/README.md): excluded S2 rain/occupancy exploration; kept only as research history.
- Current source code, contracts and checked results remain in [`analysis/v4/`](../analysis/v4/README.md); present figures/text remain in [`manuscripts/v4_s1_s2_frozen_20261025/`](../manuscripts/v4_s1_s2_frozen_20261025/).

## Required before any deletion or move

1. Re-list branches and their SHAs; gather PR head/base refs, linked issues, release tags, workflow references and published commit checksums.
2. Verify all affected import paths and markdown links, including archived CI and SI citations.
3. Record an immutable archive reference (commit SHA and source/receipt manifest) that survives branch cleanup.
4. Propose removals in a **separate** reviewable change; get explicit approval for deleting branch refs, source artifacts or entire historical directories.
5. Run frozen numerical, manuscript, navigation and code-contract CI again. Document both what remains and what was removed.

**Implemented here:** a current README, stable directory indexes, source and submission runbook, and a CI-backed navigation check. **Not implemented:** physical relocation or deletion of scientific files, closing PRs, branch pruning, changes to numerical fits, or any public data release.

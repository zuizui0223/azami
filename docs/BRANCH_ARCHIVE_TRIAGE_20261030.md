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

## Verified complete reference snapshot — 30 October 2026

The reproducible, non-secret input is in [`BRANCH_REF_INPUT_20261030.json`](BRANCH_REF_INPUT_20261030.json). It records **137 branches**, **157 PRs** (six open), one tag ref and one draft release from GitHub at PR #159 head `cf053d7f227f`. Head SHAs in this file are historical snapshots, not live branch heads. Regenerate a report with the read-only [audit script](../scripts/repository_ref_audit.py), or replay the snapshot offline:

```bash
python scripts/repository_ref_audit.py --input-json docs/BRANCH_REF_INPUT_20261030.json --output-json /tmp/azami-branch-audit.json
python -m unittest discover -s tests -p 'test_repository_ref_audit.py' -v
```

**Duplicate group:** 49 branch names point to `7d28f2340e1d1b7fed5bd7a5611f52853e09a388`; 42 names appear as heads in historical PR metadata, and seven have **no historical PR-head-name match**. None of those 49 names is the head or base of an open PR. This does **not** prove absence of external links, scripts, deployments or collaborator dependencies.

### Seven no-PR-head-name-match candidates — NOT approved for deletion

- `agent/ecological-discussion-references-20260821`
- `agent/supplement-figures-20260820`
- `cleanup/code-only-prepublication-20260902`
- `codex/ch1-spatial-scale-reframe`
- `feature/ch1-paleoniche-pilot`
- `figures/real-photo-yolo-world-map-v2`
- `tmp-ignore`

All seven point to `7d28f2340e1d1b7fed5bd7a5611f52853e09a388`.

### Active PR dependency chains — preserve every named ref

| Dependency chain | State | Why the branch must remain |
| --- | --- | --- |
| [#160](https://github.com/zuizui0223/azami/pull/160) → [#159](https://github.com/zuizui0223/azami/pull/159) → `main` | #160 exploratory V4.1, #159 current GEB article | #160 is based on #159's working branch; merging/removing either changes the review context. |
| [#148](https://github.com/zuizui0223/azami/pull/148) → [#146](https://github.com/zuizui0223/azami/pull/146) → [#144](https://github.com/zuizui0223/azami/pull/144) → `main` | Three older draft PRs | Their bases form a dependency chain, not independent pruning targets. |
| [#145](https://github.com/zuizui0223/azami/pull/145) → `main` | Older draft PR | Retain until its scope and overlap with the GEB manuscript are reviewed. |

The only tag in this snapshot points to `0476640e740e327bc86d868994d38eba33dd73f9`, not the 49-name shared tip. An existing draft release must not be published or removed during this cleanup. Annotated tags may require dereferencing before comparing commit identities.

**Completed:** archival refs and dependencies are reproducible from a frozen inventory and tested read-only code. **Not performed:** deleting branches or artifacts, closing/merging PRs, publishing releases, altering source data, renaming active import paths, or changing scientific estimates.


## Post-merge status (2026-10-30)

- [PR #159](https://github.com/zuizui0223/azami/pull/159) **merged to `main`** at `e2229a20045044977afe235255355ca3892186fa` using a merge commit, so the original scientific commit ancestry is preserved. Its old working head may still exist as a branch; this alone is not reason to delete it.
- [PR #160](https://github.com/zuizui0223/azami/pull/160) was **retargeted to `main`** without merging its V4.1 science. It remains open and currently reports merge conflicts: 12 commits ahead and 139 behind at the moment of comparison. Reconcile deliberately, not by force-pushing `main`.
- The historical 137-branch/157-PR [snapshot](BRANCH_REF_INPUT_20261030.json) remains unchanged as audit evidence. Its PR #159 state and base are a **dated snapshot**, not live metadata.
- Zenodo release preparations are in [the V4 data-only handoff](ZENODO_V4_RELEASE_HANDOFF.md). No DOI/deposit/publication was created by this merge or branch triage.

No old branch was deleted, and PRs #144/#145/#146/#148 were not closed. If another audit is requested later, refresh the snapshot first.

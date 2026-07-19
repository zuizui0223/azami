# Repository review and cleanup decisions

Review date: 2026-07-18

## Scope

This review examined the submission-facing Python utilities, workflow dependency installation, manuscript control files and the frozen numbered analysis layer. It did not reinterpret results or modify executed Chapter 1 calculations.

## Findings

1. The numbered `ch1_global/v2` scripts are provenance-bearing frozen implementation paths. They are numerous, but deleting or renaming them before the durable archive is built would break workflow and checksum traceability. They remain untouched.
2. New submission audits under `analysis/` had started to duplicate generic dataframe checks such as required-column validation, duplicate-key detection and complete-text validation.
3. GitHub Actions installed overlapping dependencies independently, leaving package versions and import availability implicit.
4. `manuscript/TARGET_JOURNAL_STRATEGY.md` conflicted with the later repository-grounded rule that a journal is not frozen until a dated current-scope matrix is completed. It was not referenced by the canonical manuscript workspace and was removed as obsolete guidance.
5. Audit scripts remain command-line applications rather than being moved into the package. This preserves their stable documented paths while allowing shared non-scientific helpers to be imported.

## Changes in this cleanup

- added `pyproject.toml` as the single Python package and dependency definition;
- added the `azami_ch1` package under `src/`;
- centralized generic table validation in `azami_ch1.tabular`;
- removed duplicate validation implementations from the spatial diagnostic preparation script;
- changed the corresponding workflow to install the repository package rather than an ad hoc dependency;
- removed the obsolete target-journal strategy file.

## Dependency groups

- base: pandas and submission-facing table utilities;
- `spatial`: NumPy and SciPy for spatial diagnostics;
- `image-audit`: OpenCV, Pillow and Matplotlib for detector and Figure 1 audits;
- `test`: numerical dependencies needed by the unit suite.

## Deletion policy

Safe to delete now:

- duplicate generic helper implementations after callers use the package;
- unreferenced manuscript strategy files that conflict with the current control hierarchy;
- generated caches and local outputs if present outside frozen artifacts.

Do not delete before durable release:

- numbered scripts used by frozen workflows;
- frozen output tables, manifests or checksums;
- scripts named in `analysis/ch1/pipeline.json`;
- evidence ledgers and claim-control files;
- old analysis implementations when their exact path is still cited by an artifact or manuscript ledger.

## Next cleanup pass

After the package installation passes CI, migrate the remaining non-scientific helpers in taxonomic and spatial audit scripts, then consolidate repeated workflow provenance writing. Scientific measurement functions should remain in their frozen provenance paths until a new analysis version is deliberately released.

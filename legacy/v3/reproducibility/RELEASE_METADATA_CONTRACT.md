# Chapter 1 Zenodo data-release metadata contract

The public Zenodo object is a **data-only analysis-input archive**. GitHub remains the authority for analysis code, runbooks, frozen reference outputs, figures and replay receipts.

Final mode validates metadata with `reproducibility.release_metadata_contract` before packaging.

## Zenodo content boundary

The archive contains only:

- the exact continuous-measurement input artifact;
- the exact nine-predictor environment input artifact;
- the exact broad-region lookup artifact;
- the exact historical-placement-tree artifact;
- the exact frozen `observation_native_status.csv`;
- a small README, input contract, checksums and release manifest.

It does **not** contain:

- analysis code;
- manuscript or Supporting Information files;
- figure exports;
- fitted/reference outputs;
- replay receipts;
- original photographs;
- detector training data or weights.

The archive manifest pins the GitHub repository and exact code commit required to interpret the data.

## Creator boundary

The existing v2 Zenodo dataset was published with creator `ZHANG, Ruiqi`. The prepared v3 metadata inherits that dataset creator for continuity of the same archive concept. This is a **dataset-deposition creator field**, not a statement about manuscript authorship or manuscript author order. ORCID and affiliation are optional and remain null unless explicitly supplied.

## Before final bundling

The final metadata must settle:

- approval of the data-only release title/description;
- approval of the archive strategy (normally a new version under existing concept DOI `10.5281/zenodo.22295790`);
- approval of the data/third-party reuse wording;
- the exact GitHub code commit associated with the inputs;
- `release_approved=true`.

The published v2 record `10.5281/zenodo.22295791` must remain unchanged. After the new version is published, it must be downloaded without owner credentials, checksums verified, and the GitHub numerical replay run cleanly against the archived inputs.

The metadata contract is not a direct Zenodo REST payload. It is the fail-closed internal approval record used by the bundle builder.

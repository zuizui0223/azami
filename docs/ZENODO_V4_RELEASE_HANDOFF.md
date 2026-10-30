# Zenodo V4 handoff — data-only deposit, not a publication

**State (2026-10-30): PREPARATION ONLY.** [GEB manuscript PR #159](https://github.com/zuizui0223/azami/pull/159) was merged to `main` at commit `e2229a20045044977afe235255355ca3892186fa`. This document does **not** register a DOI, create a Zenodo draft, publish a dataset, upload proprietary images, grant third-party data rights, or submit the paper.

## One research object, separate from the existing V2/V3 history

Create a **new Zenodo Dataset record** for *Cirsium* V4 **aggregate-only results and replay scripts**, not the old V2 data DOI `10.5281/zenodo.22295791`. The existing GitHub draft release tagged `private-v3-numerical-20260908` contains encrypted V3 working sources and is **not** the correct V4 deposit. Do not turn on automatic archival of all GitHub releases as a shortcut: a GitHub–Zenodo integration archives release snapshots, whereas the intended citable object here is a curated, rights-cleared aggregate data package.

The [draft record metadata](../manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/ZENODO_V4_METADATA_WORKING_DRAFT.json) is **deliberately non-depositable**: **CC BY 4.0 is selected for original author-owned derived aggregate data**, but creators, file-level third-party redistribution rights, review access and approval remain unresolved. No false DOI or author identity is inserted.

## Source and candidate identity

The existing [release candidate report](../manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/AGGREGATE_RELEASE_CANDIDATE_README.md) records two original deterministic ZIPs and a new CC BY 4.0 rights-scoped candidate:

| Candidate | Expected size | SHA-256 | Use |
| --- | ---: | --- | --- |
| `cirsium_v4_aggregate_release_candidate.zip` | 2,344,600 bytes | `3e2fa32f210977d3aae5a46a9f18e53f88a4ce7ed0af7e4d8057ca0aafb51414` | Provenance-preserving **internal** verification; original audit identifiers may remain |
| `cirsium_v4_review_aggregate_candidate.zip` | 2,333,208 bytes | `0c18472b67c9478a63c8eb0a926096ecf6501e200f314059dc46de80b2255c38` | Metadata-reduced review candidate; **not certified anonymous** due to required copyright attribution |
| `cirsium_v4_ccby4_scope_release_candidate.zip` | 2,347,009 bytes | `d4ba5af840f9cd65e4c8abc0d48cdc1becf97800a9baca9573f245e51d914780` | **Unpublished** license-scope candidate with 14 scientific resources byte-unchanged; third-party rights pending |

The two original candidates have 63 ZIP entries; the CC BY 4.0 candidate has 65. Its SHA-256 has been verified twice through independent deterministic builds, and offline replay completed all nine commands and five byte-identical numerical CSV outputs. These local candidate ZIPs remain outside Git and are **not public or author-cleared**.

Primary fixed inputs: source artifact `12547876704` (SHA-256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`), H3 robustness artifact `12891286463` (SHA-256 `216fd6b8b722dd615ebbce9f58412ce245477c25ec65057afe9b5962cee33a91`), verified table/figure artifact `12956327726` (SHA-256 `0c6879ebddbf8ba56363049ca09fd69c7260150d0b38bcdcf52fe59f285a5c32`). These are **historical verification receipts**, not new archive DOIs.

## How to rebuild and check

Use the [original deterministic candidate builder](../manuscripts/v4_s1_s2_frozen_20261025/build_review_aggregate_bundle.py) with the **exact** three authenticated artifact ZIPs and a fresh local output directory, as documented in the [aggregate candidate README](../manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/AGGREGATE_RELEASE_CANDIDATE_README.md). Do not silently substitute a newer artifact if one of these downloads has expired. Preserve the original source hashes.

Build the exact license-scope candidate from the verified original ZIP with [the versioned generator](../scripts/build_zenodo_ccby4_overlay.py) and [the rights-scope documentation](../manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/LICENSE_DATA_CC_BY_4_0.md):

```bash
python scripts/build_zenodo_ccby4_overlay.py \
  --source /path/to/cirsium_v4_aggregate_release_candidate.zip \
  --output /path/to/cirsium_v4_ccby4_scope_release_candidate.zip
```

The builder requires the exact original SHA-256 and checks the new expected SHA-256 without modifying the input.

The read-only pre-deposit guard works without a Zenodo account:

```bash
python scripts/validate_zenodo_v4_release.py \
  --metadata manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/ZENODO_V4_METADATA_WORKING_DRAFT.json

# Only after recreating and checking the exact ZIP bytes:
python scripts/validate_zenodo_v4_release.py \
  --metadata manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4/ZENODO_V4_METADATA_WORKING_DRAFT.json \
  --bundle /path/to/cirsium_v4_ccby4_scope_release_candidate.zip

python -m unittest discover -s tests -p 'test_v4_zenodo_release_guard.py' -v
```

Metadata-only validation checks the release state and required approvals; **it cannot inspect unpublished ZIP bytes**. Providing `--bundle` additionally checks the selected ZIP's expected SHA-256, size, member count, unsafe paths, duplicate filenames, symlinks and sensitive internal-provenance filenames. Existing bundle-specific tests and `replay_aggregates.py` provide more extensive structural and numerical checks. None of these checks grants permission to publish.

## Required decisions before creating a real Zenodo record

1. **Creators and funding:** confirm the record's authors, order, ORCID/affiliation where appropriate, and funding declarations. Do not infer them from a GitHub user name. Complete the deposit form and re-check the research article's attribution.
2. **License and redistribution:** author-controlled aggregate results are assigned **CC BY 4.0** and original code remains MIT. Verify rights and citations for any third-party photographs, observational metadata, taxonomy and raster-derived inputs before public distribution. Zenodo supports mixed-license records.
3. **Which files can be seen:** choose a rights-cleared **public dataset candidate** separately from an **anonymous reviewer copy**. Source-provenance and attribution may defeat double-anonymous review. If files are restricted, note that Zenodo's **metadata remains public**.
4. **Exact bytes:** rebuild and check the ZIP(s); read the archive manifests and the complete source privacy audit. Do not upload original image crops, precise source coordinates, uploader names, decryption keys, encrypted source assets, internal identifying crosswalks or the entire GitHub repository as the dataset.
5. **Zenodo publication:** only after explicit author approval, create a *new* Dataset upload, add the approved file, fill and save the record as a draft, and preview. Reserve a DOI only if required for manuscript files. Publishing is a **separate, explicit action**, after which record/DOI changes are constrained. Do not rely on the previous GitHub draft release or historical DOI.

## After actual deposit (not yet performed)

Record the **exact DOI**, dataset version, Zenodo record link, GitHub commit `e2229a2`, final ZIP names/sizes/SHA-256, creator/license decisions and access conditions in the manuscript's Data and Code Availability statement and `docs/`. Ensure any DOI included in the paper is for the **data** object, not a future journal article's DOI.

Official Zenodo documentation: [Create an upload](https://help.zenodo.org/docs/deposit/create-new-upload/), [DOI reservation](https://help.zenodo.org/docs/deposit/describe-records/reserve-doi/), [License choices](https://help.zenodo.org/docs/deposit/describe-records/licenses/), [GitHub integration](https://help.zenodo.org/docs/github/enable-repository/) and [Record lifecycle](https://help.zenodo.org/docs/deposit/about-records/).

**Hard stop:** CC BY 4.0 approval applies only to eligible author-owned derived aggregates. No third-party rights, final creator metadata, upload, DOI declaration, publication, release tagging or auto-archival are implied by this guide.

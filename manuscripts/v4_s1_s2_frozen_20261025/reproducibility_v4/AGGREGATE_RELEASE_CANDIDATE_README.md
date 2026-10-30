# V4 aggregate release candidate — 30 October 2026

## Outcome and supported claim

[`build_review_aggregate_bundle.py`](../build_review_aggregate_bundle.py) creates two local packages from three checksum-pinned aggregate artifacts: a source-provenance package and a review copy. Both verify the frozen aggregate results and regenerate numerical tables/figures without credentials or network access after Python dependencies are installed. **Neither package refits the statistical models from observation-level data or recreates measurements from original photographs.**

No DOI, hosted reviewer link, external publication, release visibility change or historical Zenodo approval is used. A final licence for author-owned aggregate output data remains pending. The code remains under the existing MIT licence.

## Actual candidate files

The verified local outputs are under `../submission/aggregate_rc_20261030/`; that generated directory is intentionally not versioned with these sources:

| File | Bytes | SHA-256 |
|---|---:|---|
| `cirsium_v4_aggregate_release_candidate.zip` | 2,344,600 | `3e2fa32f210977d3aae5a46a9f18e53f88a4ce7ed0af7e4d8057ca0aafb51414` |
| `cirsium_v4_review_aggregate_candidate.zip` | 2,333,208 | `0c18472b67c9478a63c8eb0a926096ecf6501e200f314059dc46de80b2255c38` |

Each package has 63 files, including its own `review_manifest.json`, `BUILD_VALIDATION.json`, exact requirements, data dictionary, legal notices and offline runner. `internal_provenance.json` stays outside the review ZIP and records original source/run/artifact identities, member hashes, copied-file crosswalks and every JSON metadata redaction. The final ZIP checksum is a sidecar rather than a circular hash inside the ZIP.

### Stable reviewer source identifiers

| ID | Resource |
|---|---|
| R001 | Complete 208 primary model tests (80 S1 + 128 S2) |
| R002 | All 28 raw and adjusted photographic covariance pairs |
| R003 | Six original native/introduced colour coefficients |
| R004 | Complete 80-test S1 raw/BH/Holm table |
| R005 | Both original S1 angle–BIO18 scale coefficients |
| R006 | Original exploratory H3c climate-analogue result |
| R007 | H3 climatic and photographic breadth |
| R008 | H3d trait/covariance context |
| R009 | S1 geographic-regressor sensitivity |
| R010 | S1 geographic block bootstrap |
| R011 | S2 six-row geographic block bootstrap |
| R012 | Full-source focal-species support audit |
| R013 | Complete post-positive H3c threshold/trait/region/bandwidth results |
| R014 | Complete post-positive H3c six-seed results |

Exact paths and hashes are in each package's manifest. The first two filenames are `si/Table_S1_full_208_BIO18_tests.csv` and `si/Table_S2_all_28_RV_pairs.csv`. No fictional anonymous access URL is supplied.

## Versioned reproducibility sources

This directory preserves the handoff documentation, [`source_resource_manifest.json`](source_resource_manifest.json), and two complete source aggregate objects in [`source_aggregates/`](source_aggregates/). The resource manifest records the review-copy paths and hashes; its paths are relative to the root of the extracted review package. It is a reviewer resource index, not the identifying internal provenance crosswalk.

The two JSON files in `source_aggregates/` are **exact original artifact-member bytes**, including their source-audit metadata. Their hashes therefore differ from the metadata-neutralized R013 and R014 copies in the review manifest. The accompanying [source preservation receipt](source_aggregates/README.md) records their original hashes and privacy-check scope. This versioned directory is not itself an anonymous review package.

The eight focused behavioral checks are versioned at [`tests/test_v4_review_aggregate_bundle.py`](../../../tests/test_v4_review_aggregate_bundle.py). From the repository root, run:

```bash
python -m unittest discover -s tests -p 'test_v4_review_aggregate_bundle.py' -v
```

### What aggregate replay does and does not establish

Replay checks the integrity of the archived results, reconstructs the five full numerical CSVs from those aggregate objects, and renders the specified numerical figures. It confirms that the supplied tables and displays agree with the archived coefficients, test statistics and result objects. It does **not** independently recompute the fitted coefficients, permutation statistics, bootstrap samples or image measurements from the original observations.

Complete model refitting would additionally require the protected observation-level inputs, their source-matched environmental values and the original fitting/measurement workflows, which are outside this package. Preserving the full robustness and six-seed objects closes the archival gap for those outputs; it does not turn aggregate replay into full source-data reproduction. The manifests and replay receipt retain `observation_level_refitting_supported: false` or the corresponding `observation_level_models_refitted: false` flag.

## Source inputs, all inspected

1. Primary aggregate artifact `12547876704`, source run `40062126006`, SHA-256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`: eight allowlisted members, including the three source-locked JSONs required by the freeze validator. Individual record IDs/coordinates/uploader fields were absent from recursive schema inspection.
2. H3 aggregate artifact `12891286463`, source run `40859868003`, SHA-256 `216fd6b8b722dd615ebbce9f58412ce245477c25ec65057afe9b5962cee33a91`: five complete JSON objects. This preserves the full sensitivity and six-seed data that were previously available only in an expiring Actions artifact. The package retains the separate original H3c frozen summary; the later sensitivity run does not replace it.
3. Verified figure/table artifact `12956327726`, SHA-256 `0c6879ebddbf8ba56363049ca09fd69c7260150d0b38bcdcf52fe59f285a5c32`: explicit five-CSV/22-image/freeze-receipt allowlist. The historical `figure2_cross_functional_trait_pairs` asset is a supplementary cross-domain display. The current geographic main Figure 2 is now source-verified separately; its pooled coverage data, renderer and provenance receipts accompany the manuscript outside these numerical aggregate ZIPs.

The builder rejects unexpected or duplicate ZIP members, path traversal, archive symlinks, oversized members, wrong archive hashes, protected JSON/CSV fields and individual observation/photo links. It never loads the original encrypted/private source archives or a whole repository directory into the payload.

## Review-copy transformation

The canonical source package retains original source bytes. The review copy replaces repository/authorship references and numerical workflow identifiers with neutral audit aliases. It preserves every scientific numeric leaf and list/dictionary structure, and records metadata changes in the internal crosswalk. Copied code has only source-identity and checksum constants adjusted to the corresponding review copies. The complete scientific checks are rerun against those review copies.

`LICENSE` and `NOTICE.md` are byte-identical to the repository originals. Their attribution to Azami contributors is legally required and remains visible. **Full author anonymity is therefore not certified.** An author/journal-approved route for legally preserving attribution while handling double-anonymous review remains unresolved; the builder does not erase copyright notices to make an anonymity check pass.

The candidate includes no raw photographs/crops, observation/uploader/photo identifiers, protected coordinates, climate raster bytes, model weights, encrypted private observations, Git history or credentials. Species-level names and taxonomic authority keys in R012 are aggregate descriptors, not individual records.

## Rebuild

The builder only reads local inputs. Use the exact three archives named above and an environment with NumPy, SciPy, pandas, statsmodels, patsy and Matplotlib:

```bash
python manuscripts/v4_s1_s2_frozen_20261025/build_review_aggregate_bundle.py \
  --source-aggregate-zip /path/to/v4-source-12547876704.zip \
  --robustness-zip /path/to/v4-h3-robustness-12891286463.zip \
  --verified-assets-zip /path/to/geb-v4-verified-assets-90c21b0.zip \
  --output-dir /path/to/new_candidate_directory
```

The output directory must be new, preventing silent overwrites. To use either extracted package, create an isolated Python 3.12 environment, install its `requirements.txt`, and run `python replay_aggregates.py`. The replay is entirely local. It checks manifest hashes, runs nine verification/rendering commands, byte-compares all five exported CSVs, and checks that the numerical figures are produced. Rendering is verified; figure pixels are not asserted identical across environments.

## Verification performed

- Eight focused behavioral tests pass, including a failing-then-passing uploader-name rejection case.
- Both complete candidate self-replays pass: nine commands and five byte-identical CSVs each; row counts 208, 28, 6, 80 and 2.
- Final ZIP inventories match the manifest exactly, all member/resource hashes verify, and the MIT legal notices match the repository byte for byte.
- The original and review versions of every primary aggregate JSON preserve all scientific numeric leaves, with only explicitly classified run/job/artifact identifiers excluded from that comparison.
- The final Figure 2 status correction changed only each packaged README and its file-inventory manifest; all R001–R014 resource paths/hashes and scientific bytes remained unchanged. Both complete self-replays passed again.
- A second full build produced byte-identical ZIPs. All eleven numerical PNGs were inspected through Pillow; their metadata contained only Software and dpi fields and no identifying content.
- An independent code reviewer reran all eight unit tests and found no concrete critical or important issue in the builder or focused tests. That review did not claim an additional independent complete replay.
- The full project `pytest` command was attempted but the runtime does not have `pytest` installed. This is not a claim that the repository-wide test suite passed.

## Remaining publication gates

This closes local aggregate preservation, source-byte inspection, a concrete candidate inventory and aggregate replay. It does not close the data-output licence decision, true anonymous reviewer access, a v4 persistent archive/DOI or complete original-source numerical reproduction. The manuscript's Figure 1 photograph rights were handled separately; none of those photographs or photographer identities are part of this aggregate package. The geographic main Figure 2 also has a separate source-aggregation workflow. WorldClim raster redistribution and item-specific iNaturalist rights remain distinct from the MIT software licence. The older v2 DOI and archived v3 approval do not establish permission or availability for this package.

Generated `submission/` contents are ignored by the repository's existing `.gitignore`. The builder, focused tests and this `reproducibility_v4/` source directory are outside that ignored output directory and can be versioned normally. Binary candidates and `internal_provenance.json` remain generated handoff artifacts; the identifying crosswalk is not included in the review ZIP or this versioned source directory. Publishing a candidate remains a separate concrete decision.

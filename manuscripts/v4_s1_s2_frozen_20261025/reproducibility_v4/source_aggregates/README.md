# Preserved full H3c source aggregates

These two files are exact JSON member bytes from the source-verified H3 aggregate artifact. They retain every diagnostic condition and all six random-sequence results. They are original aggregate outputs, not observation-level data or a substitute for refitting the original models.

| File | Bytes | Original SHA-256 |
|---|---:|---|
| `H3C_ROBUSTNESS_SENSITIVITY_AGGREGATE.json` | 25,078 | `8e3d61dc5ce265abbc021239b4bf338431f9f00fa4112b76be765af716c8a24d` |
| `H3C_MONTE_CARLO_SEED_STABILITY_AGGREGATE.json` | 13,750 | `3c5519882a34782aa28dcdb22e74143fc314bf44be6c730739792d4d5cd55f6b` |

The parent ZIP is artifact `12891286463` from run `40859868003`, SHA-256 `216fd6b8b722dd615ebbce9f58412ce245477c25ec65057afe9b5962cee33a91`. Before preservation, the builder verified the archive checksum, exact five-member inventory, safe paths, absence of duplicate members and symlinks, and size bounds. Its recursive protected-field and individual observation/photo-link checks passed for every JSON member. Each saved file was then compared byte-for-byte with its original member. No values, formatting or metadata were rewritten.

The preserved JSONs contain pooled scientific results and source-audit metadata. They do not contain original records, exact source coordinates, uploader identities, observation/photo identifiers, source photographs, raster bytes or credentials. The privacy check permits aggregate scientific descriptions such as trait coordinates and geographic sector names; it does not promise that the retained audit metadata is anonymous.

The review package has separately neutralized R013 and R014 copies. Its hashes are listed in `../source_resource_manifest.json`; those are intentionally different from the original hashes above. The identifying source/review crosswalk remains in the generated companion `internal_provenance.json`, outside the review ZIP and outside this versioned directory. No resource ID changed when these source bytes were preserved.

## Interpretation and rights

The sensitivity conditions and random sequences are post-positive exploratory checks of the same original source cohort. They do not create independent replications or redefine the original exploratory Holm-2 family. The original H3c result remains separately preserved in the repository's source-verified summary.

Full model refitting is not supported by this directory or the aggregate review package. Original software retains the existing MIT licence and required attribution to Azami contributors. Aggregate-output licensing and any external deposit remain pending; preservation does not create a DOI, grant upstream data rights, or certify anonymous access.

# Chapter 1 release metadata contract

The public-release bundle must not become `release_ready=true` merely because a JSON file was supplied. Final mode validates the metadata with `reproducibility.release_metadata_contract` before packaging.

Start from `zenodo_release_metadata.template.json`. The checked-in template is intentionally **not valid final metadata**: it has `release_approved=false` and `TBD` author-owned fields.

Before final bundling, the authors must explicitly settle:

- final ordered creators and any chosen ORCID/affiliation values;
- release title and description;
- whether the archive is a new version under the existing Zenodo concept or linked code/data records;
- the data/third-party licensing strategy while keeping original software under MIT;
- the exact final Git commit.

The validator also requires the published v2 record (`10.5281/zenodo.22295791`, concept `10.5281/zenodo.22295790`) to remain unchanged and requires anonymous redownload plus clean replay after publication.

A final metadata contract therefore needs `release_approved=true`, no placeholder values, and `final_code_commit` equal to the checkout used by the bundle builder. This is an approval gate, not a tool for choosing author order or licensing decisions automatically.

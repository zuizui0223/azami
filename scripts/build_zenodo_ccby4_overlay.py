"""Rebuild the 2026-10-30 *unpublished* rights-scoped V4 aggregate candidate.

This is a local archive transformation. It never connects to Zenodo, and does
not establish that material from third-party sources can be redistributed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

ORIGINAL_SHA = "3e2fa32f210977d3aae5a46a9f18e53f88a4ce7ed0af7e4d8057ca0aafb51414"
EXPECTED_SHA = "d4ba5af840f9cd65e4c8abc0d48cdc1becf97800a9baca9573f245e51d914780"


def digest(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def build(source: Path, target: Path, rights_dir: Path) -> dict:
    if source.resolve() == target.resolve() or target.exists():
        raise ValueError("Original must remain unchanged; output must not exist")
    if digest(source.read_bytes()) != ORIGINAL_SHA:
        raise ValueError("Unexpected original source archive SHA-256")
    with zipfile.ZipFile(source) as z:
        if z.testzip() is not None or len(z.namelist()) != 63 or len(z.namelist()) != len(set(z.namelist())):
            raise ValueError("Original ZIP structure invalid")
        files = {name: z.read(name) for name in z.namelist()}
    manifest = json.loads(files["review_manifest.json"])
    prior = {x["path"]: x["sha256"] for x in manifest["files"]}
    if set(prior) != set(files) - {"review_manifest.json"}:
        raise ValueError("Original manifest has missing/unexpected paths")
    if any(digest(files[name]) != value for name, value in prior.items()):
        raise ValueError("Original manifest contains a hash mismatch")
    original_science = {row["path"]: files[row["path"]] for row in manifest["resources"]}

    readme = files["README.md"].decode()
    old = "The licence for original aggregate output data has not been selected for a final\nexternal deposit."
    if readme.count(old) != 1:
        raise ValueError("Original README wording changed; investigate before reuse")
    readme = readme.replace(old, "The authors selected CC BY 4.0 for the original, redistributable derived aggregate\ndata only; the exact terms and exclusions are in `LICENSE_DATA_CC_BY_4_0.md`\nand `DATA_RIGHTS_SCOPE.json`. This does not establish file-by-file third-party\nrights clearance or authorize public release.")
    readme += "\n## Status of this licensed data candidate\n\nThis is an **unpublished, rights-scoped CC BY 4.0 candidate**. Numerical scientific\ncontent is byte-identical to the original 63-member candidate. The existing\nMIT software license is unchanged. The ZIP includes attribution/rights notes\nso a later Zenodo deposit can declare both code and original-data licenses.\nCreator details, anonymous-review arrangements and release approval remain open.\n"
    files["README.md"] = readme.encode()

    notices = files["THIRD_PARTY_NOTICES.md"].decode()
    old_note = "The final licence for author-owned aggregate outputs remains an author decision."
    if notices.count(old_note) != 1:
        raise ValueError("Original license notice changed; investigate before reuse")
    notices = notices.replace(old_note,"The planned licence for author-owned original derived aggregate outputs is CC BY 4.0 (author-approved policy on 2026-10-30); file-level rights clearance and deposit approvals remain outstanding. See `LICENSE_DATA_CC_BY_4_0.md` and `DATA_RIGHTS_SCOPE.json` for limitations.")
    files["THIRD_PARTY_NOTICES.md"] = notices.encode()
    for name in ("LICENSE_DATA_CC_BY_4_0.md", "DATA_RIGHTS_SCOPE.json"):
        files[name] = (rights_dir / name).read_bytes()
    rights = json.loads(files["DATA_RIGHTS_SCOPE.json"])
    if rights.get("author_approved_data_license_policy") != "CC-BY-4.0":
        raise ValueError("Rights scope does not record CC BY 4.0 decision")
    if rights.get("verified_file_by_file_redistribution_rights") is not False:
        raise ValueError("Cannot claim third-party clearance")
    manifest.update({
        "candidate_status":"UNPUBLISHED_CC_BY_4_0_POLICY_APPROVED_RIGHTS_UNCLEARED",
        "data_output_license":"CC-BY-4.0 for original author-controlled derived aggregate data only; third-party rights excluded; file-level review pending",
        "data_license_policy_approved":True,
        "third_party_file_level_rights_cleared":False,
        "public_upload_performed":False,
        "doi":None,
        "author_anonymity_certified":False,
    })
    manifest["files"] = [{"path":n,"size_bytes":len(contents),"sha256":digest(contents)} for n, contents in sorted(files.items()) if n != "review_manifest.json"]
    files["review_manifest.json"] = (json.dumps(manifest, indent=2, ensure_ascii=False)+"\n").encode()
    for name,original in original_science.items():
        if files[name] != original:
            raise ValueError("Scientific numeric bytes modified: "+name)
    for row in manifest["resources"]:
        if digest(files[row["path"]]) != row["sha256"]:
            raise ValueError("Scientific SHA mismatch: "+row["path"])
    with zipfile.ZipFile(target, "w") as z:
        for name,contents in sorted(files.items()):
            info = zipfile.ZipInfo(name, (2026,10,30,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = (0o100644 << 16)
            z.writestr(info, contents, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    with zipfile.ZipFile(target) as z:
        if z.testzip() is not None or len(z.namelist()) != 65:
            raise ValueError("Generated ZIP integrity error")
    if digest(target.read_bytes()) != EXPECTED_SHA:
        raise ValueError("Generated candidate SHA differs from verified reference")
    return {"status":"PASS_UNPUBLISHED_CCBY4_CANDIDATE", "zip_sha256":EXPECTED_SHA,"members":65,"source_sha256":ORIGINAL_SHA,"numeric_resources_preserved":len(original_science),"zenodo_publication_performed":False}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--rights-dir", type=Path, default=Path(__file__).resolve().parents[1]/"manuscripts/v4_s1_s2_frozen_20261025/reproducibility_v4")
    a = p.parse_args(argv)
    print(json.dumps(build(a.source,a.output,a.rights_dir), indent=2))

if __name__ == "__main__":
    main()

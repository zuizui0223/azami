"""Read-only V4 Zenodo release guard. It NEVER uploads/publishes records."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import stat
import sys
import zipfile
from pathlib import Path, PurePosixPath

SHA256 = re.compile(r"^[a-f0-9]{64}$")
SENSITIVE = re.compile(
    r"(^|/)(?:private|secrets?|decrypted|internal_provenance|"
    r"uploader[-_]?ids?|exact[-_]?coordinates)(/|\.|$)|"
    r"\.(?:aesgcm|key|pem)$",
    re.IGNORECASE,
)
REQUIRED_BLOCKERS = (
    "authors_confirmed",
    "third_party_rights_confirmed",
    "peer_review_anonymity_confirmed",
    "separate_zenodo_deposition_authorized",
    "publication_approved",
)


def check_metadata(metadata: dict) -> dict:
    if metadata.get("schema_version") != 1:
        raise ValueError("Unsupported metadata schema")
    if metadata.get("record_type") != "dataset":
        raise ValueError("Wrong resource type: must be an aggregate dataset")
    if metadata.get("status") != "PREPARATION_ONLY_NOT_PUBLISHED":
        raise ValueError("This tool validates unpublished preparation metadata only")
    if metadata.get("published_doi") is not None:
        raise ValueError("V4 publication DOI is not verified by this preparation tool")
    if metadata.get("code_license") != "MIT":
        raise ValueError("Unexpected code-license claim")
    for key in ("canonical_git_commit", "original_S1_S2_frozen_commit"):
        if not re.fullmatch(r"[a-f0-9]{40}", metadata.get(key, "")):
            raise ValueError(f"Invalid Git SHA in {key}")
    packages = metadata.get("candidate_files", [])
    names = [p.get("name") for p in packages]
    if len(packages) != 3 or len(set(names)) != 3:
        raise ValueError("Exactly three distinct documented candidates expected")
    for p in packages:
        if not p["name"].endswith(".zip") or not SHA256.fullmatch(p["sha256"]):
            raise ValueError("Malformed candidate name or SHA-256")
        if not isinstance(p.get("bytes"), int) or p["bytes"] < 1:
            raise ValueError("Invalid candidate size")
        if not isinstance(p.get("members"), int) or p["members"] < 1:
            raise ValueError("Invalid member count")
        if p.get("recommended_publication") is not False:
            raise ValueError("Unapproved candidate marked for automatic publication")
    if metadata.get("data_license_confirmed") is not True:
        raise ValueError("Author approval of the data license has not been recorded")
    if metadata.get("data_license") != "CC-BY-4.0":
        raise ValueError("Unexpected V4 author-controlled aggregate data license")
    if "third_party_rights_confirmed" not in metadata:
        raise ValueError("Explicit third-party rights clearance gate required")
    if not re.fullmatch(r"[0-9a-f]{64}", metadata.get("candidate_files", [{}])[-1].get("sha256", "")):
        raise ValueError("Missing licensed candidate checksum")
    for key in REQUIRED_BLOCKERS:
        if not isinstance(metadata.get(key), bool):
            raise ValueError(f"Gate must be explicit boolean: {key}")
    if not isinstance(metadata.get("creators"), list):
        raise ValueError("Creators must be an explicit list")
    if metadata["authors_confirmed"] and not metadata["creators"]:
        raise ValueError("Confirmed authors require actual creators")
    if metadata["data_license_confirmed"] and not metadata.get("data_license"):
        raise ValueError("Confirmed data license requires its identifier")
    return {
        "status": "PREPARED_NOT_PUBLISHED",
        "blockers": [key for key in REQUIRED_BLOCKERS if not metadata[key]],
        "version": metadata.get("version_label"),
        "candidate_names": names,
        "doi_assigned": False,
    }


def inspect_bundle(path: Path, metadata: dict) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    matched = [p for p in metadata["candidate_files"] if p["name"] == path.name]
    if len(matched) != 1:
        raise ValueError("Candidate ZIP filename not listed in metadata")
    expected = matched[0]
    if path.stat().st_size != expected["bytes"]:
        raise ValueError("Candidate ZIP size mismatch")
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for part in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(part)
    if digest.hexdigest() != expected["sha256"]:
        raise ValueError("Candidate ZIP SHA-256 mismatch")
    with zipfile.ZipFile(path) as archive:
        members = archive.infolist()
        if len(members) != expected["members"]:
            raise ValueError("Candidate ZIP member count mismatch")
        filenames = [x.filename for x in members]
        if len(filenames) != len(set(filenames)):
            raise ValueError("Duplicate ZIP member names")
        for member in members:
            normalized = member.filename.replace("\\", "/")
            p = PurePosixPath(normalized)
            if (
                not normalized or p.is_absolute() or
                ".." in p.parts or any(x in (".", "") for x in p.parts) or
                normalized.startswith("/") or
                re.match(r"^[A-Za-z]:", normalized) or
                stat.S_ISLNK(member.external_attr >> 16) or
                SENSITIVE.search(normalized)
            ):
                raise ValueError(f"Unsafe or identifying ZIP member: {member.filename}")
    if expected.get("licensed_overlay"):
        with zipfile.ZipFile(path) as archive:
            names = set(archive.namelist())
            necessary = {"LICENSE_DATA_CC_BY_4_0.md", "DATA_RIGHTS_SCOPE.json", "review_manifest.json"}
            if not necessary <= names:
                raise ValueError("CC BY 4.0 licensed candidate is missing rights/manifest documentation")
            rights = json.loads(archive.read("DATA_RIGHTS_SCOPE.json"))
            manifest = json.loads(archive.read("review_manifest.json"))
            if rights.get("author_approved_data_license_policy") != "CC-BY-4.0":
                raise ValueError("Rights scope contradicts chosen data license")
            if rights.get("verified_file_by_file_redistribution_rights") is not False:
                raise ValueError("Rights clearance has not been independently verified")
            if manifest.get("data_license_policy_approved") is not True:
                raise ValueError("ZIP manifest is missing the approved policy record")
            declared = {row["path"]: row["sha256"] for row in manifest["files"]}
            if set(declared) != names - {"review_manifest.json"}:
                raise ValueError("ZIP manifest does not cover all files")
            for name, expected_hash in declared.items():
                if hashlib.sha256(archive.read(name)).hexdigest() != expected_hash:
                    raise ValueError("ZIP manifest hash mismatch: " + name)
    return {
        "checked_file": path.name,
        "size": path.stat().st_size,
        "sha256": digest.hexdigest(),
        "member_count": len(filenames),
        "note": "ZIP identity and filenames checked; privacy/rights of file contents require separate builder and legal audits",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--bundle", type=Path, help="Optional local exact-candidate ZIP")
    parser.add_argument("--require-ready", action="store_true",
                        help="Fail unless explicit legal/authorship/publication gates and ZIP bytes are verified")
    args = parser.parse_args(argv)
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    report = check_metadata(metadata)
    if args.bundle:
        report["archive"] = inspect_bundle(args.bundle, metadata)
    else:
        report["archive"] = None
    if args.require_ready and (report["blockers"] or report["archive"] is None):
        report["status"] = "BLOCKED"
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 2
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

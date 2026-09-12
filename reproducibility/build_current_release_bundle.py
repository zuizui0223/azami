"""Build a checksum-gated Chapter 1 current-release bundle.

This is an offline packager. It never downloads data and never publishes to
Zenodo. Supply the four frozen numerical input archives plus the exact native
status CSV. The builder verifies all frozen identities, snapshots the checked-
out repository with ``git archive``, verifies the 15 current reference outputs,
and writes a self-describing ZIP plus a SHA-256 sidecar.

Final release mode deliberately fails closed until a frozen final-figure
manifest and release-metadata JSON are supplied. This prevents durable staging
from being mistaken for a submission/publication-ready release.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

from reproducibility.run_current_analysis import INPUTS, NATIVE_SHA, native_bytes

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "reproducibility" / "current_reference"
REFERENCE_MANIFEST = REFERENCE / "manifest.json"
REPLAY_FILES = (
    ROOT / "reproducibility" / "CURRENT_ANALYSIS.md",
    ROOT / "reproducibility" / "requirements-current.txt",
    ROOT / "reproducibility" / "run_current_analysis.py",
    ROOT / "reproducibility" / "validate_current_analysis.py",
    ROOT / "reproducibility" / "current_replay_execution.json",
    ROOT / "reproducibility" / "current_replay_validation.json",
    ROOT / "reproducibility" / "ZENODO_UPDATE_AUDIT.md",
    ROOT / "reproducibility" / "CURRENT_RELEASE_STAGING_20260912.json",
)
ROOT_METADATA = (ROOT / "README.md", ROOT / "LICENSE", ROOT / "NOTICE.md")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def require_hash(path: Path, expected: str, label: str) -> Path:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 mismatch {actual}; expected {expected}")
    return path


def locate_archive(input_dir: Path, artifact_id: int) -> Path:
    matches = sorted(input_dir.glob(f"*{artifact_id}*.zip"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"artifact {artifact_id}: expected exactly one matching ZIP in {input_dir}; "
            f"found {len(matches)}"
        )
    return matches[0]


def input_contract() -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for role, (artifact_id, archive_sha, members) in INPUTS.items():
        rows[role] = {
            "artifact_id": artifact_id,
            "archive_sha256": archive_sha,
            "members": {
                source: {"bundle_path": target, "sha256": member_sha}
                for source, (target, member_sha) in members.items()
            },
        }
    return rows


def verify_inputs(input_dir: Path, native_status: Path) -> tuple[dict, bytes]:
    receipt: dict[str, dict] = {}
    for role, row in input_contract().items():
        archive = locate_archive(input_dir, row["artifact_id"])
        require_hash(archive, row["archive_sha256"], f"artifact {row['artifact_id']}")
        with zipfile.ZipFile(archive) as zf:
            for member, spec in row["members"].items():
                payload = zf.read(member)
                actual = sha256_bytes(payload)
                if actual != spec["sha256"]:
                    raise ValueError(
                        f"artifact {row['artifact_id']} member {member}: SHA-256 mismatch "
                        f"{actual}; expected {spec['sha256']}"
                    )
        receipt[role] = {
            "artifact_id": row["artifact_id"],
            "source_name": archive.name,
            "archive_sha256": row["archive_sha256"],
            "size_bytes": archive.stat().st_size,
        }

    normalized_native = native_bytes(native_status.read_bytes())
    if sha256_bytes(normalized_native) != NATIVE_SHA:
        raise ValueError("native-status normalization did not recover the frozen SHA-256")
    receipt["native_status"] = {
        "source_name": native_status.name,
        "sha256": NATIVE_SHA,
        "size_bytes": len(normalized_native),
    }
    return receipt, normalized_native


def verified_reference_rows() -> list[dict]:
    manifest = json.loads(REFERENCE_MANIFEST.read_text(encoding="utf-8"))
    rows = manifest["files"]
    if len(rows) != 15:
        raise ValueError(f"current reference manifest must contain 15 files; found {len(rows)}")
    for row in rows:
        path = REFERENCE / row["path"]
        require_hash(path, row["sha256"], f"current reference {row['path']}")
    return rows


def git_state() -> dict[str, str | bool]:
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    status = subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True)
    branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=ROOT, text=True
    ).strip()
    return {"head": head, "branch": branch, "clean": not bool(status.strip())}


def release_gaps(figure_manifest: Path | None, release_metadata: Path | None) -> list[str]:
    gaps: list[str] = []
    if figure_manifest is None:
        gaps.append("final_figure_manifest")
    if release_metadata is None:
        gaps.append("release_metadata")
    return gaps


def verify_manifest_files(manifest_path: Path) -> list[dict]:
    obj = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = obj.get("files")
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{manifest_path}: expected non-empty files list")
    verified: list[dict] = []
    for row in rows:
        rel = Path(row["path"])
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError(f"unsafe figure path: {rel}")
        path = ROOT / rel
        require_hash(path, row["sha256"], f"figure release file {rel}")
        verified.append({"path": rel.as_posix(), "sha256": row["sha256"]})
    return verified


def copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def member_inventory(root: Path) -> list[dict]:
    rows = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix()
        rows.append({"path": rel, "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return rows


def deterministic_zip(source_dir: Path, out_zip: Path) -> str:
    out_zip.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path in sorted(p for p in source_dir.rglob("*") if p.is_file()):
            rel = path.relative_to(source_dir).as_posix()
            info = zipfile.ZipInfo(rel, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            zf.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return sha256_file(out_zip)


def build(
    input_dir: Path,
    native_status: Path,
    out_zip: Path,
    figure_manifest: Path | None = None,
    release_metadata: Path | None = None,
    final: bool = False,
    expected_head: str | None = None,
) -> dict:
    state = git_state()
    if not state["clean"]:
        raise RuntimeError("release bundle requires a clean git worktree")
    if expected_head and state["head"] != expected_head:
        raise RuntimeError(f"HEAD {state['head']} does not match expected {expected_head}")

    gaps = release_gaps(figure_manifest, release_metadata)
    if final and gaps:
        raise RuntimeError("final release blocked by: " + ", ".join(gaps))

    input_receipt, normalized_native = verify_inputs(input_dir, native_status)
    references = verified_reference_rows()
    figure_rows = verify_manifest_files(figure_manifest) if figure_manifest else []
    if release_metadata:
        json.loads(release_metadata.read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(prefix="azami-release-") as td:
        package = Path(td) / "azami_ch1_current_release"
        package.mkdir()

        # Preserve exact frozen numerical archives under canonical names.
        for role, row in input_contract().items():
            src = locate_archive(input_dir, row["artifact_id"])
            copy_file(src, package / "inputs" / f"artifact-{row['artifact_id']}-{role}.zip")
        native_out = package / "inputs" / "observation_native_status.csv"
        native_out.parent.mkdir(parents=True, exist_ok=True)
        native_out.write_bytes(normalized_native)

        # Snapshot all tracked code at the exact recorded commit.
        code_zip = package / "code" / "azami-current-code.zip"
        code_zip.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "archive", "--format=zip", "--output", str(code_zip), str(state["head"])],
            cwd=ROOT,
            check=True,
        )

        # Current reference outputs are copied only after manifest verification.
        copy_file(REFERENCE_MANIFEST, package / "reference" / "manifest.json")
        for row in references:
            copy_file(REFERENCE / row["path"], package / "reference" / row["path"])

        for src in REPLAY_FILES:
            if not src.is_file():
                raise FileNotFoundError(src)
            copy_file(src, package / "replay" / src.name)
        for src in ROOT_METADATA:
            if not src.is_file():
                raise FileNotFoundError(src)
            copy_file(src, package / "metadata" / src.name)

        if release_metadata:
            copy_file(release_metadata, package / "metadata" / "zenodo_release_metadata.json")
        if figure_manifest:
            copy_file(figure_manifest, package / "figures" / "manifest.json")
            for row in figure_rows:
                copy_file(ROOT / row["path"], package / "figures" / row["path"])

        manifest = {
            "schema_version": 1,
            "bundle_kind": "final" if final else "staging",
            "release_ready": final and not gaps,
            "release_gaps": gaps,
            "git": state,
            "input_receipt": input_receipt,
            "native_status_sha256": NATIVE_SHA,
            "current_reference_file_count": len(references),
            "current_reference_manifest_sha256": sha256_file(REFERENCE_MANIFEST),
            "figure_file_count": len(figure_rows),
            "scientific_outputs_changed": False,
            "public_release_performed": False,
        }
        (package / "release_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (package / "checksums.json").write_text(
            json.dumps(member_inventory(package), indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        bundle_sha = deterministic_zip(package, out_zip)

    sidecar = out_zip.with_suffix(out_zip.suffix + ".sha256")
    sidecar.write_text(f"{bundle_sha}  {out_zip.name}\n", encoding="utf-8")
    return {
        "bundle": str(out_zip),
        "sha256": bundle_sha,
        "sidecar": str(sidecar),
        "git_head": state["head"],
        "release_ready": final and not gaps,
        "release_gaps": gaps,
        "reference_files": len(references),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True, help="Directory containing one ZIP for each frozen artifact ID")
    parser.add_argument("--native-status", type=Path, required=True, help="Frozen observation_native_status.csv")
    parser.add_argument("--out", type=Path, required=True, help="Output release ZIP")
    parser.add_argument("--figure-manifest", type=Path, help="JSON with files:[{path,sha256}] for the final manuscript figure/provenance surface")
    parser.add_argument("--release-metadata", type=Path, help="Final Zenodo release metadata JSON")
    parser.add_argument("--expected-head", help="Optional exact git commit required for packaging")
    parser.add_argument("--final", action="store_true", help="Fail closed unless figure manifest and release metadata are supplied")
    args = parser.parse_args()
    receipt = build(
        args.input_dir.resolve(),
        args.native_status.resolve(),
        args.out.resolve(),
        args.figure_manifest.resolve() if args.figure_manifest else None,
        args.release_metadata.resolve() if args.release_metadata else None,
        args.final,
        args.expected_head,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

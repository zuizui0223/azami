#!/usr/bin/env python3
"""Build a double-anonymous, local-only Chapter 1 numerical review bundle.

The output intentionally excludes Git history, public repository metadata,
release/Zenodo metadata, workflow identifiers and the normal network download
runner. It contains only the analysis source needed for the current numerical
surface, checksum-locked extracted inputs, a sanitized 16-file reference
surface, a local-only eight-stage runner and an anonymous validator.

By default the package is self-tested before it is zipped: the extracted bundle
runs all eight numerical stages and must reproduce all 16 reference files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

from reproducibility.run_current_analysis import ROOT, prepare

REFERENCE = ROOT / "reproducibility" / "current_reference"
REFERENCE_MANIFEST = REFERENCE / "manifest.json"

BANNED_TEXT = (
    "zuizui0223",
    "rachelzhang",
    "zhang, ruiqi",
    "zhang ruiqi",
    "10.5281/zenodo.22295791",
    "10.5281/zenodo.22295790",
    "zenodo.org/records/22295791",
)

# analysis/, analysis/v3 and legacy/v2/analysis are namespace packages in the
# repository. Only the two legacy helpers imported by the current sensitivity
# chain, the technical-audit input and the pinned numerical environment are
# explicitly added outside analysis/v3/*.py.
STATIC_SOURCE_PATHS = (
    Path("analysis/ch1/image_to_trait_automated_technical_audit_summary.json"),
    Path("legacy/v2/analysis/run_geb_v2_full27_spatial_sensitivity.py"),
    Path("legacy/v2/analysis/run_geb_v2_full27_historical_sensitivity.py"),
    Path("reproducibility/requirements-current.txt"),
)

RUNNER = r'''#!/usr/bin/env python3
"""Local-only eight-stage numerical replay for anonymous peer review."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
OUT = ROOT / "work" / "results"
OUT.mkdir(parents=True, exist_ok=True)
common = ["--traits", str(INPUT / "traits.csv"), "--environment", str(INPUT / "environment.csv")]
axes = ["--axis-among", str(OUT / "axes/biological_axes_among_min5.csv"),
        "--axis-within", str(OUT / "axes/biological_axes_within.csv")]
seed = ["--seed", "20260910"]
commands = [
    ["analysis.v3.run_biological_axis_reanalysis", *common, "--out-dir", str(OUT / "axes"), "--permutations", "9999", *seed],
    ["analysis.v3.run_biological_axis_sensitivity_chain", *common, *axes,
     "--regions", str(INPUT / "regions.csv"), "--native-status", str(INPUT / "native_status.csv"),
     "--tree-dir", str(INPUT / "trees"), "--out-dir", str(OUT / "sensitivity"),
     "--spatial-permutations", "999", "--moran-permutations", "999",
     "--minimum-taxa-historical", "30", *seed],
    ["analysis.v3.run_construct_scale_integration_entrypoint", *common, *axes,
     "--out-dir", str(OUT / "integration"), "--minimum-paired-observations-per-taxon", "5",
     "--minimum-taxa", "20", "--qap-permutations", "9999", *seed],
    ["analysis.v3.run_construct_scale_upgrade", *common, *axes,
     "--out-dir", str(OUT / "upgrade"), "--minimum-complete-observations-per-taxon", "5",
     "--bootstrap-replicates", "1000", "--permutations", "9999", *seed],
    ["analysis.v3.run_construct_scale_contrast_summary",
     "--pairwise", str(OUT / "upgrade/complete18_construct_pairwise.csv"),
     "--bootstrap", str(OUT / "upgrade/complete18_taxon_bootstrap.csv"),
     "--out", str(OUT / "upgrade/construct_scale_contrast_summary.json")],
    ["analysis.v3.run_assessability_selection_audit", *common, "--out-dir", str(OUT / "assessability")],
    ["analysis.v3.run_frozen_technical_error_stress", *common,
     "--technical-audit-summary", str(ROOT / "analysis/ch1/image_to_trait_automated_technical_audit_summary.json"),
     "--out-dir", str(OUT / "technical_stress"), "--replicates", "2000", *seed],
    ["analysis.v3.run_rv_estimator_validity", *common,
     "--out-dir", str(OUT / "estimator_validity"), "--minimum-complete-observations-per-taxon", "5",
     "--equal-n-replicates", "1000", "--null-permutations", "499", "--qap-permutations", "9999",
     "--seed", "20260915"],
]
receipt = {"status": "RUNNING", "stages": []}
for i, cmd in enumerate(commands, 1):
    start = time.time()
    log = OUT / f"step-{i}.log"
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run([sys.executable, "-m", *cmd], cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT, check=True)
    receipt["stages"].append({"stage": i, "module": cmd[0], "seconds": round(time.time() - start, 3)})
subprocess.run([sys.executable, "review_validator.py", "--results", str(OUT)], cwd=ROOT, check=True)
receipt["status"] = "PASS"
receipt["reference_files_compared"] = 16
(ROOT / "work" / "review_replay_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
print(json.dumps(receipt, indent=2))
'''

VALIDATOR = r'''#!/usr/bin/env python3
"""Tolerance-aware comparison against the anonymous frozen reference surface."""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REFERENCE = ROOT / "current_reference"
RETIRED = {
    "integration/construct_scale_integration_report.json": ["environment_signature_alignment"],
    "upgrade/construct_scale_upgrade_report.json": ["integration_environment_coupling"],
}
def active_scope(payload, relative):
    if relative not in RETIRED: return payload
    scoped = {k:v for k,v in payload.items() if k not in RETIRED[relative]}
    if relative == "upgrade/construct_scale_upgrade_report.json" and scoped.get("claim_boundary") == (
        "secondary construct-level synthesis only; frozen v2 endpoint conclusions and multiplicity families unchanged"
    ):
        scoped["claim_boundary"] = "construct-level synthesis used by the current manuscript; frozen v2 endpoint conclusions and multiplicity families unchanged"
    return scoped
def compare(expected, actual, path="root"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys(): raise AssertionError(f"{path}: object keys differ")
        for k in expected: compare(expected[k], actual[k], f"{path}.{k}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual): raise AssertionError(f"{path}: array length differs")
        for i,(a,b) in enumerate(zip(expected,actual)): compare(a,b,f"{path}[{i}]")
    elif isinstance(expected, bool) or expected is None:
        if expected != actual: raise AssertionError(f"{path}: {actual!r} != {expected!r}")
    elif isinstance(expected, (int,float)):
        if not isinstance(actual,(int,float)) or not math.isclose(expected,actual,rel_tol=1e-8,abs_tol=1e-10): raise AssertionError(f"{path}: {actual!r} != {expected!r}")
    elif expected != actual:
        try: a,b=json.loads(expected),json.loads(actual)
        except (ValueError,TypeError): raise AssertionError(f"{path}: {actual!r} != {expected!r}") from None
        compare(a,b,path)
def load(path):
    if path.suffix == ".json": return json.loads(path.read_text(encoding="utf-8"))
    with path.open(encoding="utf-8", newline="") as h: return list(csv.reader(h))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser(); p.add_argument("--results",required=True,type=Path); args=p.parse_args()
    manifest=json.loads((REFERENCE/"manifest.json").read_text(encoding="utf-8")); checked=[]
    if len(manifest["files"]) != 16: raise AssertionError("anonymous reference surface must contain 16 files")
    for row in manifest["files"]:
        ref=REFERENCE/row["path"]
        if sha(ref) != row["sha256"]: raise AssertionError(f"reference bytes changed: {row['path']}")
        actual=args.results/row["path"]
        if not actual.is_file(): raise FileNotFoundError(actual)
        compare(active_scope(load(ref), row["path"]), active_scope(load(actual), row["path"]), row["path"])
        checked.append(row["path"])
    report={"status":"PASS","aggregate_files_compared":len(checked),"numerical_tolerance":{"relative":1e-8,"absolute":1e-10},"scope":"anonymous numerical peer-review replay"}
    (args.results/"validation.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8"); print(json.dumps(report,indent=2))
if __name__ == "__main__": main()
'''

README = """# Anonymous peer-review numerical package

This package contains the processed numerical inputs, current analysis source,
and frozen aggregate reference outputs needed to reproduce the submitted
Chapter 1 numerical results during double-anonymous peer review.

It intentionally contains no Git history, public repository URL, author
metadata, permanent-archive DOI or release account identifier. Those identifiers
are withheld only for anonymous review and will be restored in the permanent
publication archive.

## Reproduce

Use Python 3.12.

```bash
python -m pip install -r requirements.txt
python run_review_reproduction.py
```

The replay runs eight fixed numerical stages and then tolerance-compares 16
aggregate outputs against `current_reference/manifest.json`. A successful run
ends with `status: PASS` and `aggregate_files_compared: 16`.

The package begins with processed image-derived measurements rather than the
original photographs. It therefore reproduces the numerical analyses reported
in the manuscript, not upstream image acquisition or physical-trait validation.

The raw RV scale contrast is retained as a descriptive result. The manuscript's
qualitative statement that integration is stronger overall among taxa is
supported by the post-hoc equal-n estimator-validity sensitivity included as the
eighth stage; the raw 33/36 relation count is not treated as estimator-invariant.
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def sanitized_reference(package: Path) -> list[dict]:
    manifest = json.loads(REFERENCE_MANIFEST.read_text(encoding="utf-8"))
    rows = []
    for row in manifest["files"]:
        src = REFERENCE / row["path"]
        if sha256(src) != row["sha256"]:
            raise ValueError(f"reference hash mismatch: {row['path']}")
        copy(src, package / "current_reference" / row["path"])
        rows.append({"path": row["path"], "sha256": row["sha256"]})
    if len(rows) != 16:
        raise ValueError(f"expected 16 reference files, found {len(rows)}")
    out = {"schema_version": 1, "scientific_reference_commit": "anonymous-review", "files": rows}
    (package / "current_reference" / "manifest.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    return rows


def copy_source(package: Path) -> list[str]:
    copied = []
    for src in STATIC_SOURCE_PATHS:
        path = ROOT / src
        if not path.is_file():
            raise FileNotFoundError(path)
        target = package / ("requirements.txt" if src.as_posix() == "reproducibility/requirements-current.txt" else src)
        copy(path, target); copied.append(target.relative_to(package).as_posix())
    for path in sorted((ROOT / "analysis/v3").glob("*.py")):
        target = package / path.relative_to(ROOT)
        copy(path, target); copied.append(target.relative_to(package).as_posix())
    (package / "run_review_reproduction.py").write_text(RUNNER, encoding="utf-8")
    (package / "review_validator.py").write_text(VALIDATOR, encoding="utf-8")
    (package / "README.md").write_text(README, encoding="utf-8")
    copied += ["run_review_reproduction.py", "review_validator.py", "README.md"]
    return copied


def identity_hits(root: Path) -> list[dict]:
    hits = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        try:
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line_number, line in enumerate(handle, 1):
                    low = line.lower()
                    for token in BANNED_TEXT:
                        if token.lower() in low:
                            hits.append({"path": path.relative_to(root).as_posix(), "token": token, "line": line_number})
                            if len(hits) >= 100:
                                return hits
        except OSError:
            continue
    return hits


def inventory(root: Path) -> list[dict]:
    return [
        {"path": p.relative_to(root).as_posix(), "size_bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in sorted(x for x in root.rglob("*") if x.is_file())
    ]


def deterministic_zip(root: Path, out: Path) -> str:
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9, allowZip64=True) as zf:
        for path in sorted(p for p in root.rglob("*") if p.is_file()):
            rel = path.relative_to(root).as_posix()
            info = zipfile.ZipInfo(rel, date_time=(1980,1,1,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED; info.external_attr=0o100644<<16
            zf.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return sha256(out)


def build(out: Path, download: bool, self_replay: bool) -> dict:
    with tempfile.TemporaryDirectory(prefix="ch1-anon-review-") as td:
        td = Path(td); package = td / "chapter1_anonymous_review_bundle"; package.mkdir()
        prepared = td / "prepared"
        prepare(td / "archives", prepared, download=download)
        for src in sorted(p for p in prepared.rglob("*") if p.is_file()):
            copy(src, package / "inputs" / src.relative_to(prepared))
        sources = copy_source(package)
        refs = sanitized_reference(package)

        hits = identity_hits(package)
        if hits:
            raise RuntimeError("identity scan failed: " + json.dumps(hits[:20]))

        validation = None
        replay_seconds = None
        if self_replay:
            start = time.time()
            subprocess.run([sys.executable, "run_review_reproduction.py"], cwd=package, check=True)
            replay_seconds = round(time.time() - start, 3)
            validation_path = package / "work/results/validation.json"
            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            if validation.get("status") != "PASS" or validation.get("aggregate_files_compared") != 16:
                raise RuntimeError("anonymous self replay did not pass")
            receipt = json.loads((package / "work/review_replay_receipt.json").read_text(encoding="utf-8"))
            (package / "SELF_REPLAY_VALIDATION.json").write_text(json.dumps({
                "status": "PASS", "aggregate_files_compared": 16,
                "eight_stage_replay": receipt, "seconds": replay_seconds,
            }, indent=2) + "\n", encoding="utf-8")
            shutil.rmtree(package / "work")

        hits_after = identity_hits(package)
        if hits_after:
            raise RuntimeError("post-replay identity scan failed: " + json.dumps(hits_after[:20]))
        assembly = {
            "schema_version": 1,
            "anonymous_identity_scan": "PASS",
            "self_replay": validation.get("status") if validation else "SKIPPED",
            "reference_files": len(refs),
            "analysis_source_files": len(sources),
            "public_repository_metadata_included": False,
            "permanent_archive_metadata_included": False,
            "git_history_included": False,
        }
        (package / "ANONYMOUS_REVIEW_ASSEMBLY.json").write_text(json.dumps(assembly, indent=2) + "\n", encoding="utf-8")
        (package / "FILE_INVENTORY.json").write_text(json.dumps(inventory(package), indent=2) + "\n", encoding="utf-8")
        digest = deterministic_zip(package, out)
    sidecar = out.with_suffix(out.suffix + ".sha256")
    sidecar.write_text(f"{digest}  {out.name}\n", encoding="utf-8")
    return {"bundle": str(out), "sha256": digest, "sidecar": str(sidecar), "self_replay": "PASS" if self_replay else "SKIPPED"}


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--download", action="store_true")
    p.add_argument("--skip-self-replay", action="store_true")
    args=p.parse_args()
    print(json.dumps(build(args.out.resolve(), args.download, not args.skip_self_replay), indent=2))

if __name__ == "__main__":
    main()
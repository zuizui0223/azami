#!/usr/bin/env python3
"""Build an offline V4 aggregate-verification candidate and a review copy.

This does not publish, create a DOI, use credentials, restore confidential
observations, or refit observation-level models. Every input archive is pinned
and every copied file is allowlisted. Original MIT notices remain intact.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tempfile
import tokenize
import zipfile


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT = Path("manuscripts/v4_s1_s2_frozen_20261025")
SOURCE_SHA = "6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205"
ROBUSTNESS_SHA = "216fd6b8b722dd615ebbce9f58412ce245477c25ec65057afe9b5962cee33a91"
ASSET_SHA = "0c6879ebddbf8ba56363049ca09fd69c7260150d0b38bcdcf52fe59f285a5c32"
SOURCE_MEMBERS = {
    "S1_S2_extreme_LEGO_results.json", "V4_S1_S2_four_model_comparison.json",
    "S1_extreme_LEGO_aggregate.json", "V4_S1_S2_four_model_comparison_JA.md",
    "V4_S1_S2_four_vs_vif5_vif10_alltraits.json",
    "V4_S1_S2_full_hypothesis_FDR.json", "V4_S1_module_label_null.json",
    "environment_source_receipt.json",
}
ROBUSTNESS_MEMBERS = {
    "H3C_CONDITIONAL_KDE_TRAIT_OVERLAP_AGGREGATE.json",
    "H3C_ROBUSTNESS_SENSITIVITY_AGGREGATE.json",
    "H3C_MONTE_CARLO_SEED_STABILITY_AGGREGATE.json",
    "H3D_CLIMATE_MODERATED_TRAIT_COVARIANCE_AGGREGATE.json",
    "H3_CLIMATE_TRAIT_HYPERVOLUME_AGGREGATE.json",
}
TABLES = (
    ("R001", "Table_S1_full_208_BIO18_tests.csv", 208),
    ("R002", "Table_S2_all_28_RV_pairs.csv", 28),
    ("R003", "Table_S2_S2_original_native_introduced_colour_coefficients.csv", 6),
    ("R004", "SI_S1_all_80_original_raw_BH_Holm.csv", 80),
    ("R005", "SI_S1_BIO18_image_angle_both_scales.csv", 2),
)
JSON_RESOURCES = (
    ("R006", "analysis/v4/results/20261029_H3C_CLIMATE_ANALOGUE_TRAIT_OVERLAP_SOURCE_VERIFIED.json", "Original exploratory H3c distance and density metrics"),
    ("R007", "analysis/v4/results/20261029_H3_CLIMATE_TRAIT_HYPERVOLUME_SOURCE_VERIFIED.json", "H3 climatic and photographic breadth"),
    ("R008", "analysis/v4/results/20261029_H3D_ENVIRONMENT_MODERATED_TRAIT_COVARIANCE_SOURCE_VERIFIED.json", "H3d trait and covariance context"),
    ("R009", "analysis/v4/results/20261029_S1_BIO18_GEOGRAPHY_SOURCE_VERIFIED.json", "S1 geographic-regressor sensitivity"),
    ("R010", "analysis/v4/results/20261029_S1_BIO18_SPATIAL_BLOCK_SOURCE_VERIFIED.json", "S1 geographic block bootstrap"),
    ("R011", "analysis/v4/results/20261029_S2_COLOUR_SPATIAL_BLOCK_SOURCE_VERIFIED.json", "S2 six-row geographic block bootstrap"),
    ("R012", "analysis/v4/results/20261021_unified_abcde/S2_ALL_350_SPECIES_SOURCE_SCREEN_VERIFIED_20261021.json", "Whole-source focal-species support screen"),
)
ROBUST_RESOURCES = (
    ("R013", "H3C_ROBUSTNESS_SENSITIVITY_AGGREGATE.json", "Post-positive H3c threshold, trait, region and bandwidth diagnostics"),
    ("R014", "H3C_MONTE_CARLO_SEED_STABILITY_AGGREGATE.json", "Post-positive H3c six-seed stability"),
)
CODE = ("analysis/v4/final_v4_bio18_holm_freeze.py",) + tuple(
    (MANUSCRIPT / name).as_posix() for name in (
        "build_full_supplement.py", "verify_s2_frozen_colour_context.py",
        "plot_frozen_summary.py", "plot_bio18_angle_part2.py",
        "plot_bio18_geographic_sensitivity.py", "plot_s2_native_introduced_colour.py",
        "plot_cross_functional_pairs.py", "plot_h3_climate_morphospace.py",
    )
)
COMPANION_CSV = ("table1_holm_associations.csv", "module_summary.csv", "s2_module_summary.csv")
FIGURES = tuple(f"{folder}/{name}.{ext}" for folder, names in (
    ("figures", ("SI_H3_context_distance_overlap_breadth", "SI_S2_trait_domain_RV",
                  "SI_four_Holm_supported_coefficients", "SI_module_RV_medians",
                  "figure2_cross_functional_trait_pairs", "figure3_BIO18_among_taxa_image_angle",
                  "figure4_S2_native_introduced_colour_coefficients",
                  "figure5_H3_climatic_breadth_floral_morphospace")),
    ("si", ("Figure_S1_S1_all_28_RV", "SI_Figure_BIO18_image_angle_both_scales",
             "figure3_BIO18_geographic_sensitivity")),
) for name in names for ext in ("png", "svg"))
ASSET_MEMBERS = set(FIGURES) | {"freeze_check.json", "si/SOURCE_SHA256.txt"} | {
    f"si/{name}" for _, name, _ in TABLES
}
PROTECTED_KEYS = {
    "obs_id", "observation_id", "user_id", "accepted_key", "photo_id",
    "latitude", "longitude", "photographer_id", "observer_id", "uploader_id",
    "decimal_latitude", "decimal_longitude", "decimallatitude", "decimallongitude",
    "photo_url", "image_url", "original_photo_url", "observation_url",
    "observer_name", "photographer_name", "uploader_name", "user_login",
    "username", "user_name", "photo_filename", "image_filename",
    "original_image_url", "observation_uuid", "lat", "lon", "lng",
}
AUDIT_ID_KEYS = {
    "run_id", "job_id", "artifact_id", "github_actions_run", "original_private_job_id",
    "aggregate_artifact_id", "original_github_artifact_id", "public_artifact_id",
    "original_secure_workflow_run", "authenticated_original_run", "cache_artifact_id",
}
IDENTITY = re.compile(r"github(?:usercontent)?\.com|api\.github\.com|zuizui0223|rachelzhang|zhang\s*,?\s*ruiqi|zenodo\.(?:org|2229579[01])", re.I)
OBSERVATION_URL = re.compile(r"(?:inaturalist\.[^/\s\"']+)/(?:observations|photos)/\d+", re.I)
REPO_URL = re.compile(r"https?://(?:www\.)?(?:api\.)?(?:raw\.)?github(?:usercontent)?\.com/[^\s\"'<>)]*", re.I)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def json_bytes(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode()


def check_data(name: str, raw: bytes) -> None:
    """Reject protected record fields and individual links in data, recursively."""
    if name.endswith(".png"):
        if IDENTITY.search(raw.decode("latin1")):
            raise ValueError(f"identifying image metadata: {name}")
        return
    text = raw.decode("utf-8")
    if OBSERVATION_URL.search(text):
        raise ValueError(f"individual observation/photo link: {name}")
    def visit(value, path):
        if isinstance(value, dict):
            for key, child in value.items():
                if key.lower() in PROTECTED_KEYS:
                    raise ValueError(f"protected data field: {name}:{path}.{key}")
                visit(child, path + "." + key)
        elif isinstance(value, list):
            for child in value:
                visit(child, path + "[]")
    if name.endswith(".json"):
        visit(json.loads(text), "root")
    elif name.endswith(".csv"):
        reader = csv.DictReader(io.StringIO(text))
        for key in reader.fieldnames or []:
            if key.lower() in PROTECTED_KEYS:
                raise ValueError(f"protected CSV column: {name}:{key}")


def read_archive(path: Path, expected_sha: str, members: set[str]) -> dict[str, bytes]:
    raw = path.read_bytes()
    if sha(raw) != expected_sha:
        raise ValueError(f"archive SHA-256 mismatch: {path.name}")
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)) or set(names) != members:
            raise ValueError(f"unexpected/duplicate archive members: {path.name}")
        result = {}
        for item in archive.infolist():
            p = PurePosixPath(item.filename)
            if p.is_absolute() or ".." in p.parts or item.is_dir() or item.file_size > 20_000_000:
                raise ValueError(f"unsafe archive member: {item.filename}")
            if (item.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError(f"archive symlink: {item.filename}")
            data = archive.read(item)
            check_data(item.filename, data)
            result[item.filename] = data
        return result


def zip_bytes(files: dict[str, bytes]) -> bytes:
    out = io.BytesIO()
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, raw in sorted(files.items()):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, raw, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return out.getvalue()


def audit_ids(value, found: set[int]) -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            if key.lower() in AUDIT_ID_KEYS and isinstance(child, int) and not isinstance(child, bool):
                found.add(child)
            audit_ids(child, found)
    elif isinstance(value, list):
        for child in value:
            audit_ids(child, found)


def review_text(text: str, identifiers: dict[int, str], hashes: dict[str, str]) -> str:
    text = REPO_URL.sub("source-audit-withheld-for-review", text)
    text = re.sub(r"https?://(?:doi\.org/10\.5281/)?zenodo[^\s\"'<>)]*", "archive-identity-withheld", text, flags=re.I)
    text = re.sub(r"zuizui0223/azami|zuizui0223|rachelzhang|zhang\s*,?\s*ruiqi", "identity-withheld", text, flags=re.I)
    text = re.sub(r"\bazami\b", "study", text, flags=re.I)
    text = text.replace("azami_", "study_").replace("AZAMI_", "STUDY_")
    text = re.sub(r"\b[0-9a-f]{40}\b", "source-revision-withheld", text)
    for original, alias in identifiers.items():
        text = re.sub(rf"(?<!\d){original}(?!\d)", alias, text)
    for original, replacement in hashes.items():
        text = text.replace(original, replacement)
    return text


def assert_scientific_numbers(original, reviewed, key: str = "", path: str = "root") -> None:
    if isinstance(original, dict):
        if not isinstance(reviewed, dict) or original.keys() != reviewed.keys():
            raise ValueError(f"scientific structure changed: {path}")
        for k, value in original.items():
            assert_scientific_numbers(value, reviewed[k], k, path + "." + k)
    elif isinstance(original, list):
        if not isinstance(reviewed, list) or len(original) != len(reviewed):
            raise ValueError(f"scientific list changed: {path}")
        for a, b in zip(original, reviewed):
            assert_scientific_numbers(a, b, key, path + "[]")
    elif isinstance(original, (int, float, bool)) and key.lower() not in AUDIT_ID_KEYS:
        if type(original) is not type(reviewed) or original != reviewed:
            raise ValueError(f"scientific number changed: {path}")


def review_json(value, identifiers: dict[int, str], hashes: dict[str, str], ledger: list,
                path: str = "root", key: str = ""):
    if isinstance(value, dict):
        return {k: review_json(v, identifiers, hashes, ledger, path + "." + k, k)
                for k, v in value.items()}
    if isinstance(value, list):
        return [review_json(v, identifiers, hashes, ledger, f"{path}[{i}]", key)
                for i, v in enumerate(value)]
    result = value
    if key.lower() in AUDIT_ID_KEYS and isinstance(value, int) and value in identifiers:
        result = identifiers[value]
    elif isinstance(value, str):
        result = review_text(value, identifiers, hashes)
    if result != value:
        ledger.append({"json_path": path, "original": value, "review": result,
                       "reason": "identifying provenance or copied checksum only"})
    return result


def review_code(raw: bytes, identifiers: dict[int, str], hashes: dict[str, str]) -> bytes:
    tokens = []
    for token in tokenize.tokenize(io.BytesIO(raw).readline):
        value = token.string
        if token.type == tokenize.NUMBER and value.isdigit() and int(value) in identifiers:
            token = token._replace(type=tokenize.STRING, string=repr(identifiers[int(value)]))
        elif token.type in (tokenize.STRING, tokenize.COMMENT):
            token = token._replace(string=review_text(value, identifiers, hashes))
        tokens.append(token)
    out = tokenize.untokenize(tokens)
    compile(out, "review-copy.py", "exec")
    return out


REPLAY = r'''#!/usr/bin/env python3
"""Verify packaged aggregate bytes and regenerate all five CSVs and figures offline."""
from pathlib import Path
import argparse, csv, hashlib, json, os, subprocess, sys, tempfile

BASE = Path(__file__).resolve().parent
M = "manuscripts/v4_s1_s2_frozen_20261025/"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path)
    a=p.parse_args()
    manifest=json.loads((BASE/"review_manifest.json").read_text())
    for row in manifest["files"]:
        path=BASE/row["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row["sha256"]:
            raise ValueError("PAYLOAD_CHECKSUM_CHANGED: "+row["path"])
    out=(a.output_dir or BASE/"reproduced").resolve()
    out.mkdir(parents=True, exist_ok=True)
    si=out/"si"; figures=out/"figures"
    si.mkdir(exist_ok=True); figures.mkdir(exist_ok=True)
    archive=BASE/"inputs/frozen_primary_aggregates.zip"
    commands=[
      ["analysis/v4/final_v4_bio18_holm_freeze.py","--artifact-zip",archive,"--output",out/"freeze_check.json"],
      [M+"build_full_supplement.py","--artifact",archive,"--out",si],
      [M+"verify_s2_frozen_colour_context.py","--artifact",archive,"--output",si/"Table_S2_S2_original_native_introduced_colour_coefficients.csv"],
      [M+"plot_frozen_summary.py","--output-dir",figures],
      [M+"plot_bio18_angle_part2.py","--artifact",archive,"--output-dir",figures,"--si-dir",si],
      [M+"plot_bio18_geographic_sensitivity.py","--output-dir",si],
      [M+"plot_s2_native_introduced_colour.py","--source",si/"Table_S2_S2_original_native_introduced_colour_coefficients.csv","--output-dir",figures],
      [M+"plot_cross_functional_pairs.py","--pairs-csv",si/"Table_S2_all_28_RV_pairs.csv","--output-dir",figures],
      [M+"plot_h3_climate_morphospace.py","--output-dir",figures],
    ]
    env=os.environ.copy();env["PYTHONDONTWRITEBYTECODE"]="1"
    env["MPLCONFIGDIR"]=str(out/"matplotlib-cache")
    for cmd in commands:
        proc=subprocess.run([sys.executable,*map(str,cmd)],cwd=BASE,env=env,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        if proc.returncode:
            raise RuntimeError(str(cmd[0])+" failed:\n"+proc.stdout+proc.stderr)
    checked=[]
    for resource in manifest["resources"]:
        if resource.get("row_count") is None:
            continue
        expected=BASE/resource["path"];actual=si/expected.name
        if actual.read_bytes()!=expected.read_bytes():
            raise ValueError("REPRODUCED_TABLE_DIFFERS: "+expected.name)
        with actual.open(newline="") as f:count=sum(1 for _ in csv.DictReader(f))
        if count!=resource["row_count"]:
            raise ValueError("TABLE_ROW_COUNT_CHANGED: "+expected.name)
        checked.append({"resource_id":resource["id"],"rows":count,
                        "sha256":hashlib.sha256(actual.read_bytes()).hexdigest()})
    needed=[figures/"figure3_BIO18_among_taxa_image_angle.png",
            figures/"figure4_S2_native_introduced_colour_coefficients.png",
            figures/"figure5_H3_climatic_breadth_floral_morphospace.png",
            si/"Figure_S1_S1_all_28_RV.png"]
    if not all(p.is_file() and p.stat().st_size>1000 for p in needed):
        raise ValueError("EXPECTED_NUMERICAL_FIGURE_MISSING")
    report={"status":"PASS","scope":"aggregate verification and re-rendering only",
            "observation_level_models_refitted":False,"tables_byte_identical":checked,
            "completed_commands":len(commands),"required_numerical_figures_present":len(needed),
            "figure_pixels_byte_compared":False,
            "limitation":"Original confidential observations are absent; no image measurement or source-model refitting is performed."}
    (out/"replay_validation.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report))

if __name__=="__main__":main()
'''

README = """# Cirsium V4 aggregate verification candidate

This local candidate verifies archived aggregate results and redraws numerical
figures. It does not refit models from observation-level inputs or reproduce
image measurement. Protected observations, coordinates, uploader/photo IDs,
photographs, raster bytes, trained weights and repository history are absent.

## Use

Use Python 3.12. Install `requirements.txt` in an isolated environment, then run
`python replay_aggregates.py`. Once dependencies are available, this command
needs no network or credentials. It verifies file hashes, reruns nine aggregate
checks/renderers, and byte-compares five complete CSVs (208, 28, 6, 80 and 2 rows).
The output is written to `reproduced/`; its validation receipt explicitly states
that original-source models were not refitted. Figure rendering is checked, but
pixel identity is not promised across environments.

## Evidence and interpretation

`review_manifest.json` maps Sources R001–R014 to real files and their hashes.
The primary 24 October results use separate 80-test S1 and 128-test S2 Holm
families. The 29 October H3 family and spatial/seed sensitivities are separate,
exploratory analyses. The two post-positive H3c objects preserve every condition
and all six seeds; they are not independent source replications. The original
full aggregate archive includes exploratory context outside the primary model.
Its presence does not expand the manuscript's claims.

The historical `figure2_cross_functional_trait_pairs` file is a supplementary
cross-domain display. The current geographic main Figure 2 is source-verified
separately, using pooled sampling-coverage cells and its own provenance receipts.
That map and the Figure 1 photograph panel are outside this numerical aggregate
package.

## Rights and review status

Original software is MIT licensed; `LICENSE` and `NOTICE.md` are copied exactly.
They retain required legal attribution. Therefore this is a review candidate,
and **full author anonymity is not certified**. Repository and workflow links
are removed from the review copy, but legal attribution is not deleted.

The licence for original aggregate output data has not been selected for a final
external deposit. Third-party data/image licences remain separate; no upstream
material is relicensed by the software MIT licence. `THIRD_PARTY_NOTICES.md`
records this boundary. No DOI, public upload or reviewer-access URL is created
by this builder, and historical deposit approval is not applied to this release.
"""

NOTICES = """# Third-party and aggregate-output rights

The root MIT licence covers original software/software documentation only.
The unchanged NOTICE.md defines its scope. The archive contains pooled author
analysis outputs, numerical figure displays and code, not original photographs,
coordinate records, uploader information, climate rasters or trained models.
The final licence for author-owned aggregate outputs remains an author decision.

Source photographs and observation metadata have separate item-specific terms:
https://help.inaturalist.org/en/support/solutions/articles/151000175695
WorldClim's current terms require permission for data redistribution; published
research maps are treated separately: https://worldclim.org/about.html
The current CHELSA climatology and bioclim dataset pages state CC0; model software
has separate licensing: https://www.chelsa-climate.org/datasets/chelsa_climatologies
and https://www.chelsa-climate.org/datasets/chelsa_bioclim
Kew data licences are dataset-specific; retain the exact authority provenance:
https://www.kew.org/science/collections-and-resources/data-and-digital/terms-of-use

No figure containing an original source photograph or geographic point records
is included. This candidate is not evidence that source-data redistribution
permissions or journal data-policy requirements have been satisfied.
"""


def manifest(files: dict[str, bytes], resources: list[dict], review: bool) -> dict:
    return {
        "schema_version": 1, "candidate_status": "LOCAL_UNPUBLISHED_REVIEW_CANDIDATE",
        "scope": "aggregate verification and numerical re-rendering",
        "observation_level_refitting_supported": False,
        "public_upload_performed": False, "doi": None, "review_access_url": None,
        "author_anonymity_certified": False,
        "identity_handling": "Repository/workflow metadata neutralized; mandatory MIT legal notices unchanged" if review else "Full source provenance kept in companion internal_provenance.json",
        "data_output_license": "PENDING_AUTHOR_DECISION",
        "software_license": "MIT; see unchanged LICENSE and NOTICE.md",
        "excluded": ["original observations", "precise source coordinates", "uploader or photo identifiers",
                     "photographs and crops", "encrypted private inputs", "raster bytes", "model weights", "Git history"],
        "resources": resources,
        "files": [{"path": p, "size_bytes": len(raw), "sha256": sha(raw)}
                  for p, raw in sorted(files.items()) if p != "review_manifest.json"],
        "manifest_self_hash": "Excluded to avoid a circular digest; final ZIP SHA-256 is in its separate sidecar.",
    }


def write_files(folder: Path, files: dict[str, bytes]) -> None:
    folder.mkdir(parents=True, exist_ok=False)
    for rel, raw in files.items():
        path = folder / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)


def validate_candidate(folder: Path) -> dict:
    with tempfile.TemporaryDirectory(prefix="v4-aggregate-replay-") as tmp:
        proc = subprocess.run([sys.executable, str(folder / "replay_aggregates.py"),
                               "--output-dir", tmp], cwd=folder,
                              capture_output=True, text=True, timeout=240,
                              env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        if proc.returncode:
            raise RuntimeError("offline aggregate replay failed:\n" + proc.stdout + proc.stderr)
        result = json.loads((Path(tmp) / "replay_validation.json").read_text())
        if result["status"] != "PASS" or len(result["tables_byte_identical"]) != 5:
            raise ValueError("aggregate self-replay did not pass")
        return result


def build(args) -> dict:
    root = args.repo_root.resolve()
    out = args.output_dir.resolve()
    if out.exists():
        raise ValueError("output directory already exists; choose a fresh candidate directory")
    original = read_archive(args.source_aggregate_zip, SOURCE_SHA, SOURCE_MEMBERS)
    robust = read_archive(args.robustness_zip, ROBUSTNESS_SHA, ROBUSTNESS_MEMBERS)
    assets = read_archive(args.verified_assets_zip, ASSET_SHA, ASSET_MEMBERS)
    resources = []
    common = {}
    source_mapping = {}
    for resource_id, name, count in TABLES:
        path = f"si/{name}"
        raw = assets[path]
        rows = list(csv.DictReader(io.StringIO(raw.decode())))
        if len(rows) != count:
            raise ValueError(f"wrong table row count: {name}")
        common[path] = raw
        resources.append({"id": resource_id, "path": path, "row_count": count, "sha256": sha(raw)})
    for resource_id, path, title in JSON_RESOURCES:
        raw = (root / path).read_bytes()
        check_data(path, raw)
        common[path] = raw
        resources.append({"id": resource_id, "path": path, "description": title, "sha256": sha(raw)})
    for resource_id, name, title in ROBUST_RESOURCES:
        path = "aggregates/" + name
        common[path] = robust[name]
        resources.append({"id": resource_id, "path": path, "description": title, "sha256": sha(robust[name])})
    # Preserve all five full H3 output objects, including their exploratory context.
    for name, raw in robust.items():
        common["aggregates/" + name] = raw
    common["inputs/frozen_primary_aggregates.zip"] = args.source_aggregate_zip.read_bytes()
    for path in CODE:
        common[path] = (root / path).read_bytes()
    for name in COMPANION_CSV:
        path = (MANUSCRIPT / name).as_posix()
        common[path] = (root / path).read_bytes()
        check_data(path, common[path])
    for path in FIGURES:
        common[path] = assets[path]
    common["validation/original_freeze_check.json"] = assets["freeze_check.json"]
    for path in ("LICENSE", "NOTICE.md"):
        common[path] = (root / path).read_bytes()
    common["README.md"] = README.encode()
    common["THIRD_PARTY_NOTICES.md"] = NOTICES.encode()
    common["replay_aggregates.py"] = REPLAY.encode()
    packages = {name: importlib.metadata.version(name) for name in (
        "numpy", "scipy", "pandas", "statsmodels", "patsy", "matplotlib", "Pillow", "packaging",
        "contourpy", "cycler", "fonttools", "kiwisolver", "pyparsing", "python-dateutil", "six",
    )}
    common["requirements.txt"] = ("# Validated with Python 3.12; isolated environment recommended.\n" +
                                    "\n".join(f"{k}=={v}" for k, v in sorted(packages.items())) + "\n").encode()
    common["SOFTWARE_ENVIRONMENT.json"] = json_bytes({"python": sys.version.split()[0], "packages": packages})
    common["DATA_DICTIONARY.md"] = (
        "# Aggregate data dictionary\n\n"
        "R001: all 208 source model tests; family, group, photographic coordinate, environmental predictor, scale, test kind, n, standardized beta, pointwise Wald CI, raw/BH/Holm P, and trait-cohort VIF flag. Hue is one joint test; unavailable scalar beta/CI fields are blank. S1=80 and S2=128 are separate families.\n\n"
        "R002: 28 unordered trait pairs, shared n_photos=156683/n_taxa=174, raw/adjusted within/among RV, signed r and among-minus-within. These are descriptive image metrics.\n\n"
        "R003: six original native/introduced colour coefficients; R004: complete S1 80 tests; R005: the two S1 image-angle/BIO18 scales. R003 and R005 are subsets, not new test families.\n\n"
        "R006–R008: exploratory H3 result objects. R009–R011: spatial sensitivities. R012: aggregate support counts, including species names and taxonomic authority keys, no individual observers or positions. R013–R014: full post-positive conditions and all six random sequences, not independent source replications.\n\n"
        "BIO18 is warmest-quarter mean monthly precipitation; BIO05 is maximum warmest-month temperature; RSDS is downwelling shortwave radiation; elevation is metres. Model beta values are standardized. 'coordinate' names a photographic metric, not a geographic record. Refer to source JSON method/design fields for units and denominators of each H3 statistic.\n"
    ).encode()
    ids = {40062126006, 12547876704, 40859868003, 12891286463,
           40837107465, 12881655588, 40841007115, 12883021095,
           40847696537, 12885860903}
    for path, raw in {**common, **original}.items():
        if path.endswith(".json"):
            audit_ids(json.loads(raw), ids)
    identifier_map = {value: f"audit-{i:03d}" for i, value in enumerate(sorted(ids), 1)}
    reviewed_members = {}
    changes = []
    hash_map = {}
    for name, raw in original.items():
        ledger = []
        if name.endswith(".json"):
            parsed = json.loads(raw)
            reviewed = review_json(parsed, identifier_map, {}, ledger)
            assert_scientific_numbers(parsed, reviewed)
            result = json_bytes(reviewed)
        else:
            result = review_text(raw.decode(), identifier_map, {}).encode()
        reviewed_members[name] = result
        hash_map[sha(raw)] = sha(result)
        changes.append({"path": "inputs/frozen_primary_aggregates.zip!" + name,
                        "original_sha256": sha(raw), "review_sha256": sha(result), "changes": ledger})
    reviewed_zip = zip_bytes(reviewed_members)
    hash_map[SOURCE_SHA] = sha(reviewed_zip)
    review = {}
    for path, raw in common.items():
        ledger = []
        if path == "inputs/frozen_primary_aggregates.zip":
            result = reviewed_zip
        elif path in ("LICENSE", "NOTICE.md") or path.endswith((".csv", ".png")):
            result = raw
        elif path.endswith(".json"):
            parsed = json.loads(raw)
            reviewed = review_json(parsed, identifier_map, hash_map, ledger)
            assert_scientific_numbers(parsed, reviewed)
            result = json_bytes(reviewed)
        elif path.endswith(".py"):
            result = review_code(raw, identifier_map, hash_map)
        else:
            result = review_text(raw.decode(), identifier_map, hash_map).encode()
        review[path] = result
        source_mapping[path] = {"original_sha256": sha(raw), "review_sha256": sha(result)}
        if result != raw:
            changes.append({"path": path, **source_mapping[path], "changes": ledger,
                            "copied_code_checksum_constants_updated": path.endswith(".py")})
    # Author/project attribution is a deliberate legal exception, never removed.
    identity_hits = []
    for path, raw in review.items():
        if path.endswith(".zip"):
            content = "\n".join(v.decode() for v in reviewed_members.values())
        elif path.endswith(".png"):
            continue
        else:
            content = raw.decode()
        if path not in ("LICENSE", "NOTICE.md") and (IDENTITY.search(content) or re.search(r"\bazami\b|azami_", content, re.I)):
            identity_hits.append(path)
    if identity_hits:
        raise ValueError("identifying metadata remains outside legal notices: " + ", ".join(identity_hits))
    review_resources = [{**item, "sha256": sha(review[item["path"]])} for item in resources]
    common["review_manifest.json"] = json_bytes(manifest(common, resources, False))
    review["review_manifest.json"] = json_bytes(manifest(review, review_resources, True))
    out.mkdir(parents=True)
    internal_dir = out / "aggregate_release_candidate"
    review_dir = out / "review_aggregate_candidate"
    write_files(internal_dir, common)
    write_files(review_dir, review)
    internal_validation = validate_candidate(internal_dir)
    review_validation = validate_candidate(review_dir)
    for files, folder, res, is_review, validation in (
        (common, internal_dir, resources, False, internal_validation),
        (review, review_dir, review_resources, True, review_validation),
    ):
        files["BUILD_VALIDATION.json"] = json_bytes(validation)
        files["review_manifest.json"] = json_bytes(manifest(files, res, is_review))
        (folder / "BUILD_VALIDATION.json").write_bytes(files["BUILD_VALIDATION.json"])
        (folder / "review_manifest.json").write_bytes(files["review_manifest.json"])
    provenance = {
        "schema_version": 1, "scope": "INTERNAL_IDENTIFYING_AUDIT_MAP_NOT_IN_REVIEW_ZIP",
        "repository": "https://github.com/zuizui0223/azami",
        "code_base_commit": "90c21b0d52e6c983c6a9761d8569a7edc6f2d1a0",
        "candidate_builder_sha256": sha(Path(__file__).read_bytes()),
        "scientific_source_commit": "d2c45b00018a0c04e7c69a29857207dc1ac06d0c",
        "source_archives": [
            {"run": 40062126006, "artifact": 12547876704, "sha256": SOURCE_SHA,
             "members": {k: sha(v) for k, v in original.items()}},
            {"run": 40859868003, "artifact": 12891286463, "sha256": ROBUSTNESS_SHA,
             "members": {k: sha(v) for k, v in robust.items()}},
            {"artifact": 12956327726, "sha256": ASSET_SHA,
             "members": {k: sha(v) for k, v in assets.items()}},
        ],
        "review_identifier_crosswalk": {alias: value for value, alias in identifier_map.items()},
        "copied_file_crosswalk": source_mapping, "review_changes": changes,
        "scientific_numeric_preservation": "PASS for every JSON numeric leaf except explicitly identified workflow/run/artifact metadata",
        "privacy_scope": "13 source archive members and explicit tracked/generated allowlist checked; no observation-level data included",
        "review_repository_identity_scan": "PASS outside unchanged legal notices",
        "legal_attribution_conflict": "Unchanged MIT and NOTICE identify Azami contributors; full anonymity not certified",
        "historical_zenodo_approval_applied": False, "published": False, "doi": None,
    }
    (out / "internal_provenance.json").write_bytes(json_bytes(provenance))
    archives = []
    for name, files in (("cirsium_v4_aggregate_release_candidate.zip", common),
                        ("cirsium_v4_review_aggregate_candidate.zip", review)):
        raw = zip_bytes(files)
        (out / name).write_bytes(raw)
        (out / (name + ".sha256")).write_text(sha(raw) + "  " + name + "\n")
        archives.append({"file": name, "sha256": sha(raw), "size_bytes": len(raw), "files": len(files)})
    report = {
        "status": "LOCAL_AGGREGATE_CANDIDATE_VERIFIED", "archives": archives,
        "original_and_review_replay": "PASS: five byte-identical CSVs and nine commands each",
        "protected_observations_included": False, "observation_level_models_refitted": False,
        "review_author_anonymity_certified": False, "legal_attribution_preserved": True,
        "data_license_finalized": False, "public_upload_performed": False, "doi": None,
        "review_manifest": "review_aggregate_candidate/review_manifest.json",
        "internal_source_mapping": "internal_provenance.json",
    }
    (out / "release_candidate_manifest.json").write_bytes(json_bytes(report))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--source-aggregate-zip", required=True, type=Path)
    parser.add_argument("--robustness-zip", required=True, type=Path)
    parser.add_argument("--verified-assets-zip", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args), indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Rebuild frozen Chapter 1 numerical analyses from archived source ZIPs.

The runner deliberately accepts local ZIP files rather than downloading GitHub
Actions artifacts itself. This keeps the reconstruction path usable after
Actions retention expires and makes the exact archived bytes part of the input
contract.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = {
    "continuous_zip": "101e996b638996a0c5ae79d358bf51293c3585f0e84c4a961b91dcbedf96211e",
    "multilevel_zip": "51e7a26b5bd09e030b67b9342586699abaaf46e630f45b6bb4ee7bfc9152ced6",
    "historical_zip": "499061e7a49f9455cf8c367fe26e313b7e0e33b2280d2354717e61a90ea8c6bc",
    "traits": "d775794f2bce2dfd0c1f63c5c8e01778c518f6eeb327bf0d9944045143a02344",
    "core_environment": "2172e3570f684770d0f919ecd81265c8460574e287bc4fb057db4f719cab7bb0",
    "reference_environment": "1ab84254a80493776b4c435152ed3d2a1c1e68dd0e0342da0ea081eeb5cd3d9b",
    "full_environment": "e242aa7ce69d12b11937c1335e84b9638799c50b42ef36b95725e77190df98e7",
    "regions": "085c4e8d45ceb34d32c6c961675ce74a4f0a33580f6cdd8ecd2ff1800a6364ff",
    "native_status": "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a",
    "historical_s1": "8ef5d5ea5f4e0c2f166071244a838cae77a0fe582817d729bba0b36f6b5ccd92",
    "historical_s2_random": "82655f79297e44a6630a599d8b0a1dc6f85e792812b8b01f2234e567a478e3af",
    "historical_s3": "8ef5d5ea5f4e0c2f166071244a838cae77a0fe582817d729bba0b36f6b5ccd92",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--continuous-zip", type=Path, required=True)
    parser.add_argument("--multilevel-zip", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=("validate_sources", "rebuild_full27", "rebuild_secondary", "rebuild_all"),
        default="validate_sources",
    )
    parser.add_argument("--raster-cache-dir", type=Path)
    parser.add_argument("--with-spatial", action="store_true")
    parser.add_argument("--with-sampling", action="store_true")
    parser.add_argument("--regions", type=Path)
    parser.add_argument("--native-status", type=Path)
    parser.add_argument("--with-historical", action="store_true")
    historical = parser.add_mutually_exclusive_group()
    historical.add_argument("--historical-zip", type=Path)
    historical.add_argument("--tree-dir", type=Path)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256(path)
    if actual != expected:
        raise SystemExit(f"{label} SHA-256 mismatch: {actual}; expected {expected}")


def safe_extract(zip_path: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if root not in target.parents and target != root:
                raise SystemExit(f"Unsafe ZIP path: {member.filename}")
        archive.extractall(destination)


def run(*parts: object) -> None:
    command = [str(part) for part in parts]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def verify_sources(continuous: Path, multilevel: Path) -> dict[str, object]:
    traits = continuous / "universe" / "continuous_trait_universe_observation_long.csv"
    core_environment = continuous / "environment" / "strict_spatial_chelsa.csv"
    reference_environment = multilevel / "process_environment" / "complete18_chelsa_process.csv"
    for path in (traits, core_environment, reference_environment):
        if not path.is_file():
            raise SystemExit(f"Required archived input is missing: {path}")
    require_hash(traits, EXPECTED["traits"], "trait universe")
    require_hash(core_environment, EXPECTED["core_environment"], "strict-spatial core environment")
    require_hash(reference_environment, EXPECTED["reference_environment"], "complete-18 reference environment")
    trait_frame = pd.read_csv(
        traits,
        usecols=["obs_id", "taxon_name", "endpoint_id", "value"],
        low_memory=False,
    )
    core = pd.read_csv(core_environment, low_memory=False)
    reference = pd.read_csv(reference_environment, low_memory=False)
    checks = {
        "trait_rows": int(len(trait_frame)),
        "registered_endpoint_ids": int(trait_frame["endpoint_id"].nunique()),
        "core_environment_rows": int(len(core)),
        "core_environment_taxa": int(core["taxon_name"].nunique()),
        "reference_environment_rows": int(len(reference)),
        "reference_environment_taxa": int(reference["taxon_name"].nunique()),
    }
    expected_checks = {
        "registered_endpoint_ids": 27,
        "core_environment_rows": 46276,
        "core_environment_taxa": 259,
        "reference_environment_rows": 1874,
        "reference_environment_taxa": 124,
    }
    for key, expected in expected_checks.items():
        if checks[key] != expected:
            raise SystemExit(f"Unexpected {key}: {checks[key]}; expected {expected}")
    return checks


def resolve_tree_dir(args: argparse.Namespace) -> Path | None:
    if not args.with_historical:
        return None
    if args.historical_zip is None and args.tree_dir is None:
        raise SystemExit("--with-historical requires --historical-zip or --tree-dir")
    if args.historical_zip is not None:
        require_hash(args.historical_zip, EXPECTED["historical_zip"], "historical tree source ZIP")
        extracted = args.out_dir / "source" / "historical"
        safe_extract(args.historical_zip, extracted)
        tree_dir = extracted / "historical_trees"
    else:
        tree_dir = args.tree_dir
    assert tree_dir is not None
    expected_files = {
        "gbotb_lcvp_scenario1.tre": EXPECTED["historical_s1"],
        "gbotb_lcvp_scenario2_randomized.trees": EXPECTED["historical_s2_random"],
        "gbotb_lcvp_scenario3.tre": EXPECTED["historical_s3"],
    }
    for name, expected in expected_files.items():
        path = tree_dir / name
        if not path.is_file():
            raise SystemExit(f"Required historical tree resource is missing: {path}")
        require_hash(path, expected, name)
    return tree_dir


def reconstruct_environment(core_environment: Path, out_dir: Path, raster_cache_dir: Path | None) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / "strict_spatial_chelsa_process.csv"
    report = out_dir / "chelsa_process_extraction_report.json"
    command: list[object] = [
        sys.executable,
        "analysis/sample_chelsa_process_environment.py",
        "--environment",
        core_environment,
        "--source-contract",
        "analysis/ch1/chelsa_process_environment_sources.json",
        "--out-csv",
        output,
        "--report",
        report,
    ]
    if raster_cache_dir is not None:
        command.extend(["--raster-cache-dir", raster_cache_dir])
    run(*command)
    require_hash(output, EXPECTED["full_environment"], "reconstructed nine-predictor environment")
    return output


def rebuild_full27(traits: Path, environment: Path, reference_environment: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    run(
        sys.executable,
        "analysis/run_geb_v2_full27_environment_atlas.py",
        "--contract",
        "ch1_global/v2/ontology/ch1_continuous_trait_contract.csv",
        "--analysis-contract",
        "analysis/ch1/v2_full27_environment_atlas_contract.json",
        "--traits-long",
        traits,
        "--environment",
        environment,
        "--out-dir",
        out_dir,
    )
    validation = out_dir / "v2_full27_environment_validation.json"
    run(
        sys.executable,
        "analysis/validate_geb_v2_full27_environment_atlas.py",
        "--contract",
        "ch1_global/v2/ontology/ch1_continuous_trait_contract.csv",
        "--traits-long",
        traits,
        "--environment",
        environment,
        "--reference-environment",
        reference_environment,
        "--results-dir",
        out_dir,
        "--report",
        validation,
    )
    payload = json.loads(validation.read_text(encoding="utf-8"))
    if payload.get("status") != "PASS" or payload.get("n_checks") != payload.get("n_passed"):
        raise SystemExit("Rebuilt full-27 atlas failed fail-closed validation")


def rebuild_secondary(traits: Path, environment: Path, out_dir: Path) -> None:
    run(
        sys.executable,
        "analysis/run_capitulum_functional_space.py",
        "--traits-long",
        traits,
        "--out-dir",
        out_dir / "capitulum_space",
    )
    run(
        sys.executable,
        "analysis/run_capitulum_environment_blocks.py",
        "--traits-long",
        traits,
        "--environment",
        environment,
        "--block-contract",
        "analysis/ch1/capitulum_environment_blocks_contract.json",
        "--out-dir",
        out_dir / "capitulum_environment",
    )
    run(
        sys.executable,
        "analysis/run_capitulum_environment_incremental.py",
        "--traits-long",
        traits,
        "--environment",
        environment,
        "--block-contract",
        "analysis/ch1/capitulum_environment_blocks_contract.json",
        "--out-dir",
        out_dir / "capitulum_environment_incremental",
    )


def rebuild_optional_sensitivities(
    args: argparse.Namespace,
    traits: Path,
    environment: Path,
    atlas_dir: Path,
    tree_dir: Path | None,
) -> None:
    spatial_dir = args.out_dir / "spatial"
    if args.with_sampling:
        if args.regions is None or args.native_status is None:
            raise SystemExit("--with-sampling requires --regions and --native-status")
        require_hash(args.regions, EXPECTED["regions"], "frozen broad-region lookup")
        require_hash(args.native_status, EXPECTED["native_status"], "frozen native-status table")
        run(
            sys.executable,
            "analysis/run_geb_v2_full27_sampling_composition_sensitivity.py",
            "--analysis-contract",
            "analysis/ch1/v2_full27_environment_atlas_contract.json",
            "--sensitivity-contract",
            "analysis/ch1/v2_full27_sampling_composition_sensitivity_contract.json",
            "--traits-long",
            traits,
            "--environment",
            environment,
            "--atlas-within",
            atlas_dir / "v2_full27_environment_within.csv",
            "--atlas-among",
            atlas_dir / "v2_full27_environment_among.csv",
            "--regions",
            args.regions,
            "--native-status",
            args.native_status,
            "--out-dir",
            args.out_dir / "sampling",
        )
    if args.with_spatial or args.with_historical:
        run(
            sys.executable,
            "analysis/run_geb_v2_full27_spatial_sensitivity.py",
            "--analysis-contract",
            "analysis/ch1/v2_full27_environment_atlas_contract.json",
            "--traits-long",
            traits,
            "--environment",
            environment,
            "--atlas-dir",
            atlas_dir,
            "--out-dir",
            spatial_dir,
        )
    if args.with_historical:
        assert tree_dir is not None
        run(
            sys.executable,
            "analysis/run_geb_v2_full27_historical_sensitivity.py",
            "--analysis-contract",
            "analysis/ch1/v2_full27_environment_atlas_contract.json",
            "--traits-long",
            traits,
            "--environment",
            environment,
            "--spatial-dir",
            spatial_dir,
            "--tree-dir",
            tree_dir,
            "--out-dir",
            args.out_dir / "historical",
        )


def main() -> int:
    args = parse_args()
    require_hash(args.continuous_zip, EXPECTED["continuous_zip"], "continuous source ZIP")
    require_hash(args.multilevel_zip, EXPECTED["multilevel_zip"], "multilevel source ZIP")
    source = args.out_dir / "source"
    continuous = source / "continuous"
    multilevel = source / "multilevel"
    safe_extract(args.continuous_zip, continuous)
    safe_extract(args.multilevel_zip, multilevel)
    checks = verify_sources(continuous, multilevel)
    tree_dir = resolve_tree_dir(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "source_validation.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                **checks,
                "historical_tree_source_verified": bool(tree_dir is not None),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    reference_environment = multilevel / "process_environment" / "complete18_chelsa_process.csv"
    run(
        sys.executable,
        "analysis/audit_v2_environment_collinearity.py",
        "--environment",
        reference_environment,
        "--out-dir",
        args.out_dir / "environment_collinearity",
    )
    if args.mode == "validate_sources":
        print(json.dumps({"status": "PASS", "mode": args.mode, **checks}, indent=2))
        return 0

    traits = continuous / "universe" / "continuous_trait_universe_observation_long.csv"
    core_environment = continuous / "environment" / "strict_spatial_chelsa.csv"
    environment = reconstruct_environment(core_environment, args.out_dir / "environment", args.raster_cache_dir)
    atlas_dir = args.out_dir / "full27"
    if args.mode in {"rebuild_full27", "rebuild_all"} or args.with_sampling or args.with_spatial or args.with_historical:
        rebuild_full27(traits, environment, reference_environment, atlas_dir)
    if args.mode in {"rebuild_secondary", "rebuild_all"}:
        rebuild_secondary(traits, environment, args.out_dir / "secondary")
    if args.with_sampling or args.with_spatial or args.with_historical:
        rebuild_optional_sensitivities(args, traits, environment, atlas_dir, tree_dir)
    print(json.dumps({"status": "PASS", "mode": args.mode, "output": str(args.out_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

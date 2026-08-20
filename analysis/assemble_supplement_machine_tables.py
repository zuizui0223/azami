#!/usr/bin/env python3
"""Assemble machine-readable Supplement tables from pinned frozen workflow artifacts.

This is a packaging step only. It copies exact archived output tables, renames them
for Supplement presentation and writes a checksum/provenance manifest. It never
refits a model or changes a scientific result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spde-root", type=Path, required=True)
    p.add_argument("--precision-root", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--spde-artifact-id", required=True)
    p.add_argument("--precision-artifact-id", required=True)
    p.add_argument("--git-sha", default="")
    return p.parse_args()


def unique(root: Path, name: str) -> Path:
    hits = sorted(p for p in root.rglob(name) if p.is_file())
    if len(hits) != 1:
        raise SystemExit(f"Expected exactly one {name} under {root}, found {len(hits)}")
    return hits[0]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy(src: Path, dst: Path) -> dict[str, object]:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return {
        "source_name": src.name,
        "output": str(dst),
        "size_bytes": dst.stat().st_size,
        "sha256": sha256(dst),
    }


def main() -> None:
    a = parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    spde_map = {
        "spde_fixed_effects_by_group.csv": "tables/S06_spde_fixed_effects_full.csv",
        "spde_effect_stability_across_groups.csv": "tables/S06b_spde_effect_stability.csv",
        "spde_hyperparameters_by_group.csv": "tables/S06c_spde_hyperparameters.csv",
        "spde_predictor_selection_by_group.csv": "tables/S06d_spde_predictor_selection.csv",
        "spde_run_metadata.csv": "tables/S06e_spde_run_metadata.csv",
    }
    precision_map = {
        "table_s2_species_specific_coefficients.csv": "release_only/S11_species_specific_coefficients_archived.csv",
        "reviewer_precision_summary.json": "release_only/S11_reviewer_precision_summary.json",
    }

    files: list[dict[str, object]] = []
    for src_name, rel in spde_map.items():
        files.append(copy(unique(a.spde_root, src_name), a.out_dir / rel))
    for src_name, rel in precision_map.items():
        files.append(copy(unique(a.precision_root, src_name), a.out_dir / rel))

    manifest = {
        "purpose": "machine-readable Supplement packaging from frozen artifacts",
        "scientific_effect": "none; exact archived outputs copied without refitting",
        "spde_artifact_id": str(a.spde_artifact_id),
        "precision_artifact_id": str(a.precision_artifact_id),
        "git_sha": a.git_sha,
        "files": files,
    }
    target = a.out_dir / "SUPPLEMENT_MACHINE_TABLE_PROVENANCE.json"
    target.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

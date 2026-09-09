#!/usr/bin/env python3
"""Run the existing v2 full-27 atlas on a native-only v2 cohort.

This is intentionally a thin wrapper around analysis/run_geb_v2_full27_environment_atlas.py.
It changes only cohort membership and restores five fields already measured in the historical
head table but omitted from the old observation aggregation. It does not change the nine
predictors, inferential statistics, BH correction, endpoint contract, or v2 baseline outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FROZEN_NATIVE_SHA = "c01eeb9ff245d7f73da1a12fa4eede904dd9770467655f20e3d85de2ac8dd84a"
REGENERATED_NATIVE_SHA = "9686b8f515deef3b3aa9311b5137317a094af195e0f2414f2b1b70d7c72b5021"
EXPECTED_NATIVE_ROWS = 27066
EXPECTED_TOTAL_ROWS = 46276
EXPECTED_STATUS_COUNTS = {
    "native": 27066,
    "introduced": 10554,
    "unresolved_taxon": 5491,
    "unmapped_tdwg": 2100,
    "unlisted": 1065,
}
RECOVERED = {
    "visible_floret_fraction": "corolla_visible_fraction",
    "corolla_white_pixel_fraction": "corolla_white_fraction",
    "corolla_redmagenta_pixel_fraction": "corolla_redmagenta_fraction",
    "corolla_purple_pixel_fraction": "corolla_purple_fraction",
    "corolla_yellow_pixel_fraction": "corolla_yellow_fraction",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path, dict]:
    native_sha = sha256(args.native_status)
    if native_sha == FROZEN_NATIVE_SHA:
        native_source_mode = "exact_frozen_lfs_csv"
    elif native_sha == REGENERATED_NATIVE_SHA:
        native_source_mode = "deterministic_regenerated_from_frozen_v2_auxiliary_tables"
    else:
        raise SystemExit(f"Unrecognized native-status input hash: {native_sha}")

    native_all = pd.read_csv(
        args.native_status,
        usecols=["obs_id", "taxon_name", "native_range_status"],
        low_memory=False,
    )
    native_all["obs_id"] = native_all["obs_id"].astype(str)
    if len(native_all) != EXPECTED_TOTAL_ROWS or native_all["obs_id"].duplicated().any():
        raise SystemExit("Native-status table must contain 46,276 unique frozen observation IDs")
    status_counts = {
        str(k): int(v) for k, v in native_all["native_range_status"].value_counts().items()
    }
    if status_counts != EXPECTED_STATUS_COUNTS:
        raise SystemExit(f"Native-status counts differ from frozen v2: {status_counts}")

    environment = pd.read_csv(args.environment, low_memory=False)
    environment["obs_id"] = environment["obs_id"].astype(str)
    if len(environment) != EXPECTED_TOTAL_ROWS or environment["obs_id"].duplicated().any():
        raise SystemExit("Environment input is not the frozen 46,276-row universe")
    if set(environment["obs_id"]) != set(native_all["obs_id"]):
        raise SystemExit("Native-status IDs and frozen environment IDs do not match")

    native = native_all[native_all["native_range_status"].eq("native")].copy()
    native_ids = set(native["obs_id"])
    if len(native) != EXPECTED_NATIVE_ROWS:
        raise SystemExit(f"Expected {EXPECTED_NATIVE_ROWS} native observations, found {len(native)}")

    environment = environment[environment["obs_id"].isin(native_ids)].copy()
    if len(environment) != EXPECTED_NATIVE_ROWS:
        raise SystemExit("Native environment membership differs from the frozen v2 universe")

    traits = pd.read_csv(args.traits_long, low_memory=False)
    traits["obs_id"] = traits["obs_id"].astype(str)
    traits = traits[traits["obs_id"].isin(native_ids)].copy()

    recovered = pd.read_csv(args.recovered_display, low_memory=False)
    recovered["obs_id"] = recovered["obs_id"].astype(str)
    recovered = recovered[recovered["obs_id"].isin(native_ids)].set_index("obs_id")
    if recovered.index.duplicated().any():
        raise SystemExit("Recovered display table must be unique by obs_id")

    for endpoint, column in RECOVERED.items():
        values = pd.to_numeric(recovered[column], errors="coerce")
        mapping = values.dropna().to_dict()
        mask = traits["endpoint_id"].eq(endpoint)
        mapped = traits.loc[mask, "obs_id"].map(mapping)
        usable = mapped.notna()
        target = traits.loc[mask].index[usable]
        traits.loc[target, "value"] = mapped[usable].to_numpy(float)
        traits.loc[target, "measurement_available"] = True

    args.work_dir.mkdir(parents=True, exist_ok=True)
    trait_out = args.work_dir / "native_traits_long.csv"
    env_out = args.work_dir / "native_environment.csv"
    contract_out = args.work_dir / "native_analysis_contract.json"
    traits.to_csv(trait_out, index=False)
    environment.to_csv(env_out, index=False)

    contract = json.loads(args.analysis_contract.read_text(encoding="utf-8"))
    contract["analysis_id"] = "geb_v2_full27_native_primary_v1"
    contract["status"] = "retrospective_native_primary_reanalysis_reusing_frozen_v2_model"
    contract["purpose"] = "Run the unchanged v2 full-27 x nine-predictor atlas on the native-only v2 cohort; restore five already-measured fields omitted by the historical aggregator."
    contract["environment"]["cohort"] = "27066 native observations selected from the frozen 46276 strict-spatial v2 universe"
    contract["environment"]["expected_observations"] = EXPECTED_NATIVE_ROWS
    contract["claim_boundary"] = list(contract.get("claim_boundary", [])) + [
        "Native-only membership is a scope change, not a new acquisition universe.",
        "The five restored display/composition endpoints were already measured at head level and are not new image measurements.",
        "When the native-status source mode is regenerated, the original historical LFS bytes are unavailable and no byte-identity claim is made.",
    ]
    contract_out.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

    counts = {
        "native_source_mode": native_source_mode,
        "native_status_sha256": native_sha,
        "native_observations": len(environment),
        "native_taxa_environment": int(environment["taxon_name"].nunique()),
        "trait_rows_native": int(len(traits)),
        "measured_endpoint_counts": {
            endpoint: int(
                pd.to_numeric(
                    traits.loc[traits["endpoint_id"].eq(endpoint), "value"],
                    errors="coerce",
                ).notna().sum()
            )
            for endpoint in RECOVERED
        },
    }
    return trait_out, env_out, contract_out, counts


def signal_summary(path: Path, scope: str | None = None) -> dict:
    frame = pd.read_csv(path, low_memory=False)
    if scope is not None:
        frame = frame[frame["scope"].eq(scope)].copy()
    ok = frame["status"].eq("ok")
    q = pd.to_numeric(frame["q_fdr_bh_global_family"], errors="coerce")
    sig = frame[ok & q.lt(0.05)]
    return {
        "ok_rows": int(ok.sum()),
        "fdr_signals": int(len(sig)),
        "signals_by_module": {
            str(k): int(v) for k, v in sig["module"].value_counts().sort_index().items()
        },
    }


def _joined_column(frame: pd.DataFrame, candidates: tuple[str, ...], suffix: str) -> str:
    for base in candidates:
        column = base + suffix
        if column in frame.columns:
            return column
    raise KeyError(f"None of {candidates} found with suffix {suffix}")


def compare_direction(new_path: Path, old_path: Path, scope: str | None = None) -> dict:
    new = pd.read_csv(new_path, low_memory=False)
    old = pd.read_csv(old_path, low_memory=False)
    if scope is not None:
        new = new[new["scope"].eq(scope)].copy()
        old = old[old["scope"].eq(scope)].copy()
    keys = ["unit_id", "predictor", "inferential_unit"]
    both = new.merge(old, on=keys, suffixes=("_new", "_old"))
    both = both[both["status_new"].eq("ok") & both["status_old"].eq("ok")].copy()
    stable = []
    for row in both.to_dict("records"):
        if row["inferential_unit"] == "linear_endpoint":
            a = float(row[_joined_column(both, ("beta_std", "beta_std_among"), "_new")])
            b = float(row[_joined_column(both, ("beta_std", "beta_std_among"), "_old")])
            stable.append(bool(a != 0 and b != 0 and np.sign(a) == np.sign(b)))
        else:
            sin_new = _joined_column(both, ("beta_sine_std", "beta_sine_std_among"), "_new")
            cos_new = _joined_column(both, ("beta_cosine_std", "beta_cosine_std_among"), "_new")
            sin_old = _joined_column(both, ("beta_sine_std", "beta_sine_std_among"), "_old")
            cos_old = _joined_column(both, ("beta_cosine_std", "beta_cosine_std_among"), "_old")
            an = np.array([float(row[sin_new]), float(row[cos_new])])
            ao = np.array([float(row[sin_old]), float(row[cos_old])])
            denom = float(np.linalg.norm(an) * np.linalg.norm(ao))
            stable.append(bool(denom > 0 and float(np.dot(an, ao) / denom) > 0))
    return {
        "comparable_ok_rows": int(len(both)),
        "same_direction_or_positive_circular_alignment": int(sum(stable)),
        "fraction_direction_stable": float(np.mean(stable)) if stable else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--endpoint-contract", type=Path, required=True)
    parser.add_argument("--analysis-contract", type=Path, required=True)
    parser.add_argument("--traits-long", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--native-status", type=Path, required=True)
    parser.add_argument("--recovered-display", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    traits, environment, contract, counts = build_inputs(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "analysis/run_geb_v2_full27_environment_atlas.py"),
            "--contract",
            str(args.endpoint_contract),
            "--analysis-contract",
            str(contract),
            "--traits-long",
            str(traits),
            "--environment",
            str(environment),
            "--out-dir",
            str(args.out_dir),
        ],
        cwd=ROOT,
        check=True,
    )

    old = ROOT / "analysis_outputs/v2_full27_environment_atlas_2026-08-27"
    report = {
        "status": "NATIVE_PRIMARY_V2_REUSE_REANALYSIS_COMPLETE",
        "scope": counts,
        "within": signal_summary(args.out_dir / "v2_full27_environment_within.csv"),
        "among_min5": signal_summary(
            args.out_dir / "v2_full27_environment_among.csv", "among_taxon_min5"
        ),
        "among_min2": signal_summary(
            args.out_dir / "v2_full27_environment_among.csv", "among_taxon_min2"
        ),
        "direction_vs_full46276": {
            "within": compare_direction(
                args.out_dir / "v2_full27_environment_within.csv",
                old / "v2_full27_environment_within.csv",
            ),
            "among_min5": compare_direction(
                args.out_dir / "v2_full27_environment_among.csv",
                old / "v2_full27_environment_among.csv",
                "among_taxon_min5",
            ),
        },
        "unchanged": [
            "27-endpoint contract",
            "nine environmental predictors",
            "within/among estimators",
            "9999 permutations",
            "BH family rule",
        ],
        "changed": [
            "primary cohort restricted to native observations",
            "five historically omitted measured fields restored to the trait table",
        ],
    }
    (args.out_dir / "native_primary_comparison_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

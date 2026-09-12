#!/usr/bin/env python3
"""Prepare an outcome-blind WCVP accepted-taxon sensitivity cohort.

The primary Chapter 1 analyses intentionally use source-assigned taxon labels.
This sensitivity addresses the named taxonomic-lumping/synonym threat without
changing that primary definition. Observations are retained only when the
frozen WCVP audit resolved the source canonical name to exactly one accepted
key; all retained taxon labels are then replaced by that accepted key, so
source synonyms collapse before the existing v3 analyses are rerun.

No trait or environmental outcome is used to decide retention or grouping.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path

import pandas as pd

from reproducibility.run_current_analysis import NATIVE_SHA, native_bytes

RESOLVED = "resolved_unique_accepted_key"


def canonical_key(value) -> str:
    if pd.isna(value):
        raise ValueError("resolved WCVP row has missing accepted_key")
    try:
        return f"wcvp:{int(float(value))}"
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid accepted_key {value!r}") from exc


def prepare(traits_path: Path, environment_path: Path, native_status_path: Path, out_dir: Path) -> dict:
    normalized = native_bytes(native_status_path.read_bytes())
    # native_bytes itself fails unless one of the permitted newline forms has
    # the frozen hash; preserve the identity explicitly in the report too.
    import hashlib
    actual_native_sha = hashlib.sha256(normalized).hexdigest()
    if actual_native_sha != NATIVE_SHA:
        raise ValueError(f"native status hash mismatch: {actual_native_sha} != {NATIVE_SHA}")

    native = pd.read_csv(io.BytesIO(normalized), dtype={"obs_id": str}, low_memory=False)
    required_native = {"obs_id", "taxon_name", "resolution_status", "accepted_key", "accepted_name"}
    missing = required_native.difference(native.columns)
    if missing:
        raise ValueError(f"native-status table missing columns: {sorted(missing)}")
    if native.obs_id.duplicated().any():
        raise ValueError("native-status table must contain one row per observation")

    traits = pd.read_csv(traits_path, dtype={"obs_id": str}, low_memory=False)
    env = pd.read_csv(environment_path, dtype={"obs_id": str}, low_memory=False)
    for name, frame in (("traits", traits), ("environment", env)):
        if not {"obs_id", "taxon_name"}.issubset(frame.columns):
            raise ValueError(f"{name} input lacks obs_id/taxon_name")
        frame["taxon_name"] = frame.taxon_name.astype(str)
    if env.obs_id.duplicated().any():
        raise ValueError("environment input must contain one row per observation")

    map_cols = native[["obs_id", "taxon_name", "resolution_status", "accepted_key", "accepted_name"]].copy()
    map_cols = map_cols.rename(columns={"taxon_name": "source_taxon_name"})

    # Every environment observation should have the same source label in the
    # frozen authority table. This is a provenance check, not a filtering rule.
    check = env[["obs_id", "taxon_name"]].merge(map_cols, on="obs_id", how="left", validate="one_to_one")
    if check.source_taxon_name.isna().any():
        raise ValueError("environment observations are missing from frozen native-status table")
    mismatch = check.taxon_name.astype(str) != check.source_taxon_name.astype(str)
    if mismatch.any():
        raise ValueError(f"source taxon mismatch for {int(mismatch.sum())} environment observations")

    resolved = map_cols[map_cols.resolution_status == RESOLVED].copy()
    if resolved.empty:
        raise ValueError("no authority-resolved observations")
    resolved["authority_taxon"] = resolved.accepted_key.map(canonical_key)
    resolved_ids = set(resolved.obs_id)

    env_out = env[env.obs_id.isin(resolved_ids)].copy()
    env_out = env_out.drop(columns=["taxon_name"]).merge(
        resolved[["obs_id", "authority_taxon"]], on="obs_id", how="inner", validate="one_to_one"
    )
    env_out = env_out.rename(columns={"authority_taxon": "taxon_name"})

    traits_out = traits[traits.obs_id.isin(resolved_ids)].copy()
    traits_out = traits_out.drop(columns=["taxon_name"]).merge(
        resolved[["obs_id", "authority_taxon"]], on="obs_id", how="inner", validate="many_to_one"
    )
    traits_out = traits_out.rename(columns={"authority_taxon": "taxon_name"})

    # One source canonical name must not resolve to multiple accepted keys.
    resolved_taxon = resolved.groupby("source_taxon_name").agg(
        accepted_keys=("authority_taxon", lambda s: sorted(set(s))),
        accepted_names=("accepted_name", lambda s: sorted({str(v) for v in s if pd.notna(v) and str(v)})),
        n_observations=("obs_id", "nunique"),
    ).reset_index()
    bad = resolved_taxon.accepted_keys.map(len) != 1
    if bad.any():
        raise ValueError("a source taxon maps to multiple accepted WCVP keys")
    resolved_taxon["authority_taxon"] = resolved_taxon.accepted_keys.map(lambda x: x[0])
    resolved_taxon["accepted_name"] = resolved_taxon.accepted_names.map(lambda x: x[0] if len(x) == 1 else " | ".join(x))
    resolved_taxon = resolved_taxon.drop(columns=["accepted_keys", "accepted_names"])

    authority_groups = resolved_taxon.groupby("authority_taxon").agg(
        n_source_taxa=("source_taxon_name", "nunique"),
        source_taxa=("source_taxon_name", lambda s: " | ".join(sorted(set(s)))),
        accepted_name=("accepted_name", lambda s: " | ".join(sorted(set(s)))),
        n_observations=("n_observations", "sum"),
    ).reset_index()

    unresolved_taxa = sorted(set(native.taxon_name.astype(str)) - set(resolved_taxon.source_taxon_name.astype(str)))
    status_counts = native.resolution_status.fillna("missing").value_counts().sort_index().to_dict()

    out_dir.mkdir(parents=True, exist_ok=True)
    traits_out.to_csv(out_dir / "traits_authority.csv", index=False)
    env_out.to_csv(out_dir / "environment_authority.csv", index=False)
    resolved_taxon.to_csv(out_dir / "source_to_authority_taxon.csv", index=False)
    authority_groups.to_csv(out_dir / "authority_taxon_groups.csv", index=False)

    report = {
        "analysis_id": "ch1_v3_wcvp_authority_taxonomy_sensitivity_20260912",
        "selection_rule": (
            "retain only frozen native-status rows with resolution_status=resolved_unique_accepted_key; "
            "replace source taxon with accepted WCVP key; collapse source labels sharing a key; no outcome-based filtering"
        ),
        "native_status_sha256": NATIVE_SHA,
        "environment_observations_start": int(env.obs_id.nunique()),
        "environment_observations_resolved": int(env_out.obs_id.nunique()),
        "source_taxa_start": int(env.taxon_name.nunique()),
        "source_taxa_resolved": int(resolved_taxon.source_taxon_name.nunique()),
        "authority_taxa_after_collapse": int(authority_groups.authority_taxon.nunique()),
        "authority_taxa_with_multiple_source_labels": int((authority_groups.n_source_taxa > 1).sum()),
        "unresolved_source_taxa": unresolved_taxa,
        "resolution_status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "scientific_outputs_changed": False,
        "claim_boundary": (
            "taxonomy sensitivity only; primary source-assigned taxon analyses remain canonical unless manuscript scope is explicitly revised"
        ),
    }
    (out_dir / "taxonomy_preparation_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--traits", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--native-status", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.traits, args.environment, args.native_status, args.out_dir), indent=2))


if __name__ == "__main__":
    main()

"""Unified native Capitulum Gate A and Gate E: exact-source, private-row QC and concept joins.

The native source is selected BEFORE trait eligibility. All private identifiers stay
inside the Actions runner. Public output contains counts and the seven already
published East Asian taxon concepts, never observation/observer/cell identifiers.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration_audit import eligible
from analysis.v4.native_integration import FEATURES
from analysis.v4.native_spatial import cells

EXPECTED_N = 313350
EXPECTED_TAXA = 306
CORE9 = tuple(x for parts in FEATURES.values() for x in parts)
EXPECTED_CORE9 = (157779, 242)
MAPPING_GIT_BLOB_SHA1 = "b662f9687f11d4069965531d5214962b07f6df8c"
EA_FOCAL = (
    "Cirsium japonicum var. albescens", "Cirsium japonicum var. australe",
    "Cirsium brevicaule", "Cirsium irumtiense", "Cirsium kawakamii",
    "Cirsium pengii", "Cirsium tatakaense",
)


def authority_mapping(path: Path) -> list[dict]:
    raw = path.read_bytes()
    # The previously executed phylogenetic mapping is a fixed observation-derived
    # concept map; it must not be silently replaced by a fresh name search.
    blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if blob != MAPPING_GIT_BLOB_SHA1:
        raise ValueError("Tree/taxon concept source changed: stop for expert reconciliation")
    rows = json.loads(raw)
    if len(rows) != EXPECTED_TAXA or len({str(r["accepted_key"]) for r in rows}) != EXPECTED_TAXA:
        raise ValueError("Original WCVP accepted concepts changed")
    return rows


def source_and_observer(frame: pd.DataFrame, observers: pd.DataFrame | None):
    if len(frame) != EXPECTED_N or frame["obs_id"].isna().any() or not frame["obs_id"].is_unique:
        raise ValueError("Native roster count or observation keys changed")
    frame = frame.copy()
    frame["obs_id"] = frame["obs_id"].astype(str)
    frame["accepted_key"] = frame["accepted_key"].astype(str)
    if frame["accepted_key"].nunique() != EXPECTED_TAXA:
        raise ValueError("Native taxon count changed")
    if observers is None:
        return frame, False
    required = {"obs_id", "user_id"}
    if not required <= set(observers):
        raise ValueError("Original photographer join unavailable")
    o = observers.copy()
    o["obs_id"] = o["obs_id"].astype(str)
    o["user_id"] = o["user_id"].astype("string")
    if not o.obs_id.is_unique or o.user_id.isna().any() or len(o) != EXPECTED_N:
        raise ValueError("Original native observer roster incomplete or duplicated")
    if set(o.obs_id) != set(frame.obs_id):
        raise ValueError("Original observers do not exactly match native observation IDs")
    original_cols = ["obs_id", "user_id"]
    if "quality_grade" in o:
        original_cols.append("quality_grade")
    frame = frame.merge(o[original_cols], on="obs_id", how="left", validate="one_to_one")
    if frame.user_id.isna().any():
        raise ValueError("Missing observer after one-to-one source join")
    return frame, True


def _monthly(frame, mask):
    if "observed_on" not in frame:
        return {"status": "not_recorded"}
    month = pd.to_datetime(frame["observed_on"], errors="coerce").dt.month
    return {
        "status": "source_vs_matched_observation_month_not_anthesis",
        "observed_source": int(month.notna().sum()),
        "observed_matched": int(month[mask].notna().sum()),
        "monthly": [
            {"month": k, "source": int((month == k).sum()),
             "matched": int((month[mask] == k).sum())}
            for k in range(1, 13)
        ],
    }


def gate_a(
    frame: pd.DataFrame, mapping: list[dict], observers: pd.DataFrame | None,
    *, min_cells: int = 5, min_observers: int = 10, min_balanced: int = 50,
    max_taxa: int = 40
) -> tuple[dict, pd.DataFrame]:
    if min_cells < 2 or min_observers < 2 or min_balanced < 10 or max_taxa < 2:
        raise ValueError("Unreasonable spatial/observer/rareness sample thresholds")
    frame, have_observers = source_and_observer(frame, observers)
    expected_keys = {str(r["accepted_key"]) for r in mapping}
    if set(frame.accepted_key) != expected_keys:
        raise ValueError("Native accepted taxonomy differs from frozen mapping")
    core_mask = eligible(frame, list(CORE9)).to_numpy(bool)
    core = frame.loc[core_mask].copy()
    if (len(core), core.accepted_key.nunique()) != EXPECTED_CORE9:
        raise ValueError("Same-observation nine-coordinate cohort changed")
    collisions = {str(r["accepted_key"]): r for r in mapping if len(set(r["source_names"])) > 1}
    uncertain_tree = {str(r["accepted_key"]): r for r in mapping if r["status"] != "matched"}
    # A source-name collision does not become taxonomically valid just because
    # one of its source names happens to match a synthetic tree tip.
    no_collision = core.loc[~core.accepted_key.isin(collisions)].copy()
    no_collision["_cell"] = cells(
        no_collision.latitude.to_numpy(float), no_collision.longitude.to_numpy(float))
    # Identity-stable draw from within a species×cell×human-source stratum.
    # No image feature, ecological response or fitness outcome enters this choice.
    no_collision["_sort_hash"] = pd.util.hash_pandas_object(
        no_collision.obs_id, index=False).to_numpy("uint64")
    group = ["accepted_key", "_cell"] + (["user_id"] if have_observers else [])
    distinct = (no_collision.sort_values([*group, "_sort_hash", "obs_id"])
                .drop_duplicates(group, keep="first"))
    groups = distinct.groupby("accepted_key", sort=True).agg(
        records=("obs_id", "size"),
        independent_cells=("_cell", "nunique"),
    )
    if have_observers:
        groups["independent_observers"] = distinct.groupby("accepted_key").user_id.nunique()
    else:
        groups["independent_observers"] = 0
    # Scientific inclusion is observer-aware: no observer identity means *no*
    # supposedly observer-balanced hypervolume or definitive cohort.
    qualifiers = groups.loc[
        (groups.records >= min_balanced)
        & (groups.independent_cells >= min_cells)
        & (groups.independent_observers >= min_observers)
    ].copy()
    qualifiers["_key"] = qualifiers.index.astype(str)
    qualifiers = qualifiers.sort_values(
        ["independent_cells", "independent_observers", "records", "_key"],
        ascending=[False, False, False, True]).head(max_taxa)
    retained = distinct.loc[distinct.accepted_key.isin(qualifiers.index)].copy()
    if have_observers:
        original_observers = int(frame.user_id.nunique())
        matched_observers = int(core.user_id.nunique())
    else:
        original_observers = matched_observers = None
    all_cells = no_collision.groupby("accepted_key")["_cell"].nunique()
    quality = {}
    if "quality_grade" in frame:
        a = frame.quality_grade.fillna("unknown").astype(str).value_counts()
        b = core.quality_grade.fillna("unknown").astype(str).value_counts()
        quality = {str(k): {"source": int(a[k]), "matched": int(b.get(k, 0))}
                   for k in a.index}
    report = {
        "status": "observer_cell_balanced_selection_constructed" if have_observers else "blocked_missing_original_observers",
        "source_native_observations": len(frame),
        "source_accepted_keys": EXPECTED_TAXA,
        "matched_nine_output_observations": len(core),
        "matched_nine_output_taxa": int(core.accepted_key.nunique()),
        "source_original_observers": original_observers,
        "matched_nine_output_original_observers": matched_observers,
        "repeated_name_accepted_key_collisions": len(collisions),
        "source_merged_name_concepts": int(sum(len(r["source_names"]) for r in collisions.values())),
        "matched_rows_in_collided_accepted_keys": int(core.accepted_key.isin(collisions).sum()),
        "matched_taxa_after_collision_exclusion": int(no_collision.accepted_key.nunique()),
        "tree_mapping_not_exactly_matched_keys": len(uncertain_tree),
        "spatial_method": "5-degree longitude by equal-sin(latitude) spherical cells, as used in primary Gate D spatial diagnostic",
        "taxa_with_min_5_occupied_cells": int((all_cells >= 5).sum()),
        "taxa_with_min_10_occupied_cells": int((all_cells >= 10).sum()),
        "taxa_with_min_20_occupied_cells": int((all_cells >= 20).sum()),
        "predeclared_selection": {"minimum_distinct_cells": min_cells,
             "minimum_original_observers": min_observers,
             "minimum_distinct_taxon_cell_observer_observations": min_balanced,
             "maximum_taxa": max_taxa,
             "selected_taxa": int(len(qualifiers)),
             "selected_nonredundant_observations": int(len(retained))},
        "observation_month": _monthly(frame, core_mask),
        "quality_grade": quality if quality else {"status": "unavailable"},
        "image_quality_and_camera_pose": "not_externally_validated",
        "not_interpretable_as": ["independent individual plants", "genetic variation",
                                 "source-corrected calibrated pigment/angle", "complete 27-trait head"],
    }
    return report, retained


def gate_e(mapping: list[dict], *, currently_selected_keys: set[str] | None = None) -> dict:
    def normalize(s): return " ".join(s.replace("_", " ").split())
    focal = []
    for name in EA_FOCAL:
        matches = [r for r in mapping if name in [normalize(t) for t in r["source_names"]]]
        if len(matches) > 1:
            raise ValueError("One source concept assigned multiple WCVP keys")
        if not matches:
            status, source, phylo = "not_in_native_source_as_exact_concept", None, None
        else:
            row = matches[0]
            source = str(row["accepted_key"])
            if len(row["source_names"]) > 1:
                status = "COLLISION_DISTINCT_TAXA_SHARE_ACCEPTED_KEY_DO_NOT_POOL"
            else:
                status = "single_source_name_exact"
            phylo = str(row["status"])
        focal.append({
            "published_concept": name, "status": status,
            "native_key_present": source is not None,
            "historical_tree_mapping": phylo,
            "selected_native_comparison": bool(
                source in currently_selected_keys) if currently_selected_keys is not None else None,
        })
    collisions = [r for r in mapping if len(r["source_names"]) > 1]
    return {
        "status": "exact_source_taxon_concepts_checked_not_functional_verification",
        "eazami_pilot_original_source_taxa": len(EA_FOCAL),
        "eazami_source_concepts": focal,
        "collision_keys": len(collisions),
        "collision_n_source_names": sum(len(r["source_names"]) for r in collisions),
        "ryukyu_pair_share_wcvp_key": (
            focal[2]["status"].startswith("COLLISION")
            and focal[3]["status"].startswith("COLLISION")
        ),
        "cross_repo_aza3_evidence": {
            "published_cirsium_arvense_pollinator_exclusion": "reproductive treatment, NOT phenotype-specific causal effect",
            "north_american_63_tip_hummingbird_syndrome": "syndrome coding is not independent observed pollination; exact concepts not silently joined",
        },
        "functional_gate": "NOT_CLOSED: image proxies require same-individual pigment, gravity, anatomy, pollination/herbivory and fitness validation",
        "phylogenetic_gate": "NOT_CLOSED: no new dated molecular tree, cytotype or independent lineage contrasts",
    }

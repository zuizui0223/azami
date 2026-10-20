"""Gate A: read-only native-first measurement intersection and support audit.

Source identity: the verified pre-detector native subset prepared from the
665,115-observation 2026-07-03 snapshot. No original image/raster download,
taxonomic update, phenotype remeasurement or ecological model is performed.

The public report contains ONLY aggregate counts, never observation, observer,
taxon, photo or locality identifiers. Scientific endpoint eligibility remains
distinct from genuine biological absence.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

PREPARED_SHA256 = "8ceed49e4ac73b31f9a5c726ca84f4198bfc3c51e9ee43303e29524e915a5c29"
# Environmental rasters in this ORIGINAL table had undecoded scale/offset and
# invalid GSP sentinel values. For any new joint Gate A-D execution, use the
# separately hash-verified PR143 corrected table; keep the original fingerprint
# exclusively as a fail-closed diagnostic for superseded inputs.
CORRECTED_PREPARED_SHA256 = "342d78e9a8a960e8d9bdcf0bc2cdacd13bee48214a8b962419d5a71d5aef46ed"
SOURCE_RECEIPT_SHA256 = "246a5ea9643617c6dc766e0b8d27a444a5496afe150d7768032740cf6a53f0a6"
NATIVE_N = 313350
MEASURED_N = 217104
HEAD_ARCHIVE_SHA256 = "b3d13ba7014100d4e44bf63d7732e21361c863db949f304e513be140a5011aa8"
HEAD_ARCHIVE_N = 1255791
CORE3 = (
    "orientation_image_vertical_angle",
    "corolla_lab_chroma",
    "capitulum_outline_aspect_ratio",
)
COLOUR4 = (
    "corolla_lab_lightness",
    "corolla_lab_chroma",
    "corolla_hue_sin",
    "corolla_hue_cos",
)
OUTLINE4 = (
    "capitulum_outline_aspect_ratio",
    "capitulum_outline_circularity",
    "capitulum_outline_solidity",
    "capitulum_width_profile_cv",
)
CORE9 = (CORE3[0], *COLOUR4, *OUTLINE4)
MODULES = {
    "presentation": (CORE3[0],),
    "colour": COLOUR4,
    "outline": OUTLINE4,
    "primary_core3": CORE3,
    "complete_primary_core9": CORE9,
}
RECOVERED_FIVE = {
    "visible_floret_fraction": "corolla_visible_fraction",
    "corolla_white_pixel_fraction": "corolla_white_fraction",
    "corolla_redmagenta_pixel_fraction": "corolla_redmagenta_fraction",
    "corolla_purple_pixel_fraction": "corolla_purple_fraction",
    "corolla_yellow_pixel_fraction": "corolla_yellow_fraction",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def load_ledger(path: Path) -> dict:
    ledger = json.loads(path.read_text(encoding="utf-8"))
    d = ledger["denominators"]
    if d["native_primary_observations"] != NATIVE_N:
        raise ValueError("Native roster denominator does not match the frozen source")
    if d["native_with_inherited_measurement_row"] != MEASURED_N:
        raise ValueError("Historical measured-row universe changed")
    if ledger["provenance"]["receipt_sha256"] != SOURCE_RECEIPT_SHA256:
        raise ValueError("Original measurement aggregate receipt identity changed")
    endpoints = ledger["endpoints"]
    if len(endpoints) != 27 or len({x["endpoint_id"] for x in endpoints}) != 27:
        raise ValueError("27-endpoint ontology changed")
    return ledger


def guaranteed_bounds(ledger: dict) -> dict:
    """Sharp Frechet bounds for intersecting subsets of the SAME measured-row universe."""
    rows = {r["endpoint_id"]: r for r in ledger["endpoints"]}
    universe = ledger["denominators"]["native_with_inherited_measurement_row"]
    sets = {
        "orientation_x_chroma": (CORE3[0], CORE3[1]),
        "orientation_x_outline_aspect": (CORE3[0], CORE3[2]),
        "chroma_x_outline_aspect": (CORE3[1], CORE3[2]),
        "matched_core3": CORE3,
        "complete_colour4": COLOUR4,
        "complete_outline4": OUTLINE4,
        "complete_core9": CORE9,
    }
    result = {}
    for key, endpoints in sets.items():
        if len(set(endpoints)) != len(endpoints):
            raise ValueError("Repeated endpoint within a mathematical intersection")
        vals = []
        for endpoint in endpoints:
            r = rows[endpoint]
            if r["source"] != "full_archived_observation_table":
                raise ValueError("Frechet bound cannot use distinct measurement universes")
            n = int(r["n_eligible_native"])
            if n > universe or n < 0:
                raise ValueError("Eligible count outside inherited measurement universe")
            vals.append(n)
        result[key] = {
            "endpoints": list(endpoints),
            "lower": max(0, sum(vals) - (len(vals) - 1) * universe),
            "upper": min(vals),
            "status": "aggregate_frechet_bound_not_exact_join",
        }
    return {
        "source_receipt_sha256": SOURCE_RECEIPT_SHA256,
        "native_source_observations": ledger["denominators"]["native_primary_observations"],
        "shared_inherited_measurement_universe": universe,
        "method": "sharp inclusion-exclusion limits, same inherited observation-row universe",
        "intersections": result,
        "five_historical_fields": {
            "frozen_27_endpoint_atlas_eligible_observations": 0,
            "historical_head_measurement_rows_expected": HEAD_ARCHIVE_N,
            "historical_head_source_sha256": HEAD_ARCHIVE_SHA256,
            "head_to_observation_recovery_code": "PR92 analysis/v3/recover_display_composition.py",
            "native_obs_recovery_status": "NOT_YET_VERIFIED",
            "interpretation": "five zero atlas columns do not prove absence of legacy head-level pixels",
        },
        "caution": "No exact native core3/core9 intersection, regional support, photographer QC, ecological hypervolume or functional trait inference is established.",
    }




def validate_existing_overlap_receipt(ledger: dict, receipt: dict) -> dict:
    """Check an independently archived, NON-IDENTIFYING trait intersection receipt.

    This uses an existing frozen source aggregate and cannot substitute for
    the observer/cell QC or a fresh join of individual native source rows.
    """
    if receipt.get("schema_version") != "native_gate_a_existing_aggregate_overlap_20261020_v1":
        raise ValueError("Unrecognized historical native intersection receipt schema")
    source = receipt["original_source"]
    if source["native_observations"] != NATIVE_N or source["inherited_measurement_rows"] != MEASURED_N:
        raise ValueError("Existing intersection source denominator changed")
    if source["accepted_taxa"] != ledger["denominators"]["native_wcvp_accepted_keys"]:
        raise ValueError("Existing intersection taxon universe changed")
    p = receipt["provenance"]
    if p["workflow_run_id"] != 39428065666 or p["artifact_id"] != 12266812186:
        raise ValueError("Existing intersection artifact origin changed")
    if p["artifact_member_sha256"] != "39d5613b611cce30be989682b3e47cdb3153aa880a0c2744cd40bb0d72b9edd4":
        raise ValueError("Existing intersection member hash changed")
    counts = {e["endpoint_id"]: e["n_eligible_native"] for e in ledger["endpoints"]}
    by = receipt["module_census"]
    expected = {
        "orientation": ("orientation_image_vertical_angle",),
        "colour": ("corolla_lab_lightness", "corolla_lab_chroma", "corolla_hue_sin", "corolla_hue_cos"),
        "outline": ("capitulum_outline_aspect_ratio", "capitulum_outline_circularity", "capitulum_outline_solidity", "capitulum_width_profile_cv"),
        "bract_si": ("bract_projection_p95", "bract_spread_fraction"),
    }
    for name, fields in expected.items():
        n = int(by[name]["n"])
        if not (0 < n <= MEASURED_N) or not all(n <= counts[k] for k in fields):
            raise ValueError("Invalid module marginal or missing registered endpoint: " + name)
        if name in ("orientation", "outline", "bract_si") and not all(counts[k] == n for k in fields):
            raise ValueError("Module marginal disagrees with identical original endpoint eligibility: " + name)
        if name == "colour" and n != min(counts[k] for k in fields):
            raise ValueError("Colour4 module does not have the expected common-hue gate")
        if not 0 < by[name]["taxa"] <= source["accepted_taxa"]:
            raise ValueError("Invalid module taxon count")
    pair = receipt["module_overlap"]
    for pair_name, value in pair.items():
        parts = pair_name.split("/")
        if len(parts) != 2 or parts[0] == parts[1] or any(x not in by for x in parts):
            raise ValueError("Invalid module pair key")
        n = value["n"]
        if not max(0, by[parts[0]]["n"] + by[parts[1]]["n"] - MEASURED_N) <= n <= min(by[x]["n"] for x in parts):
            raise ValueError("Historical module pair violates original source intersection bounds")
        if not 0 < value["taxa"] <= min(by[x]["taxa"] for x in parts):
            raise ValueError("Historical module pair taxon count inconsistent")
    core = receipt["complete_core9"]
    nc = core["n"]
    if nc != 157779 or core["accepted_taxa"] != 242:
        raise ValueError("Existing complete-nine intersection changed")
    if not all(nc <= pair[name]["n"] for name in ("orientation/colour", "orientation/outline", "colour/outline")):
        raise ValueError("Complete-nine cohort exceeds a necessary pairwise joint cohort")
    bounds = guaranteed_bounds(ledger)["intersections"]
    if not bounds["complete_core9"]["lower"] <= nc <= bounds["complete_core9"]["upper"]:
        raise ValueError("Complete-nine overlap violates original source count bounds")
    extra_colour = counts["corolla_lab_chroma"] - by["colour"]["n"]
    if extra_colour != 49:
        raise ValueError("Colour4-versus-chroma difference unexpected")
    # The original four outline endpoints share identical eligibility here;
    # the COMPLETE colour module is nested in the chroma endpoint set.
    interval = {"lower": nc, "upper": min(pair["orientation/outline"]["n"], nc + extra_colour)}
    if interval != {
        "lower": receipt["core3_inferred_bound"]["lower"],
        "upper": receipt["core3_inferred_bound"]["upper"],
    }:
        raise ValueError("Derived core3 observation intersection range changed")
    return {
        "status": "VERIFIED_EXISTING_AGGREGATE_NOT_REPLAY",
        "complete_core9_n": nc,
        "complete_core9_taxa": core["accepted_taxa"],
        "core3_observation_range": interval,
        "photographer_join_required": not receipt["sampling_availability"]["user_id"],
        "no_new_trait_values_estimated": True,
    }



def validate_existing_sampling_context(receipt: dict) -> dict:
    """Validate scope of historical SOURCE-WIDE effort and single-trait cells.

    The source summary must never be mistaken for photographer/cell support
    of the exact nine-trait matched observation cohort.
    """
    c = receipt["source_sampling_context"]
    if c["scope"] != "ORIGINAL_313350_NATIVE_RECORDS_NOT_MATCHED_CORE9_OBSERVER_DATA":
        raise ValueError("Incorrect native sample scope")
    e = c["evidence"]
    if (e["spatial_run_id"], e["spatial_artifact_id"]) != (39428065666, 12266812186):
        raise ValueError("Historical spatial run identity changed")
    if e["spatial_member_sha256"] != "d62e01810179334303d0d07930e070feb86af0c4cd6c08ea59fb41ea091e0430":
        raise ValueError("Historical spatial output digest changed")
    if e["integration_member_sha256"] != "1545037041b62b2494934d022c7b21e7a7109f54f3582579c6557dfadf998752":
        raise ValueError("Historical model-cohort output digest changed")
    a = c["original_observation_effort"]
    if (a["native_observations"], a["taxa"], a["observers"]) != (313350, 306, 92839):
        raise ValueError("Native source effort denominator changed")
    fields = ("top2_taxa_record_fraction", "top5_taxa_record_fraction",
              "top1pct_observer_record_fraction", "top10pct_observer_record_fraction",
              "same_taxon_cell_observer_repeated_record_fraction")
    if any(not 0 < a[k] < 1 for k in fields):
        raise ValueError("Invalid source posting fraction")
    if a["top2_taxa_record_fraction"] > a["top5_taxa_record_fraction"]:
        raise ValueError("Taxon posting fraction nonmonotone")
    if a["top1pct_observer_record_fraction"] > a["top10pct_observer_record_fraction"]:
        raise ValueError("Observer posting fraction nonmonotone")
    cells = c["individual_endpoint_coarse_spatial_support"]
    if not cells["not_joint_core9"] or not cells["not_independent_plants"]:
        raise ValueError("Single-endpoint cells falsely promoted to joint plants")
    known = {
        "orientation_image_vertical_angle": (170066, 244, 1421),
        "corolla_lab_chroma": (186439, 244, 1437),
        "corolla_hue_sin": (186390, 244, 1437),
        "capitulum_outline_aspect_ratio": (176048, 243, 1421),
        "bract_projection_p95": (1303, 111, 509),
    }
    for key, (n, taxa, ncells) in known.items():
        v = cells["endpoints"][key]
        if (v["main_n"], v["taxa"], v["one_per_taxon_cell_n"]) != (n, taxa, ncells):
            raise ValueError("Historical coarse-cell reference changed: " + key)
    m = c["previously_fitted_complete_core9"]
    if (m["measured_native_n"],m["measured_native_taxa"],m["fitted_n"],m["fitted_taxa"]) != (157779,242,157654,174):
        raise ValueError("Model-admitted taxa are not original measurement support")
    if (receipt["complete_core9"]["n"], receipt["complete_core9"]["accepted_taxa"]) != (157779,242):
        raise ValueError("Core9 measured cohort changed")
    old = receipt["legacy_five_prior_native_recovery"]
    if old["historical_pr"] != 92 or old["historical_pr_head_sha"] != "6061befaf56ee545bb3b2b4a80ed220d2def1190":
        raise ValueError("Legacy five restored trait source branch changed")
    if old["historical_receipt_blob_sha"] != "92498fb6a20418918e1d5e9b79e168c13c17197c":
        raise ValueError("Legacy five restored trait receipt SHA changed")
    if (old["historical_workflow_run_id"],old["historical_artifact_id"]) != (34334322764,10097484190):
        raise ValueError("Legacy five recovery execution identity changed")
    if (old["historical_native_cohort_observations"],old["restored_native_observations_per_endpoint"]) != (27066,23881):
        raise ValueError("Legacy five recovered observation counts changed")
    if len(set(old["restored_endpoints"])) != 5 or set(old["restored_endpoints"]) != set(RECOVERED_FIVE):
        raise ValueError("Legacy five restored endpoint names changed")
    if old["new_313350_native_roster_join_status"] != "NOT_VERIFIED" or old["new_313350_native_restored_observations_per_endpoint"] is not None:
        raise ValueError("Do not treat old 27k native sample as new 313k native measurements")
    coarse_ceiling = coarse_cell_retention_ceiling(receipt)
    a5 = receipt["legacy_five_archive_metadata"]
    if (a5["documented_artifact_id"],a5["parts_recorded"],
        a5["parts_with_existing_drive_metadata_found"]) != (8269246732,20,20):
        raise ValueError("Legacy head archive metadata inventory changed")
    if a5["archive_raw_bytes_reassembled_and_verified"] or a5["original_1255791_head_row_fields_rechecked"]:
        raise ValueError("Cannot certify original head-level content from Drive metadata")
    if a5["native_joined_five_endpoint_counts"] != "NOT_VERIFIED":
        raise ValueError("Historical five traits not source-certified")
    for value in c["pending_exact_core9_sampling"].values():
        if value != "NOT_VERIFIED":
            raise ValueError("Source-wide totals cannot close matched-core9 sampling gate")
    return {
        "status":"HISTORICAL_SOURCE_EFFORT_AND_SPATIAL_SUPPORT_VERIFIED_ONLY",
        "source_observations":313350,
        "source_observers":92839,
        "core9_measured_observations":157779,
        "core9_fitted_previous_observations":157654,
        "core9_fitted_previous_taxa":174,
        "core9_original_observer_strata":"NOT_VERIFIED",
        "legacy_head_parts_with_metadata":20,
        "legacy_five_native_values":"NOT_VERIFIED",
        "legacy_five_restored_observations_per_endpoint_on_historical_27066_frame":23881,
        "historical_coarse_cell_retention_ceiling":coarse_ceiling
    }


def coarse_cell_retention_ceiling(receipt: dict) -> dict:
    """Sharp *upper* counts from prior single-trait coarse cells, not actual support.

    If all 242 complete-core9 taxa are required to have >=k distinct cells,
    they collectively consume >=242+(k-1)*qualifying_taxa occupied taxon-cells.
    The known 244-taxon orientation margin has 1421 occupied cells, including
    at least one cell for each of its two taxa not in the 242-taxon core9.
    This bound applies to that exact archived coarse grid ONLY and cannot
    constrain a finer 0.5-degree grid or prove individual-plant independence.
    """
    evidence = receipt["source_sampling_context"]["individual_endpoint_coarse_spatial_support"]
    if not evidence["not_joint_core9"] or not evidence["not_independent_plants"]:
        raise ValueError("Invalid source spatial provenance")
    source = evidence["endpoints"]["orientation_image_vertical_angle"]
    n_taxa = int(receipt["complete_core9"]["accepted_taxa"])
    source_taxa, source_cells = source["taxa"],source["one_per_taxon_cell_n"]
    if source_taxa < n_taxa or source_cells < source_taxa:
        raise ValueError("Impossible prior spatial source/target taxon counts")
    # All core9 observations are also eligible orientation observations.
    maximum_core9_cells = source_cells - (source_taxa - n_taxa)
    available_excess = maximum_core9_cells - n_taxa
    if available_excess < 0:
        raise ValueError("No one-cell-per-core-taxon allocation is possible")
    maximum = {
        str(k): min(n_taxa, available_excess // (k-1))
        for k in (5,10,20,50)
    }
    return {
        "status":"MATHEMATICAL_UPPER_BOUNDS_NOT_ACTUAL_CELL_COVERAGE",
        "reference":"historical one-per-accepted-taxon-per-coarse-sine-latitude-cell",
        "source_orientation_taxon_cell_count":source_cells,
        "source_orientation_accepted_taxa":source_taxa,
        "matched_nine_accepted_taxa":n_taxa,
        "upper_bound_matched_nine_taxon_cells":maximum_core9_cells,
        "maximum_matched_nine_taxa_with_at_least_k_cells":maximum,
        "caution":"These upper bounds apply only to the historical coarse cell grid. Different/finer grid and observer independence are NOT inferred."
    }


def _strict_endpoint_masks(frame: pd.DataFrame, ledger: dict | None) -> tuple[dict, dict]:
    endpoint_ids = [r["endpoint_id"] for r in ledger["endpoints"]] if ledger else list(
        dict.fromkeys((*CORE9, *RECOVERED_FIVE))
    )
    missing = sorted((set(endpoint_ids) | {x + "__eligible" for x in endpoint_ids}) - set(frame.columns))
    if missing:
        raise ValueError("Missing endpoint columns: " + ", ".join(missing))
    coverage, masks = {}, {}
    expected = {r["endpoint_id"]: r for r in ledger["endpoints"]} if ledger else {}
    for name in endpoint_ids:
        raw = frame[name]
        values = pd.to_numeric(raw, errors="coerce")
        finite = np.isfinite(values.to_numpy(float))
        flags = frame[name + "__eligible"].astype("string").str.strip().str.lower()
        if not flags.isin(("true", "false")).all():
            raise ValueError("Invalid/missing endpoint-eligibility flag: " + name)
        true = flags.eq("true").to_numpy(dtype=bool)
        if np.any(true & ~finite):
            raise ValueError("Eligible nonfinite measurement: " + name)
        used = true & finite
        record = {
            "n_finite": int(finite.sum()),
            "n_eligible": int(used.sum()),
            "n_finite_but_ineligible": int((finite & ~used).sum()),
            "n_unavailable_or_nonfinite": int((~finite).sum()),
            "qc_failure_reason": "not_disaggregated_in_prepared_table",
        }
        if name in expected:
            a = expected[name]
            if (record["n_finite"], record["n_eligible"]) != (
                a["n_finite_native"], a["n_eligible_native"]
            ):
                raise ValueError("Endpoint coverage differs from archived receipt: " + name)
        coverage[name] = record
        masks[name] = used
    return coverage, masks


def _spatial_index(frame: pd.DataFrame) -> pd.Series:
    lat = pd.to_numeric(frame["latitude"], errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame["longitude"], errors="coerce").to_numpy(float)
    valid = np.isfinite(lat) & np.isfinite(lon) & (lat >= -90) & (lat <= 90)
    valid &= (lon >= -180) & (lon <= 180)
    result = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    if valid.any():
        bands = np.floor((np.sin(np.deg2rad(lat[valid])) + 1) * 45).clip(0, 89).astype(int)
        longitude = np.floor((lon[valid] + 180) / 2).clip(0, 179).astype(int)
        result.iloc[np.flatnonzero(valid)] = bands * 180 + longitude
    return result



def sampling_selection_diagnostics(
    frame: pd.DataFrame,
    matched: np.ndarray,
    spatial: pd.Series,
    dates: pd.Series | None,
    observers: pd.Series | None,
) -> dict:
    """Exact cohort sampling-effort audit; emits no individual/observer/taxon IDs.

    Source and matched comparisons use the same *original native roster*. The
    cells are the coarse approximate equal-area Gate A cells, NOT individual
    plants, independent regions, or the Gate B 0.5-degree angular grid.
    """
    if len(matched) != len(frame) or len(spatial) != len(frame):
        raise ValueError("Source/matched cohort length mismatch")
    mask = np.asarray(matched, dtype=bool)
    source = frame[["accepted_key"]].copy()
    source["cell"] = spatial.to_numpy()
    eligible = source.loc[mask].dropna(subset=["cell"]).copy()
    group_n = eligible.groupby(["accepted_key", "cell"], sort=False).size()
    source_cell_n = source.dropna(subset=["cell"]).groupby(
        ["accepted_key", "cell"], sort=False
    ).size()
    result = {
        "scope": "EXACT_MATCHED_OBSERVATIONS_NOT_INDEPENDENT_PLANTS",
        "source_observations": int(len(frame)),
        "matched_observations": int(mask.sum()),
        "matched_taxon_cells": int(len(group_n)),
        "source_taxon_cells": int(len(source_cell_n)),
        "matched_taxon_cells_with_repeated_observations": int(group_n.gt(1).sum()),
        "matched_observations_in_repeated_taxon_cells": int(group_n.loc[group_n.gt(1)].sum()),
        "coarse_cell_method": "2-degree longitude by equal-sin-latitude bands, NOT 0.5-degree Gate B cells",
        "observer_status": "NOT_EVALUABLE_NO_ORIGINAL_OBSERVER_JOIN",
        "season_status": "NOT_EVALUABLE_NO_OBSERVATION_DATE",
        "photographic_quality_status": "NOT_EVALUABLE_NO_SOURCE_IMAGE_QUALITY_OR_CAMERA_POSE_FLAGS",
        "endpoint_failure_reason_status": "NOT_DISAGGREGATED_IN_PREPARED_TABLE",
        "source_vs_matched_monthly": None,
        "matched_observer_effort": None,
    }
    if dates is not None:
        months = dates.dt.month
        source_counts = months.value_counts(dropna=True)
        matched_counts = months.loc[mask].value_counts(dropna=True)
        result["season_status"] = "SOURCE_VS_MATCHED_RECORD_MONTH_ONLY_NOT_FLOWERING_STAGE"
        result["source_n_observations_with_date"] = int(months.notna().sum())
        result["matched_n_observations_with_date"] = int(months.loc[mask].notna().sum())
        result["source_vs_matched_monthly"] = [
            {"month": k,
             "source_n": int(source_counts.get(k, 0)),
             "matched_n": int(matched_counts.get(k, 0))}
            for k in range(1, 13)
        ]
    if observers is not None:
        if len(observers) != len(frame) or observers.isna().any():
            raise ValueError("Incomplete original observer roster")
        selected_observers = observers.loc[mask]
        observer_counts = selected_observers.value_counts()
        top = max(1, int(np.ceil(0.01 * len(observer_counts)))) if len(observer_counts) else 0
        with_cell = eligible.copy()
        with_cell["observer"] = observers.loc[with_cell.index].to_numpy()
        triple = with_cell.groupby(["accepted_key", "cell", "observer"], sort=False).size()
        result["observer_status"] = "EXACT_MATCHED_ORIGINAL_OBSERVER_SUPPORT"
        result["matched_observer_effort"] = {
            "n_original_observers": int(selected_observers.nunique()),
            "n_taxon_observer_pairs": int(
                source.loc[mask].assign(observer=selected_observers.to_numpy())
                .drop_duplicates(["accepted_key", "observer"]).shape[0]
            ),
            "n_taxon_cell_observer_strata": int(len(triple)),
            "n_taxon_cell_observer_strata_with_repeats": int(triple.gt(1).sum()),
            "matched_observations_in_repeated_taxon_cell_observer_strata": int(
                triple.loc[triple.gt(1)].sum()
            ),
            "matched_top_1pct_observer_fraction": (
                float(observer_counts.iloc[:top].sum() / int(mask.sum()))
                if mask.any() else None
            ),
        }
    return result


def audit_gate_a(
    frame: pd.DataFrame,
    *,
    ledger: dict | None = None,
    observers: pd.DataFrame | None = None,
    expected_n: int | None = None,
) -> dict:
    """Exact observation-level intersection; only deidentified aggregate outputs."""
    essentials = {"obs_id", "accepted_key", "latitude", "longitude"}
    if not essentials.issubset(frame):
        raise ValueError("Native observation table missing required metadata")
    if expected_n is not None and len(frame) != expected_n:
        raise ValueError("Wrong native source count")
    if frame.obs_id.isna().any() or not frame.obs_id.is_unique:
        raise ValueError("Missing/duplicated observation ID")
    if frame.accepted_key.isna().any():
        raise ValueError("Native accepted taxon must be resolved")

    coverage, masks = _strict_endpoint_masks(frame, ledger)
    if ledger is not None:
        if "archived_measurement_row" not in frame:
            raise ValueError("Historical measurement-row flag missing")
        measured = frame.archived_measurement_row.astype("string").str.lower()
        if not measured.isin(("true", "false")).all():
            raise ValueError("Historical measurement flag invalid")
        count = int(measured.eq("true").sum())
        if count != ledger["denominators"]["native_with_inherited_measurement_row"]:
            raise ValueError("Inherited measurement-row count does not match source receipt")

    if observers is not None:
        if not {"obs_id", "user_id"}.issubset(observers):
            raise ValueError("Observer join must contain obs_id and user_id")
        if observers.obs_id.isna().any() or observers.obs_id.duplicated().any():
            raise ValueError("Duplicate or missing original observer ID key")
        if set(observers.obs_id.astype(str)) != set(frame.obs_id.astype(str)):
            raise ValueError("Observer roster must match exactly the native observation roster")
        by_id = observers.set_index("obs_id")["user_id"]
        observer_ids = frame.obs_id.map(by_id)
        if observer_ids.isna().any():
            raise ValueError("Observer IDs missing after source matching")
    else:
        observer_ids = None

    spatial = _spatial_index(frame)
    dates = (
        pd.to_datetime(frame["observed_on"], errors="coerce")
        if "observed_on" in frame else None
    )
    source_counts = frame["accepted_key"].value_counts()
    source_top5 = sorted(source_counts.index, key=lambda x: (-int(source_counts[x]), str(x)))[:5]
    result = {
        "status": "actual_matched_phenotype_eligibility_only" if ledger else "synthetic_qa_only",
        "native_observations": int(len(frame)),
        "accepted_taxa": int(frame.accepted_key.nunique()),
        "endpoint_coverage": coverage,
        "modules": {},
        "intersections": {},
        "observer_audit": "admitted_private_original_observers" if observer_ids is not None else "NOT_EVALUABLE_NO_ORIGINAL_OBSERVER_JOIN",
        "spatial_method": "approximately equal-area 2-degree-longitude x uniform-sin(latitude) cells; not independent plants",
        "disclosure": "Aggregates only. No observation, taxon, observer, locality or photo identifier is exported.",
    }
    for name, members in MODULES.items():
        mask = np.logical_and.reduce([masks[x] for x in members])
        sample = frame.loc[mask, ["accepted_key"]].copy()
        sample["cell"] = spatial.loc[mask]
        group = sample.dropna(subset=["cell"]).drop_duplicates(["accepted_key", "cell"])
        ncells = group.groupby("accepted_key", sort=True).size()
        record = {
            "n_observations_exact": int(mask.sum()),
            "n_accepted_taxa": int(sample.accepted_key.nunique()),
            "n_with_usable_coordinates": int(sample.cell.notna().sum()),
            "n_unique_taxon_cells": int(len(group)),
            "n_taxa_with_min_5_cells": int((ncells >= 5).sum()),
            "n_taxa_with_min_10_cells": int((ncells >= 10).sum()),
            "n_taxa_with_min_20_cells": int((ncells >= 20).sum()),
            "n_taxa_with_min_50_cells": int((ncells >= 50).sum()),
            "top2_source_taxa_fraction_of_module_observations": (
                float(frame.loc[mask, "accepted_key"].isin(source_top5[:2]).mean())
                if mask.any() else None
            ),
            "top5_source_taxa_fraction_of_module_observations": (
                float(frame.loc[mask, "accepted_key"].isin(source_top5).mean())
                if mask.any() else None
            ),
            "n_observations_with_observed_on_date": (
                int(dates.loc[mask].notna().sum()) if dates is not None else None
            ),
            "n_distinct_original_observers": (
                int(observer_ids.loc[mask].nunique()) if observer_ids is not None else None
            ),
        }
        result["modules"][name] = record
    result["matched_core9_sampling_diagnostics"] = sampling_selection_diagnostics(
        frame, np.logical_and.reduce([masks[x] for x in CORE9]),
        spatial, dates, observer_ids
    )
    for names in itertools.combinations(CORE3, 2):
        mask = masks[names[0]] & masks[names[1]]
        result["intersections"]["__".join(names)] = int(mask.sum())
    result["intersections"]["all_core3"] = result["modules"]["primary_core3"]["n_observations_exact"]
    result["intersections"]["all_core9"] = result["modules"]["complete_primary_core9"]["n_observations_exact"]
    if ledger is not None:
        bounds = guaranteed_bounds(ledger)["intersections"]
        check = {
            "orientation_x_chroma": result["intersections"]["__".join(CORE3[:2])],
            "orientation_x_outline_aspect": result["intersections"]["__".join((CORE3[0], CORE3[2]))],
            "chroma_x_outline_aspect": result["intersections"]["__".join(CORE3[1:])],
            "matched_core3": result["intersections"]["all_core3"],
            "complete_core9": result["intersections"]["all_core9"],
        }
        for name, val in check.items():
            lim = bounds[name]
            if not lim["lower"] <= val <= lim["upper"]:
                raise ValueError("Actual observation join violates source-count bounds: " + name)
        result["receipt_consistency"] = "PASS_27_ENDPOINT_COUNTS_AND_CORE_INTERSECTION_BOUNDS"
    return result


def audit_recovered_five(
    native_frame: pd.DataFrame, observation_csv: Path, report_json: Path
) -> dict:
    """Optional reuse of PR92 recovered HEAD data; never assume it was executed."""
    prior = json.loads(report_json.read_text(encoding="utf-8"))
    if prior.get("status") != "FIVE_PREVIOUSLY_MEASURED_FIELDS_RECOVERED_TO_VERSIONED_OBSERVATION_VIEW":
        raise ValueError("Unverified historical head aggregation")
    if prior.get("source_head_sha256") != HEAD_ARCHIVE_SHA256:
        raise ValueError("Historical head-level source identity mismatch")
    if prior.get("counts", {}).get("source_heads") != HEAD_ARCHIVE_N:
        raise ValueError("Historical head sample count mismatch")
    if prior.get("counts", {}).get("head_records_deleted") != 0:
        raise ValueError("Historical head loss")
    if prior.get("new_image_operations") is not False:
        raise ValueError("Unexpected remeasurement")
    if sha256(observation_csv) != prior.get("observation_values_sha256"):
        raise ValueError("Recovered observation CSV hash mismatch")
    needed = ["obs_id", "n_visible_photos", "n_composition_photos", *RECOVERED_FIVE.values()]
    recovered = pd.read_csv(observation_csv, usecols=needed, dtype={"obs_id": str})
    if recovered.obs_id.isna().any() or not recovered.obs_id.is_unique:
        raise ValueError("Duplicate/missing historical source observation")
    linked = recovered.loc[recovered.obs_id.isin(native_frame.obs_id)].copy()
    v = pd.to_numeric(linked["n_visible_photos"], errors="coerce").fillna(0).gt(0)
    comp = pd.to_numeric(linked["n_composition_photos"], errors="coerce").fillna(0).gt(0)
    vals = linked[list(RECOVERED_FIVE.values())].apply(pd.to_numeric, errors="coerce")
    visible = v & vals.iloc[:, 0].between(0, 1)
    fractions = vals.iloc[:, 1:]
    closure = (fractions.notna().all(axis=1) & fractions.ge(0).all(axis=1)
               & fractions.le(1).all(axis=1)
               & fractions.sum(axis=1).sub(1).abs().le(1e-6))
    if (comp & ~closure).any():
        raise ValueError("Recomputed observation-level colour fractions are not closed")
    return {
        "status": "historical_recovered_head_data_linked_to_native_source_candidate",
        "historical_observation_rows": int(len(recovered)),
        "historical_observations_matching_native_source": int(len(linked)),
        "native_observations_with_valid_visible_fraction": int(visible.sum()),
        "native_observations_with_closed_four_colour_fractions": int((comp & closure).sum()),
        "n_native_accepted_taxa_with_recovered_colour_composition": int(
            native_frame.loc[native_frame.obs_id.isin(linked.loc[comp & closure, "obs_id"]),
                             "accepted_key"].nunique()
        ),
        "taxon_label_cross_check": "unverified: historical recovered file lacks original source taxon",
        "interpretation": "Optional historical observation aggregate. Not a replacement for frozen zero-eligible V2 endpoint atlas.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=Path(
        __file__).with_name("native_first_27_endpoint_coverage_20261020.json"))
    parser.add_argument("--native-input", type=Path, help="Existing private prepared native observation CSV.gz")
    parser.add_argument("--observer-input", type=Path, help="Optional original private obs_id,user_id CSV")
    parser.add_argument("--recovered-five-csv", type=Path)
    parser.add_argument("--recovered-five-report", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="Aggregate-only JSON")
    args = parser.parse_args()
    ledger = load_ledger(args.ledger)
    existing_path = Path(__file__).with_name("native_gate_a_existing_overlap_20261020.json")
    existing = json.loads(existing_path.read_text(encoding="utf-8"))
    verified = validate_existing_overlap_receipt(ledger, existing)
    existing_sampling = validate_existing_sampling_context(existing)
    output = {"bounds_only": guaranteed_bounds(ledger),
              "existing_aggregate_audit": verified,
              "historical_sampling_context": existing_sampling}
    if bool(args.recovered_five_csv) != bool(args.recovered_five_report):
        parser.error("Historical five requires both file and provenance report")
    if (args.observer_input or args.recovered_five_csv) and not args.native_input:
        parser.error("Private joins require --native-input")
    if args.native_input:
        source_hash = sha256(args.native_input)
        if source_hash == PREPARED_SHA256:
            raise ValueError(
                "Uncorrected environment input: original GSP sentinels and "
                "CHELSA affine decoding are superseded by PR143"
            )
        if source_hash != CORRECTED_PREPARED_SHA256:
            raise ValueError("Unknown native input bytes; expected corrected PR143 table")
        native = pd.read_csv(args.native_input, dtype={"obs_id": str}, low_memory=False)
        observers = None
        if args.observer_input:
            observers = pd.read_csv(
                args.observer_input, usecols=["obs_id", "user_id"],
                dtype={"obs_id": str, "user_id": str}
            )
        output["actual_gate_a"] = audit_gate_a(
            native, ledger=ledger, observers=observers, expected_n=NATIVE_N
        )
        verified_count = output["actual_gate_a"]["modules"]["complete_primary_core9"]["n_observations_exact"]
        if verified_count != verified["complete_core9_n"]:
            raise ValueError("Current source exact core9 does not match archived aggregate receipt")
        exact_core3 = output["actual_gate_a"]["modules"]["primary_core3"]["n_observations_exact"]
        interval = verified["core3_observation_range"]
        if not interval["lower"] <= exact_core3 <= interval["upper"]:
            raise ValueError("Current source core3 falls outside archived aggregate-derived interval")
        output["actual_gate_a"]["source_prepared_table_sha256"] = source_hash
        output["actual_gate_a"]["environment_source_status"] = "PR143_CORRECTED_SHA256_LOCKED"
        if args.recovered_five_csv:
            output["historical_five_linkage"] = audit_recovered_five(
                native, args.recovered_five_csv, args.recovered_five_report
            )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({
        "actual_source_analysed": "actual_gate_a" in output,
        "core3_lower": output["existing_aggregate_audit"]["core3_observation_range"]["lower"],
        "core3_upper": output["existing_aggregate_audit"]["core3_observation_range"]["upper"],
        "existing_core9": output["existing_aggregate_audit"]["complete_core9_n"],
        "exact_core3": output.get("actual_gate_a", {}).get("intersections", {}).get("all_core3"),
    }))


if __name__ == "__main__":
    main()

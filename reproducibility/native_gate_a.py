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
    output = {"bounds_only": guaranteed_bounds(ledger),
              "existing_aggregate_audit": verified}
    if bool(args.recovered_five_csv) != bool(args.recovered_five_report):
        parser.error("Historical five requires both file and provenance report")
    if (args.observer_input or args.recovered_five_csv) and not args.native_input:
        parser.error("Private joins require --native-input")
    if args.native_input:
        source_hash = sha256(args.native_input)
        if source_hash != PREPARED_SHA256:
            raise ValueError("Private prepared native-observation checksum mismatch")
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
        output["actual_gate_a"]["source_prepared_table_sha256"] = PREPARED_SHA256
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

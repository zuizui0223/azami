"""Optional 44-axis RF comparison on the original private S1 photographic traits.

Only paired, held-out aggregate diagnostics may leave the protected runner. No
new result or provenance is asserted by the public source-only tests. Existing
V4 manuscript/accepted S1-S2/H3c numerical families are never touched.

Real source extraction is external to this module: a 32-column extended table
and a per-layer acquisition manifest must be supplied by the authorised runner.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_multitrait_predictability import (
    CORE4, FULL12, SEED, aggregate_endpoint, score_models, sha256,
    verified_source,
)

# Columns are declared, not claimed to have been downloaded or measured.
BIO9 = tuple("chelsa_bio" + x for x in
             ("02", "03", "07", "08", "09", "10", "11", "16", "19"))
WATER5 = (
    "tc_vpd_mean_1981_2010", "tc_pet_annual_1981_2010",
    "tc_aet_annual_1981_2010", "tc_climatic_deficit_annual_1981_2010",
    "tc_soil_moisture_mean_1981_2010",
)
WIND2 = ("tc_wind_mean_1981_2010", "tc_wind_seasonality_1981_2010")
TERRAIN5 = (
    "topo_slope_deg", "topo_roughness", "topo_tpi",
    "topo_northness", "topo_eastness",
)
SOIL6 = (
    "soil_phh2o_0_5cm", "soil_soc_0_5cm", "soil_clay_0_5cm",
    "soil_sand_0_5cm", "soil_bdod_0_5cm", "soil_cec_0_5cm",
)
VEG5 = (
    "veg_ndvi_mean_2001_2010", "veg_ndvi_seasonal_sd_2001_2010",
    "landcover_tree_fraction_2015", "landcover_crop_fraction_2015",
    "landcover_grass_fraction_2015",
)
NEW32 = BIO9 + WATER5 + WIND2 + TERRAIN5 + SOIL6 + VEG5
STAGES = (
    ("core4", CORE4),
    ("archived12", FULL12),
    ("complete_bioclim21", FULL12 + BIO9),
    ("water_balance26", FULL12 + BIO9 + WATER5),
    ("wind28", FULL12 + BIO9 + WATER5 + WIND2),
    ("terrain33", FULL12 + BIO9 + WATER5 + WIND2 + TERRAIN5),
    ("soil39", FULL12 + BIO9 + WATER5 + WIND2 + TERRAIN5 + SOIL6),
    ("vegetation44", FULL12 + NEW32),
)
TIERS = dict(STAGES)
PAIRS = tuple((STAGES[i - 1][0], STAGES[i][0])
              for i in range(1, len(STAGES)))
GROUPS44 = {
    "temperature": tuple(n for n in FULL12 + BIO9 if
                         n in ("chelsa_bio01", "chelsa_bio02", "chelsa_bio03",
                               "chelsa_bio04", "chelsa_bio05", "chelsa_bio06",
                               "chelsa_bio07", "chelsa_bio08", "chelsa_bio09",
                               "chelsa_bio10", "chelsa_bio11")),
    "precipitation": tuple(n for n in FULL12 + BIO9 if
                           n in ("chelsa_bio12", "chelsa_bio13", "chelsa_bio14",
                                 "chelsa_bio15", "chelsa_bio16", "chelsa_bio17",
                                 "chelsa_bio18", "chelsa_bio19")),
    "radiation": ("chelsa_rsds_mean",),
    "elevation": ("topo_elevation",),
    "water_balance": WATER5,
    "wind": WIND2,
    "terrain": TERRAIN5,
    "soil": SOIL6,
    "vegetation_landcover": VEG5,
}
DATASET_FAMILIES = {
    **{n: "CHELSA-bioclim-v2.1" for n in BIO9},
    **{n: "TerraClimate" for n in WATER5 + WIND2},
    **{n: "WorldClim-elevation-derived" for n in TERRAIN5},
    **{n: "SoilGrids-v2.0" for n in SOIL6},
    **{n: "MOD13A3-v061" for n in VEG5[:2]},
    **{n: "Copernicus-LC100-C3-2015" for n in VEG5[2:]},
}
EXPECTED_PERIODS = {
    **{n: "1981-2010" for n in BIO9 + WATER5 + WIND2},
    **{n: "static-topography" for n in TERRAIN5},
    **{n: "static-soil-model" for n in SOIL6},
    **{n: "2001-2010" for n in VEG5[:2]},
    **{n: "2015" for n in VEG5[2:]},
}
HEX = re.compile(r"[a-f0-9]{64}\Z")

assert len(FULL12) == 12
assert len(NEW32) == 32 and len(TIERS["vegetation44"]) == 44
assert len(set(TIERS["vegetation44"])) == 44
assert set().union(*map(set, GROUPS44.values())) == set(TIERS["vegetation44"])
assert sum(map(len, GROUPS44.values())) == 44
assert all(set(left) < set(right) for (_, left), (_, right) in
           zip(STAGES, STAGES[1:]))


def added_columns(upper_tier: str) -> tuple[str, ...]:
    if upper_tier not in TIERS or upper_tier in ("core4", "archived12"):
        raise ValueError("EXPANDED_TIER_FROM_21_TO_44_REQUIRED")
    return tuple(n for n in TIERS[upper_tier] if n not in FULL12)


def verify_extended_source(frame: pd.DataFrame, table_path: Path,
                           table_sha256: str, manifest_path: Path,
                           *, required_columns: tuple[str, ...] = NEW32):
    """Join 32 *declared* environmental extractions to the exact native S1 IDs.

    The manifest checks provenance declarations, not unprovided raster bytes.
    The extraction receipt must be verified separately in the secure environment.
    """
    if (not required_columns or len(required_columns) != len(set(required_columns)) or
            not set(required_columns).issubset(NEW32)):
        raise ValueError("INVALID_EXTRA_LAYER_SCOPE")
    if not HEX.fullmatch(table_sha256) or sha256(table_path) != table_sha256:
        raise ValueError("EXPANDED_32_AXIS_TABLE_SHA_MISMATCH")
    manifest_bytes = Path(manifest_path).read_bytes()
    manifest = json.loads(manifest_bytes)
    if (manifest.get("schema") != "azami_rf_extra32_source_manifest_v1" or
            manifest.get("table_sha256") != table_sha256 or
            set(manifest.get("layers", {})) != set(required_columns)):
        raise ValueError("EXTENDED_LAYER_MANIFEST_INCOMPLETE")
    for name in required_columns:
        spec = manifest["layers"][name]
        if (spec.get("dataset") != DATASET_FAMILIES[name] or
                spec.get("temporal_support") != EXPECTED_PERIODS[name] or
                not HEX.fullmatch(str(spec.get("source_bundle_sha256", ""))) or
                (name in SOIL6 and (spec.get("depth_interval") != "0-5cm" or
                                      spec.get("statistic") != "Q0.50")) or
                not all(isinstance(spec.get(k), str) and len(spec[k].strip()) > 2
                        for k in ("unit", "derivation", "sampling_method", "source_url"))):
            raise ValueError("EXTENDED_LAYER_UNVERIFIED_DECLARATION: " + name)
    if not frame.obs_id.is_unique or frame.obs_id.isna().any():
        raise ValueError("ORIGINAL_NATIVE_DUPLICATE_OR_MISSING_IDS")
    extra = pd.read_csv(table_path, dtype={"obs_id": str}, low_memory=False)
    if (set(extra.columns) != {"obs_id", *required_columns} or
            len(extra) != len(frame) or extra.obs_id.isna().any() or
            not extra.obs_id.is_unique):
        raise ValueError("EXPANDED_TABLE_32_AXIS_OR_IDENTITY_INVALID")
    # No silent row loss, duplicate matching, or trait/geography leaks.
    joined = frame.merge(extra, on="obs_id", how="left", validate="one_to_one",
                         indicator="_join_status")
    if (len(joined) != len(frame) or
            not joined._join_status.eq("both").all()):
        raise ValueError("EXTENDED_TABLE_MISSING_SOURCE_OBSERVATIONS")
    joined.drop(columns=["_join_status"], inplace=True)
    raw = joined[list(required_columns)].apply(pd.to_numeric, errors="coerce")
    for name in required_columns:
        # Detect erroneous strings and non-finite extraction values. NaN itself
        # is allowed; each tier subsequently uses complete cases without fill.
        old = joined[name]
        if ((old.notna() & raw[name].isna()).any() or
                np.isinf(raw[name].to_numpy(float)).any()):
            raise ValueError("EXTENDED_SOURCE_NONNUMERIC_OR_INFINITE: " + name)
        joined[name] = raw[name]
    coverage = {"schema": manifest["schema"],
                "table_sha256": table_sha256,
                "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
                "declaration_only_raw_raster_bytes_not_rechecked_here": True,
                "source_identifiers_not_exported": True,
                "source_rows": int(len(joined)),
                "actual_new_columns": list(required_columns),
                "nonmissing_counts_by_axis": {
                    name: int(joined[name].notna().sum()) for name in required_columns}}
    return joined, coverage


def _support_record(frame: pd.DataFrame, cols: tuple[str, ...]):
    x = frame[list(cols)].apply(pd.to_numeric, errors="coerce")
    complete = np.isfinite(x.to_numpy(dtype=float)).all(axis=1)
    return {"complete_environment_source_rows": int(complete.sum()),
            "environment_axis_count": len(cols)}


def run_extended(native: Path, receipt: Path, archived_seven: Path,
                 archived_seven_sha256: str, extra32: Path, extra32_sha256: str,
                 layer_manifest: Path, output: Path, *,
                 upper_tier: str = "complete_bioclim21",
                 trees: int = 160, repeats: int = 3,
                 constructs: tuple[str, ...] | None = None,
                 min_per_cell: int = 3,
                 splits: tuple[str, ...] = ("unseen_taxa", "unseen_longitude_sectors")) -> dict:
    """Run one nested comparison, or explicit posthoc support/split sensitivity.

    One tier per run avoids hiding a prohibitive 7x analysis compute bill, and
    keeps the missingness support for each contrast explicitly inspectable.
    """
    if upper_tier not in {higher for _, higher in PAIRS if higher != "archived12"}:
        raise ValueError("UNKNOWN_OR_BASE_TIER: " + upper_tier)
    lower_tier = dict((higher, lower) for lower, higher in PAIRS)[upper_tier]
    frame, archived_provenance = verified_source(
        native, receipt, archived_seven, archived_seven_sha256)
    frame, expanded_provenance = verify_extended_source(
        frame, extra32, extra32_sha256, layer_manifest,
        required_columns=added_columns(upper_tier))
    predictors = TIERS[upper_tier]
    paired_sets = {lower_tier: TIERS[lower_tier], upper_tier: predictors}
    requested = tuple(FEATURES) if constructs is None else tuple(constructs)
    allowed = {"unseen_taxa", "unseen_taxa_balanced", "unseen_longitude_sectors"}
    if (not requested or len(requested) != len(set(requested)) or
            not set(requested).issubset(FEATURES)):
        raise ValueError("INVALID_OR_DUPLICATED_REQUESTED_CONSTRUCTS")
    if (min_per_cell not in (1, 3) or not splits or len(set(splits)) != len(splits)
            or not set(splits).issubset(allowed)):
        raise ValueError("UNREGISTERED_POSTHOC_SENSITIVITY_DESIGN")
    results = []
    for construct in requested:
        try:
            cells = aggregate_endpoint(frame, construct, predictors=predictors,
                                       min_per_cell=min_per_cell)
        except ValueError as exc:
            # Support failure is a legitimate result, not grounds for silently
            # relaxing QC or changing the predictor set after seeing the data.
            if "INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS" not in str(exc):
                raise
            results.append({"construct": construct, "status": "insufficient_support",
                            "reason": "too_few_qualified_taxa_or_regions",
                            "both_model_cohorts_not_relaxed": True})
            continue
        for heldout in splits:
            if (heldout == "unseen_longitude_sectors" and
                    cells.groupby("region").taxon.nunique().min() < 5):
                results.append({"construct": construct, "generalisation": heldout,
                                "status": "insufficient_sector_taxon_support",
                                "min_taxa_per_heldout_sector": 5})
                continue
            out = score_models(cells, construct, heldout,
                               trees=trees, repeats=repeats,
                               predictor_sets=paired_sets,
                               predictor_groups=GROUPS44)
            out["status"] = "source_estimated_exploratory"
            results.append(out)
    report = {
        "schema": "azami_v4_exploratory_rf_extended44_v1",
        "status": "source_fitted_exploratory_not_original_Holm",
        "chosen_pair": [lower_tier, upper_tier],
        "requested_constructs": list(requested),
        "minimum_source_photos_per_taxon_cell": min_per_cell,
        "split_schemes": list(splits),
        "posthoc_support_or_folding_sensitivity": bool(
            min_per_cell != 3 or "unseen_taxa_balanced" in splits),
        "paired_predictor_universes": {n: list(p) for n, p in paired_sets.items()},
        "tier_data_support": {n: _support_record(frame, vars_) for n, vars_ in paired_sets.items()},
        "original_source": archived_provenance,
        "extension_source": expanded_provenance,
        "registered_44_candidate_axes": list(TIERS["vegetation44"]),
        "groups": {k: list(v) for k, v in GROUPS44.items()},
        "selection_scope": "Nested, source-declared, paired complete-case cohorts; no outcome-driven feature selection",
        "time_scope_warning": "1981-2010 climate, static soil/terrain, 2001-2010 NDVI and 2015 land cover are different temporal constructs",
        "prediction_interpretation": "out-of-group predictive associations; not cause, adaptation, humidity mechanism or phylogenetic convergence",
        "not_raw_raster_provenance_verification": True,
        "original_S1_S2_Holm_and_H3c_unmodified": True,
        "only_aggregate_metrics_exported": True,
        "results": results,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("native-table", "native-receipt", "archived-seven-table",
                "archived-seven-sha256", "extra32-table", "extra32-sha256",
                "layer-manifest", "output"):
        p.add_argument("--" + arg, type=Path if "sha256" not in arg else str,
                       required=True)
    p.add_argument("--upper-tier", choices=[i[1] for i in PAIRS if i[1] != "archived12"],
                   default="complete_bioclim21")
    p.add_argument("--trees", type=int, default=160)
    p.add_argument("--permutation-repeats", type=int, default=3)
    p.add_argument("--constructs", default=",".join(FEATURES))
    p.add_argument("--min-per-cell", type=int, choices=[1, 3], default=3)
    p.add_argument("--splits", default="unseen_taxa,unseen_longitude_sectors")
    args = p.parse_args()
    if args.trees < 8 or args.permutation_repeats < 1:
        p.error("at least 8 trees and one holdout permutation are required")
    run_extended(args.native_table, args.native_receipt,
                 args.archived_seven_table, args.archived_seven_sha256,
                 args.extra32_table, args.extra32_sha256,
                 args.layer_manifest, args.output, upper_tier=args.upper_tier,
                 trees=args.trees, repeats=args.permutation_repeats,
                 constructs=tuple(x.strip() for x in args.constructs.split(",") if x.strip()),
                 min_per_cell=args.min_per_cell,
                 splits=tuple(x.strip() for x in args.splits.split(",") if x.strip()))


if __name__ == "__main__":
    main()
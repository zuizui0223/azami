"""Source-gated exploratory WorldClim-terrain28 vs terrain28+SoilGrids-soil34 RF.

This module fits **only after** protected photo-coordinate raster sampling for
all six 0-5cm SoilGrids Q0.5 variables and audit of actual source TIFF bytes.
It does not obtain raster layers, invent soil, or change frozen GEB science.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_archived_atmosphere import verify_original_1981_2010_airfields
from analysis.v4.rf_extended_environment import (BIO9, TERRAIN5, SOIL6,
    verify_extended_source)
from analysis.v4.rf_multitrait_predictability import verified_source, score_models, aggregate_endpoint, sha256
from analysis.v4.rf_worldclim_terrain import MODEL_SET, GROUPS
from analysis.v4.rf_soilgrids_webdav_probe import CONVERSIONS, official_url
from analysis.v4.rf_soilgrids_private_pilot import VRT_PINNED_SHA256

TERRAIN28 = tuple(MODEL_SET["chelsa_vpd_wind23_plus_worldclim_terrain28"])
SOIL34 = TERRAIN28 + SOIL6
PAIRED = {"source_real_terrain28": TERRAIN28,
          "terrain28_plus_soilgrids_0to5cm_Q05_34": SOIL34}
PREDICTOR_GROUPS = {**GROUPS, "SoilGrids_0to5cm_Q05": SOIL6}
# Generic original-source provenance requires unit string >=3 characters.
# The official pH value itself is dimensionless: spell out water extraction.
SOIL_MANIFEST_UNITS = {
    name: ("pH_H2O" if name == "soil_phh2o_0_5cm" else CONVERSIONS[name][2])
    for name in SOIL6
}
CONSTRUCTS = ("angle", "chroma", "aspect")
SPLITS = ("unseen_taxa_balanced", "unseen_longitude_sectors")
assert len(SOIL34) == 34 and len(set(SOIL34)) == 34
assert len(TERRAIN28) == 28 and set(SOIL6).isdisjoint(TERRAIN28)


def verify_soil_source_manifest(path: Path, *, source_rows: int) -> dict:
    """Require actual-original-pixel source contracts, not a source VRT probe.

    This validates provenance declarations; the *private raster extractor* must
    independently have verified every used TIFF byte hash before writing it.
    """
    receipt = json.loads(Path(path).read_text(encoding="utf-8"))
    if (receipt.get("schema") != "azami_rf_extra32_source_manifest_v1" or
            set(receipt.get("layers", {})) != set(SOIL6) or
            receipt.get("source_observations") != source_rows or
            not receipt.get("original_S1_point_pixels_uninterpolated") or
            not receipt.get("all_used_tiff_tile_byte_sha_verified") or
            not receipt.get("raw_tiff_tile_receipt_private") or
            not receipt.get("original_pixel_nodata_preserved")):
        raise ValueError("SOILGRIDS_UNVERIFIED_FULL_ORIGINAL_PIXEL_MANIFEST")
    for col in SOIL6:
        entry = receipt["layers"][col]
        if (entry.get("dataset") != "SoilGrids-v2.0" or
                entry.get("temporal_support") != "static-soil-model" or
                entry.get("depth_interval") != "0-5cm" or
                entry.get("statistic") != "Q0.50" or
                entry.get("source_url") != official_url(col) or
                entry.get("source_bundle_sha256") != VRT_PINNED_SHA256[col] or
                entry.get("unit") != SOIL_MANIFEST_UNITS[col] or
                entry.get("conversion_factor_divisor") != CONVERSIONS[col][1] or
                not isinstance(entry.get("source_tiff_tiles_hashed"), int) or
                entry["source_tiff_tiles_hashed"] < 1 or
                not isinstance(entry.get("derivation"), str) or
                not isinstance(entry.get("sampling_method"), str)):
            raise ValueError("SOILGRIDS_UNREGISTERED_OR_UNVERIFIED_SOIL_RASTER: " + col)
    return {"manifest_sha256": sha256(Path(path)),
            "source_rows": source_rows,
            "soil6_source_vrt_pinned": True,
            "soil6_tiff_tile_hashes_declaration_checked_not_raw_bytes_rechecked_here": True,
            "no_photo_ids_coordinates_or_taxon_names_exported": True}


def run_soil34(native: Path, receipt: Path,
               original_seven: Path, seven_sha: str,
               bio_nine: Path, nine_sha: str, nine_manifest: Path,
               terrain_five: Path, terrain_manifest: Path,
               soil_six: Path, soil_manifest: Path,
               original_env_metadata: Path, output: Path, *,
               degrees: int=1, min_per_cell: int=3,
               constructs: tuple[str,...]=CONSTRUCTS,
               trees: int=80, repeats: int=2) -> dict:
    if (degrees not in (1,5) or min_per_cell != 3 or
            constructs != CONSTRUCTS or trees < 8 or repeats < 1):
        raise ValueError("UNDECLARED_SOIL34_EXPLORATORY_DESIGN")
    df, original = verified_source(native, receipt, original_seven, seven_sha)
    df, ext9 = verify_extended_source(df, bio_nine, nine_sha, nine_manifest,
                                      required_columns=BIO9)
    verify_original_1981_2010_airfields(df, receipt, original_env_metadata)
    df, terrain = verify_extended_source(df, terrain_five, sha256(terrain_five),
                                          terrain_manifest, required_columns=TERRAIN5)
    source_contract = verify_soil_source_manifest(soil_manifest, source_rows=len(df))
    df, soil = verify_extended_source(df, soil_six, sha256(soil_six),
                                      soil_manifest, required_columns=SOIL6)
    results = []
    for construct in CONSTRUCTS:
        try:
            cells = aggregate_endpoint(df, construct, predictors=SOIL34,
                                       degrees=degrees, min_per_cell=min_per_cell)
        except ValueError as exc:
            if "INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS" not in str(exc):
                raise
            results.append({"construct":construct, "status":"insufficient_qualified_taxa_or_three_regions",
                            "no_posthoc_cohort_relaxation":True})
            continue
        for split in SPLITS:
            if split=="unseen_longitude_sectors" and cells.groupby("region").taxon.nunique().min()<5:
                results.append({"construct":construct, "generalisation":split,
                                "status":"insufficient_five_taxa_per_heldout_sector"})
                continue
            result=score_models(cells,construct,split,trees=trees,repeats=repeats,
                                predictor_sets=PAIRED,predictor_groups=PREDICTOR_GROUPS)
            result["status"]="source_fitted_soil34_exploratory"
            results.append(result)
    out={"schema":"azami_v4_terrain28_SoilGrids6_soil34_predictability_v1",
         "status":"source_fitted_soil34_exploratory_not_original_GEB",
         "data":{"original_s1":original, "extra9":ext9,
                 "terrain5":terrain, "soil6":soil,
                 "soil_original_source_contract":source_contract},
         "predictors":{k:list(v) for k,v in PAIRED.items()},
         "design":{"degrees":degrees,"min_photos_per_cell":min_per_cell,
                   "min_photos_per_taxon":40,"min_cells_per_taxon":8,
                   "constructs":list(CONSTRUCTS),"split_targets":list(SPLITS),
                   "trees":trees,"permutation_repeats":repeats},
         "note":"SoilGrids climate/terrain-trained predictions are not independently measured ecological causes",
         "no_TerraClimate_7_or_soil39_assumption":True,
         "frozen_original_GEB_S1_S2_H3c_and_earlier_RF_unchanged":True,
         "only_aggregate_not_source_ids_coordinates_taxon_names_predictions_exported":True,
         "results":results}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("native", "receipt", "original-seven", "seven-sha",
                 "bio-nine", "nine-sha", "nine-manifest",
                 "terrain-five", "terrain-manifest", "soil-six", "soil-manifest",
                 "original-env-metadata", "output"):
        p.add_argument("--"+name, required=True,
                       type=str if name.endswith("-sha") else Path)
    p.add_argument("--degrees",type=int,choices=(1,5),default=1)
    a=p.parse_args()
    run_soil34(a.native,a.receipt,a.original_seven,a.seven_sha,
               a.bio_nine,a.nine_sha,a.nine_manifest,
               a.terrain_five,a.terrain_manifest,a.soil_six,a.soil_manifest,
               a.original_env_metadata,a.output,degrees=a.degrees)

if __name__=="__main__":main()
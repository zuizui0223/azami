"""Additional 1-degree geographic/taxon coverage sensitivity for terrain RF.

Consumes ONLY private S1 and previously checksum-verified 5 DEM-derived traits;
does not download or resample a DEM. Reuses the exact original 5-degree source
receipt and the same paired 23 vs 28 models, heldout taxa and longitude sectors.
The original manuscript and original 5-degree estimates are untouched.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_archived_atmosphere import verify_original_1981_2010_airfields
from analysis.v4.rf_extended_environment import BIO9, TERRAIN5, verify_extended_source
from analysis.v4.rf_multitrait_predictability import (
    aggregate_endpoint, score_models, sha256, verified_source,
)
from analysis.v4.rf_worldclim_terrain import GROUPS, MODEL_SET

CONSTRUCTS = ("angle", "chroma", "aspect")
SPLITS = ("unseen_taxa_balanced", "unseen_longitude_sectors")


def score_new_geographic_grain(frame, *, degrees=1, min_photos_per_cell=3,
                               constructs=CONSTRUCTS, trees=80, repeats=2):
    """Outcome-free grid rule frozen to 1°/3 photos before the RF rerun."""
    if (degrees != 1 or min_photos_per_cell != 3 or
            not constructs or len(set(constructs)) != len(constructs) or
            not set(constructs).issubset(FEATURES) or trees < 8 or repeats < 1):
        raise ValueError("UNDECLARED_1_DEGREE_RF_SENSITIVITY_SCOPE")
    result=[]
    for construct in constructs:
        try:
            cells = aggregate_endpoint(frame, construct,
                        predictors=MODEL_SET["chelsa_vpd_wind23_plus_worldclim_terrain28"],
                        degrees=1,min_per_cell=3)
        except ValueError as error:
            if "INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS" not in str(error):
                raise
            result.append({"construct":construct,"status":"insufficient_support"})
            continue
        for scheme in SPLITS:
            if scheme=="unseen_longitude_sectors" and cells.groupby("region").taxon.nunique().min()<5:
                result.append({"construct":construct,"generalisation":scheme,
                               "status":"insufficient_five_taxa_per_sector"})
                continue
            fit=score_models(cells,construct,scheme,trees=trees,repeats=repeats,
                             predictor_sets=MODEL_SET,predictor_groups=GROUPS)
            fit["status"]="posthoc_one_degree_more_taxa_not_original_GEB"
            result.append(fit)
    return result


def run(native:Path, receipt:Path, archived_seven:Path, archived_seven_sha256:str,
        bio_nine:Path,bio_nine_sha256:str,nine_manifest:Path,
        terrain_table:Path,terrain_manifest:Path,
        prior_five_degree:Path,original_metadata:Path,output:Path,
        *,trees:int=80,repeats:int=2):
    """Source-verified 1° sensitivity, same DEM derived in private 5° run."""
    reference=json.loads(prior_five_degree.read_text())
    if (reference.get("schema") !=
        "azami_original_s1_chelsa_atmosphere23_vs_worldclim_terrain28_exploratory_v1" or
        reference.get("status") != "original_source_terrain_derived_and_evaluated" or
        reference.get("design",{}).get("degrees") != 5 or
        reference.get("design",{}).get("min_per_cell") != 3 or
        not reference.get("frozen_original_GEB_S1_S2_H3c_unchanged")):
        raise ValueError("SOURCE_VERIFIED_FIVE_DEGREE_PRIMARY_EXPLORATORY_RECEIPT_REQUIRED")
    provenance=reference["terrain_provenance"]
    if (sha256(terrain_table) != provenance.get("terrain_private_table_sha256") or
        sha256(terrain_manifest) != provenance.get("terrain_private_manifest_sha256") or
        provenance.get("mismatch_n") != 0 or
        provenance.get("matching_original_topo_elevation_records",0) < 250000):
        raise ValueError("ORIGINAL_WORLDCLIM_TERRAIN_SHA_OR_PIXEL_IDENTITY_DRIFT")
    frame,source=verified_source(native,receipt,archived_seven,archived_seven_sha256)
    if source["source_s1_sha256"] != reference["source"]["source_s1_sha256"]:
        raise ValueError("S1_SOURCE_HASH_NOT_PRIOR_TERRAIN_RUN")
    frame,_=verify_extended_source(frame,bio_nine,bio_nine_sha256,nine_manifest,
                                    required_columns=BIO9)
    verify_original_1981_2010_airfields(frame,receipt,original_metadata)
    frame,terrain=verify_extended_source(frame,terrain_table,sha256(terrain_table),
                                         terrain_manifest,required_columns=TERRAIN5)
    results=score_new_geographic_grain(frame,trees=trees,repeats=repeats)
    out={"schema":"azami_S1_posthoc_one_degree_terrain23_vs28_sensitivity_v1",
         "status":"one_degree_source_authenticated_sensitivity_not_primary_GEB",
         "source":{"s1_observations":source["s1_observations"],
                   "source_s1_sha256":source["source_s1_sha256"]},
         "terrain_source":{"dem_tif_sha256":provenance["dem_tif_sha256"],
                           "DEM_rechecked_against_original_elevation_mismatch_n":0,
                           "private_terrain_table_sha256":terrain["table_sha256"]},
         "design":{"cell_degrees":1,"min_photos_per_cell":3,
                   "min_photos_per_taxon":40,"min_cells_per_taxon":8,
                   "constructs":list(CONSTRUCTS),"splits":list(SPLITS),
                   "comparison_models":{k:list(v) for k,v in MODEL_SET.items()},
                   "trees":trees,"permutation_repeats":repeats},
         "sampling_gates_chosen_after_original_5_degree_RF":True,
         "results":results,
         "ecological_scope":"paired out-of-taxon and out-of-sector predictive association, not causal terrain adaptation",
         "original_5_degree_RF_and_frozen_GEB_S1_S2_H3c_unchanged":True,
         "only_aggregates_not_raw_photo_obs_taxon_coords_or_predictions_exported":True}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("native","receipt","archived-seven","archived-seven-sha256",
                 "bio-nine","bio-nine-sha256","nine-manifest",
                 "terrain-table","terrain-manifest","prior-five-degree",
                 "original-metadata","output"):
        p.add_argument("--"+name,type=str if name.endswith("sha256") else Path,required=True)
    p.add_argument("--trees",type=int,default=80)
    p.add_argument("--repeats",type=int,default=2)
    a=p.parse_args()
    run(a.native,a.receipt,a.archived_seven,a.archived_seven_sha256,
        a.bio_nine,a.bio_nine_sha256,a.nine_manifest,a.terrain_table,
        a.terrain_manifest,a.prior_five_degree,a.original_metadata,a.output,
        trees=a.trees,repeats=a.repeats)

if __name__=="__main__":
    main()
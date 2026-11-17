"""Source-real, aggregate-only geography/photographic support diagnostic.

A diagnostic of eligible taxa, NOT a posthoc refit of frozen GEB or proof that
relaxed photo-cell thresholds are ecologically preferable. No taxon or photo
identifiers are emitted. This is independent of the DEM fitness test.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_archived_atmosphere import ATMOSPHERE23, verify_original_1981_2010_airfields
from analysis.v4.rf_extended_environment import BIO9, verify_extended_source
from analysis.v4.rf_multitrait_predictability import (
    aggregate_endpoint, verified_source,
)

GRAINS = (1, 2, 5)
MIN_PHOTOS_PER_CELL = (1, 3)
CONSTRUCTS = ("angle", "chroma", "aspect")


def support_lattice(frame, *, predictors=ATMOSPHERE23, constructs=CONSTRUCTS):
    """Freeze outcome-independent thresholds and output aggregate support only."""
    if not constructs or not set(constructs).issubset(FEATURES):
        raise ValueError("SUPPORT_ATLAS_CONSTRUCTS_INVALID")
    rows=[]
    for name in constructs:
        for deg in GRAINS:
            for minimum in MIN_PHOTOS_PER_CELL:
                try:
                    tab=aggregate_endpoint(frame,name,predictors=predictors,
                                          min_taxa=1,min_per_cell=minimum,
                                          min_photos=40,min_cells=8,degrees=deg)
                except ValueError as exc:
                    if "INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS" not in str(exc):
                        raise
                    rows.append({"construct":name,"cell_degrees":deg,
                                "minimum_photos_per_cell":minimum,
                                "status":"insufficient_three_sector_support",
                                "n_taxa":0,"n_taxon_cells":0,"n_contributing_photos":0})
                    continue
                sector=tab.groupby("region").taxon.nunique().to_dict()
                rows.append({"construct":name,"cell_degrees":deg,
                            "minimum_photos_per_cell":minimum,
                            "min_total_photos_per_taxon":40,
                            "min_geographic_cells_per_taxon":8,
                            "status":"eligible_counts_only",
                            "n_taxa":int(tab.taxon.nunique()),
                            "n_taxon_cells":int(len(tab)),
                            "n_contributing_photos":int(tab.photos.sum()),
                            "sector_taxa_count":{str(k):int(sector.get(k,0)) for k in range(3)},
                            "eligible_for_grouped_region_holdout":bool(min(sector.get(k,0) for k in range(3))>=5),
                            "no_outcome_prediction_run":True})
    return {"schema":"azami_rf_S1_cross_grid_coverage_sensitivity_v1",
            "status":"source_eligibility_counts_only_no_selected_grid",
            "registered_grains":list(GRAINS),"minimum_photos_per_cell":list(MIN_PHOTOS_PER_CELL),
            "constructs":list(constructs),"n_environment_variables":len(predictors),
            "posthoc_sampling_support_diagnostic":True,
            "no_source_ids_coordinates_trait_means_or_predictions_exported":True,
            "records":rows}


def run(native: Path, receipt: Path, seven: Path, seven_sha: str,
        nine: Path, nine_sha: str, nine_manifest: Path,
        metadata: Path, output: Path):
    frame,origin=verified_source(native,receipt,seven,seven_sha)
    frame,_=verify_extended_source(frame,nine,nine_sha,nine_manifest,required_columns=BIO9)
    air=verify_original_1981_2010_airfields(frame,receipt,metadata)
    result=support_lattice(frame)
    result.update(original_source={"s1_observations":origin["s1_observations"],
                                   "source_s1_sha256":origin["source_s1_sha256"]},
                  atmospheric_verified_source_hashes={k:v["source_sha256"] for k,v in air.items()},
                  original_GEB_frozen_results_unchanged=True)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("native","receipt","archived-seven","archived-seven-sha256",
                "bio-nine","bio-nine-sha256","nine-manifest","metadata","output"):
        p.add_argument("--"+arg,type=str if arg.endswith("sha256") else Path,required=True)
    a=p.parse_args()
    run(a.native,a.receipt,a.archived_seven,a.archived_seven_sha256,
        a.bio_nine,a.bio_nine_sha256,a.nine_manifest,a.metadata,a.output)

if __name__=="__main__":
    main()
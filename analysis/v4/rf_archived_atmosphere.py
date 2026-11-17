"""Source-audited post-BIO21 RF add-on: existing CHELSA VPD and wind fields.

Unlike a new soil or vegetation covariate, these two 1981--2010 atmospheric
fields were already measured and SHA/scale-decoded in the immutable native S1
archive. Their addition is exploratory and must not modify frozen GEB models.
Only aggregate held-out skill, fold counts and domain permutation loss leave
this runner. No original observation/species/photographer/geographic rows leave.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_extended_environment import (
    BIO9, GROUPS44, TIERS, verify_extended_source,
)
from analysis.v4.rf_multitrait_predictability import (
    aggregate_endpoint, score_models, verified_source,
)

ATMOSPHERE2 = ("chelsa_vpd_mean", "chelsa_sfcwind_mean")
BASE21 = tuple(TIERS["complete_bioclim21"])
ATMOSPHERE23 = BASE21 + ATMOSPHERE2
MODEL_PAIR = {"complete_bioclim21": BASE21,
              "bio21_plus_chelsa_vpd_wind23": ATMOSPHERE23}
PREDICTOR_GROUPS = {**GROUPS44,
                    "atmospheric_demand_CHELSA": (ATMOSPHERE2[0],),
                    "wind": (*GROUPS44["wind"], ATMOSPHERE2[1])}
assert len(ATMOSPHERE23)==23 and len(set(ATMOSPHERE23))==23


def verify_original_1981_2010_airfields(frame: pd.DataFrame,
                                 receipt_path: Path, metadata_path: Path) -> dict:
    """Require the exact original nine-field CHELSA decoding receipts.

    This does not invent source units: VPD/wind units follow the archived
    source raster contract, with original scale and SHA kept explicitly.
    """
    receipt=json.loads(Path(receipt_path).read_text(encoding="utf-8"))
    metadata=json.loads(Path(metadata_path).read_text(encoding="utf-8"))
    decoding=receipt.get("environment_decoding")
    if (not isinstance(decoding,dict) or
            receipt.get("environment_decoding_status") !=
            "physical_units_gsp_invalid_encoding_masked_v1"):
        raise ValueError("ORIGINAL_NINE_ENVIRONMENT_DECODE_NOT_VERIFIED")
    coverage={}
    for name in ATMOSPHERE2:
        if name not in frame or name not in metadata or name not in decoding:
            raise ValueError("ARCHIVED_ATMOSPHERIC_AXIS_NOT_PRESENT_"+name)
        source=metadata[name]
        saved=decoding[name]
        if (saved.get("source_sha256") != source.get("source_sha256") or
                saved.get("scale") != source.get("scale") or
                saved.get("offset") != source.get("offset")):
            raise ValueError("ARCHIVED_ATMOSPHERIC_SOURCE_SHA_SCALE_DRIFT_"+name)
        values=pd.to_numeric(frame[name],errors="coerce").to_numpy(float)
        if np.isinf(values).any():
            raise ValueError("ARCHIVED_ATMOSPHERIC_AXIS_INFINITE_"+name)
        n=int(np.isfinite(values).sum())
        if n < 0.75*len(frame):
            raise ValueError("ARCHIVED_ATMOSPHERIC_AXIS_LOW_SUPPORT_"+name)
        coverage[name] = {
            "source_sha256": source["source_sha256"],
            "source_stored_scale":source["scale"],
            "source_stored_offset":source["offset"],
            "original_source_variable":name,
            "environment_original_period":"1981-2010",
            "finite_s1_records":n,
            "nonmissing_fraction":n/len(frame),
            "biologically_interpreted_units_not_reinferred":True,
        }
    return coverage


def run_atmosphere(native: Path, receipt: Path, archived_seven: Path,
                   archived_seven_sha256: str, extra_nine: Path,
                   extra_nine_sha256: str, nine_manifest: Path,
                   environment_metadata: Path, output: Path, *,
                   constructs: tuple[str,...]=tuple(FEATURES),
                   min_per_cell: int=3,
                   splits: tuple[str,...] = ("unseen_taxa", "unseen_longitude_sectors"),
                   trees: int=80, repeats:int=3) -> dict:
    """21 vs 23 paired, same photos/cells within each construct and split."""
    if (not constructs or len(set(constructs))!=len(constructs) or
            not set(constructs) <= set(FEATURES) or
            min_per_cell not in (1,3) or
            not splits or len(set(splits))!=len(splits) or
            not set(splits) <= {"unseen_taxa", "unseen_taxa_balanced",
                                "unseen_longitude_sectors"} or
            trees < 8 or repeats < 1):
        raise ValueError("INVALID_PREDECLARED_ATMOSPHERE_COMPARISON_SCOPE")
    frame, original = verified_source(native,receipt,archived_seven,
                                      archived_seven_sha256)
    frame, extension=verify_extended_source(frame,extra_nine,
                                             extra_nine_sha256,nine_manifest,
                                             required_columns=BIO9)
    air=verify_original_1981_2010_airfields(frame,receipt,environment_metadata)
    results=[]
    for construct in constructs:
        try:
            cells=aggregate_endpoint(frame,construct,predictors=ATMOSPHERE23,
                                     min_per_cell=min_per_cell)
        except ValueError as exc:
            if "INSUFFICIENT_QUALIFIED_TAXA_OR_THREE_REGIONS" not in str(exc):
                raise
            results.append({"construct":construct,"status":"insufficient_source_support"})
            continue
        for split in splits:
            if (split=="unseen_longitude_sectors" and
                    cells.groupby("region").taxon.nunique().min()<5):
                results.append({"construct":construct,"generalisation":split,
                                "status":"insufficient_sector_taxa"})
                continue
            x=score_models(cells,construct,split,trees=trees,repeats=repeats,
                           predictor_sets=MODEL_PAIR,
                           predictor_groups=PREDICTOR_GROUPS)
            x["status"]="original_s1_atmosphere23_exploratory"
            results.append(x)
    report={
        "schema":"azami_v4_original_s1_chelsa_vpd_wind_bio21_to23_v1",
        "status":"source_fitted_exploratory_not_original_S1_S2_Holm",
        "source":original,
        "bio21_extension":extension,
        "archived_atmospheric_variables":air,
        "comparison_models":{name:list(p) for name,p in MODEL_PAIR.items()},
        "requested_constructs":list(constructs),
        "minimum_source_photos_per_taxon_cell":min_per_cell,
        "split_schemes":list(splits),
        "posthoc_support_or_folding_sensitivity":bool(min_per_cell!=3 or
                                                     "unseen_taxa_balanced" in splits),
        "n_estimators":trees,
        "only_climate_1981_2010_archived_sources":True,
        "no_new_external_raster_downloads_required_for_VPD_or_wind":True,
        "evidence":"Source-authenticated out-of-group predictability; not original frozen regression or causal response",
        "original_S1_S2_Holm_and_H3c_unmodified":True,
        "no_raw_source_record_taxon_coordinates_or_model_weights_exported":True,
        "results":results,
    }
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n",
                      encoding="utf-8")
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("native-table","native-receipt","archived-seven-table",
                 "archived-seven-sha256","extra-nine-table","extra-nine-sha256",
                 "nine-manifest","environment-metadata","output"):
        parser.add_argument("--"+name,type=str if "sha256" in name else Path,required=True)
    parser.add_argument("--constructs",default=",".join(FEATURES))
    parser.add_argument("--min-per-cell",type=int,choices=[1,3],default=3)
    parser.add_argument("--splits",default="unseen_taxa,unseen_longitude_sectors")
    parser.add_argument("--trees",type=int,default=80)
    parser.add_argument("--permutation-repeats",type=int,default=3)
    args=parser.parse_args()
    run_atmosphere(args.native_table,args.native_receipt,args.archived_seven_table,
                   args.archived_seven_sha256,args.extra_nine_table,args.extra_nine_sha256,
                   args.nine_manifest,args.environment_metadata,args.output,
                   constructs=tuple(x.strip() for x in args.constructs.split(",") if x.strip()),
                   min_per_cell=args.min_per_cell,
                   splits=tuple(x.strip() for x in args.splits.split(",") if x.strip()),
                   trees=args.trees,repeats=args.permutation_repeats)

if __name__=="__main__":
    main()
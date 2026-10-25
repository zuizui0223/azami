"""All-Cirsium native/introduced range and orientation-phenotype support census.

Computes, BEFORE inspecting any orientation angle or climatic association, the
number of original species-rank names with both WCVP-native and EXPLICIT
WCVP-introduced wild/source observations, archived orientation availability,
independent 5-degree equal-area cells and original human observers.

One source species name is one possible analysis unit ONLY if its accepted
WCVP key is uniquely linked to that one source name. Synonym/concept
collisions are quarantined, not silently pooled as within-species variation.

Source level IDs, coordinates, users, observation or photo IDs are never
exported; output contains only aggregate taxon-level counts and continent
support where privacy-reasonable source size is available. No climatic slope
is fitted and no unmeasured image is imputed as a straight head.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from analysis.v4.native_introduced_opportunity import (
    NATIVE_N, NATIVE_ACCEPTED, SOURCE_N, CONTINENTS, angular_cells,
)

PRE_SCREEN={
    "min_archived_orientation_per_range":30,
    "min_5degree_equal_area_cells_per_range":5,
    "min_original_observers_per_range":10,
}
STUDY_SCREEN={
    "min_archived_orientation_per_range":100,
    "min_5degree_equal_area_cells_per_range":15,
    "min_original_observers_per_range":10,
}
CONTINENT_PUBLIC_MIN_RECORDS=20
CONTINENT_PUBLIC_MIN_OBSERVERS=5
EXPECTED_PHOTO_ID_UNIT="source species-rank observation (not individual plant or independent head)"


def records_by_source_species(frame:pd.DataFrame)->dict:
    required={
        "obs_id","taxon_name","taxon_rank","accepted_key",
        "native_primary","eligibility_reason",
        "orientation_eligible","colour_eligible",
        "has_archived_measurement_row","user_id",
        "latitude","longitude","tdwg1","source_photo_count",
        "quality_grade",
    }
    if not required <= set(frame):
        raise ValueError("Original full-source eligibility/observer fields absent")
    if len(frame)!=SOURCE_N or not frame.obs_id.is_unique:
        raise ValueError("Original 665115 observation census changed")
    native=frame.native_primary.eq(True)
    intro=frame.eligibility_reason.eq("range_introduced")
    if native.sum()!=NATIVE_N or (native&intro).any() or not intro.any():
        raise ValueError("Source wild native/introduced partition no longer matches")
    eligible=frame.loc[(native|intro)&frame.taxon_rank.eq("species")].copy()
    if eligible.accepted_key.isna().any():
        raise ValueError("Source species accepted-key determination incomplete")
    if (eligible.user_id.isna().any() or eligible.latitude.isna().any() or
        eligible.longitude.isna().any()):
        raise ValueError("Original observer/geocoordinates absent on accepted wild species")
    eligible["accepted_key"]=eligible.accepted_key.astype(str)
    eligible["_range"]=np.where(eligible.native_primary,"native","introduced")
    eligible["_cell"]=angular_cells(eligible.latitude,eligible.longitude)
    if eligible._cell.isna().any():
        raise ValueError("Equal-sin source cell calculation failed")
    # Taxonomic clashes are assessed on all accepted wild data, not just the
    # native fraction, because an invasion scenario must not pool two names.
    key_to_names=eligible.groupby("accepted_key",sort=True).taxon_name.nunique()
    collision_keys=set(key_to_names[key_to_names.gt(1)].index)
    # This isn't a tree; name resolution can be a synonym, annotation error,
    # or a legitimate rank conflict. Quarantine before molecular conclusions.
    source_names=sorted(frame.loc[frame.taxon_rank.eq("species"),
                                  "taxon_name"].dropna().astype(str).unique())
    all_source=frame.loc[frame.taxon_rank.eq("species")]
    observed_wild=eligible.taxon_name.value_counts()
    output=[]
    by_continent=[]
    for name in source_names:
        allname=all_source.loc[all_source.taxon_name.eq(name)]
        wild=eligible.loc[eligible.taxon_name.eq(name)]
        keys=set(wild.accepted_key)
        if len(keys)>1:
            raise ValueError("One source species mapped to multiple accepted WCVP keys")
        accepted=next(iter(keys)) if keys else None
        boundary=("ACCEPTED_KEY_COLLISION" if accepted in collision_keys else
                  "RESOLVED_SINGLE_SOURCE_NAME" if accepted is not None else
                  "NO_ACCEPTED_WILD_SOURCE_RECORDS")
        one={"source_taxon_name":name,
             "WCVP_accepted_key":accepted,
             "taxon_concept":boundary,
             "original_all_range_observations":len(allname),
             "original_other_unresolved_or_captive_observations":len(allname)-len(wild)}
        both={}
        for state in ("native","introduced"):
            subset=wild.loc[wild._range.eq(state)]
            matched=subset.loc[subset.orientation_eligible]
            counts={
                "original_source_observations":len(subset),
                "source_photo_links":int(subset.source_photo_count.sum()),
                "original_observers_source":int(subset.user_id.nunique()),
                "original_source_5deg_cells":int(subset._cell.nunique()),
                "archived_image_measurement_rows":int(subset.has_archived_measurement_row.sum()),
                "orientation_eligible_observations":len(matched),
                "orientation_eligible_original_observers":int(matched.user_id.nunique()),
                "orientation_eligible_5deg_cells":int(matched._cell.nunique()),
                "colour_eligible_observations":int(subset.colour_eligible.sum()),
                "orientation_research_grade":int(matched.quality_grade.eq("research").sum())
            }
            both[state]=counts
            for continent,regional in matched.groupby("tdwg1",sort=True):
                if pd.isna(continent):
                    raise ValueError("Original WCVP wild source continent missing")
                group=int(continent)
                if group not in CONTINENTS:
                    raise ValueError("Unexpected WCVP Level-1 region")
                if len(regional) < CONTINENT_PUBLIC_MIN_RECORDS or regional.user_id.nunique() < CONTINENT_PUBLIC_MIN_OBSERVERS:
                    # Explicitly count source-suppressed groups, never publish
                    # a 1-observation continent cell.
                    continue
                by_continent.append({
                    "source_taxon_name":name,
                    "range_state":state,
                    "TDWG_lvl1_region":CONTINENTS[group],
                    "orientation_eligible_observations":int(len(regional)),
                    "orientation_original_observers":int(regional.user_id.nunique()),
                    "orientation_equal_sin_cells":int(regional._cell.nunique())
                })
        one["range_support"]=both
        def qualifies(threshold):
            return boundary=="RESOLVED_SINGLE_SOURCE_NAME" and all(
                both[r]["orientation_eligible_observations"]>=threshold["min_archived_orientation_per_range"]
                and both[r]["orientation_eligible_5deg_cells"]>=threshold["min_5degree_equal_area_cells_per_range"]
                and both[r]["orientation_eligible_original_observers"]>=threshold["min_original_observers_per_range"]
                for r in ("native","introduced")
            )
        one["both_ranges_pass_basic_measurement_support"]=qualifies(PRE_SCREEN)
        one["both_ranges_pass_stronger_cell_support"]=qualifies(STUDY_SCREEN)
        output.append(one)
    if len(output)!=len(source_names):
        raise AssertionError("Source full species-rank roster lost a name")
    if sum(row["range_support"]["native"]["original_source_observations"] for row in output)!=int(
       ((frame.taxon_rank=="species")&native).sum()):
        raise ValueError("Species-native source tally lost records")
    if sum(row["range_support"]["introduced"]["original_source_observations"] for row in output)!=int(
       ((frame.taxon_rank=="species")&intro).sum()):
        raise ValueError("Species-introduced source tally lost records")
    screened=[r for r in output if r["both_ranges_pass_basic_measurement_support"]]
    better=[r for r in output if r["both_ranges_pass_stronger_cell_support"]]
    return {
        "schema_version":"v4_all_source_species_native_vs_introduced_orientation_support_v1",
        "status":"EXECUTED_WCVP_ALL_SPECIES_MEASUREMENT_AVAILABILITY_NOT_CLIMATE_EFFECT",
        "full_original_source_observations":SOURCE_N,
        "native_original_observations":NATIVE_N,
        "introduced_original_observations":int(intro.sum()),
        "original_species_rank_observations":int((frame.taxon_rank=="species").sum()),
        "source_species_rank_names":len(source_names),
        "source_species_with_any_accepted_wild_record":int(sum(
            r["WCVP_accepted_key"] is not None for r in output)),
        "source_species_with_both_wild_regions":int(sum(
            r["range_support"]["native"]["original_source_observations"]>0 and
            r["range_support"]["introduced"]["original_source_observations"]>0
            for r in output)),
        "mixed_accepted_key_count_in_all_wild_source_species":len(collision_keys),
        "frozen_native_collision_count_must_not_be_reused_for_full_range":True,
        "predeclared_basic_support":PRE_SCREEN,
        "predeclared_strong_support":STUDY_SCREEN,
        "species_passing_basic_observer_cell_measurement_support":len(screened),
        "species_passing_stronger_observer_cell_measurement_support":len(better),
        "contiguous_source_species_counts":output,
        "continent_coverage_for_measurement_supported_groups":by_continent,
        "suppressed_low_support_species_x_continent_entries":True,
        "source_observation_photo_id_unit":EXPECTED_PHOTO_ID_UNIT,
        "claim_boundary":[
            "Species-rank originals are author identifications; WCVP accepted keys can merge different names, which are excluded from multilevel inference.",
            "Source image count or archived orientation availability is not independent flowering-head identity, native/introduced genetic origin or measurement accuracy.",
            "This is phenotype-availability selection only; it does not establish original CHELSA pixel availability in introduced-only cells.",
            "Climate common-support and camera quality must be audited for every candidate before cross-range effect estimation.",
            "Original observers, coordinates and observation IDs remain private; not released in this aggregate."
        ]
    }

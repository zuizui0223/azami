"""EAzami Japan38 exact-name candidate joining preserves original-source exclusions."""
from __future__ import annotations

import csv
from pathlib import Path
import pytest
from analysis.v4.t2_japan38_all_s1_nuclear_bm_ou import join_original_244
from analysis.v4.t2_full_species_bm_ou import FIXED_S1


def test_japan38_uses_original_taxon_rank_and_no_duplicate_tip(tmp_path: Path):
    path=tmp_path/"s1_244.csv"
    rows=[]
    for i in range(244):
        source=("Cirsium alpicola" if i==0 else
                "Cirsium japonicum var. horridum" if i==1 else
                "Cirsium dummy_species_"+str(i))
        row={
            "accepted_key":str(i),
            "normalized_exact_source_name_candidate":source,
            "source_concept_collision":"False",
            "eligible_photo_observations":170066-243 if i==0 else 1,
            "S1_mean_EXIF_angle_deg":88.0,
            "S1_median_EXIF_angle_deg":88.0,
        }
        row.update({name+"_photo_eligible_mean":float(k+1) for k,name in enumerate(FIXED_S1)})
        rows.append(row)
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    cmap=[
      {"paper_japan_member_id":"JPN_03","paper_taxon_concept":"Cirsium alpicola Nakai",
       "tip_ids":"J38S003","n_biological_samples":"1","trait_asr_primary_allowed":"true"},
      {"paper_japan_member_id":"JPN_12","paper_taxon_concept":"Cirsium japonicum Fisch. ex DC. var. horridum Nakai",
       "tip_ids":"J38S012","n_biological_samples":"1","trait_asr_primary_allowed":"true"},
      {"paper_japan_member_id":"JPN_20","paper_taxon_concept":"Cirsium nipponicum var. incomptum",
       "tip_ids":"J38S020|J38S021","n_biological_samples":"2","trait_asr_primary_allowed":"true"},
      {"paper_japan_member_id":"JPN_31","paper_taxon_concept":"Cirsium yuki-uenoanum",
       "tip_ids":"J38S032","n_biological_samples":"1","trait_asr_primary_allowed":"false"}
    ]
    used,ledger=join_original_244(path,cmap)
    assert len(used)==2
    assert {r["EAzami_JPN_concept_id"] for r in used}=={"JPN_03","JPN_12"}
    assert {r["S1_accepted_key"] for r in used}=={"0","1"}
    assert dict((x["EAzami_JPN_concept_id"],x["status"]) for x in ledger)["JPN_20"]=="EAzami_nonmonophyletic_replications_not_collapsed"
    assert dict((x["EAzami_JPN_concept_id"],x["status"]) for x in ledger)["JPN_31"]=="EAzami_disallowed_primary_identity_or_concept"


def test_nonmonophyletic_and_taxon_identity_are_not_silently_accepted(tmp_path: Path):
    # Full-244 source cohort is a prerequisite. Fewer rows must never be
    # promoted into a "completed Japan38 nuclear phylogenomic comparison".
    path=tmp_path/"incomplete.csv"
    path.write_text("accepted_key,source_concept_collision,eligible_photo_observations\n1,False,170066\n")
    with pytest.raises(ValueError,match="Authentic all-244"):
        join_original_244(path,[])

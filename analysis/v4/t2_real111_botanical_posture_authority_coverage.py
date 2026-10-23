"""Independent botanical capitulum posture coverage for all 251 original S1 names.

Join original photographs (244 accepted keys, 251 exact source-name units,
170066 images) and exact 111 Moreyra public nuclear sample names with
independently scored, source-authenticated U/D head positions:
Siniscalchi2023 NA63; EAzami NMNS Japan38; EAzami eastern20 exact original
descriptions/treatments; held-out Flora of China source-evidence seed.

Taxon-concept authority scoring is NOT same-voucher phenotype or
gravity-referenced original photographic angle, nor does it certify the
unresolved 111-tip nucleotide tree, independent origins, or rain adaptation.

Keep every original photo unit (including seven collapsed-key pairs). Exact
taxon/rank only: do not strip infraspecific ranks or accept a binomial
as a score for an unspecified variety. Contradictory sources remain ?.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

from analysis.v4.t2_na63_posture_real_sra_sampling_gate import download_published_states

ORIGINAL_UNITS=251
ORIGINAL_KEYS=244
ORIGINAL_IMAGES=170066
MEASURED_NUCLEAR=111
MIXED_KEYS=7
NA_BOTANICAL=63
NA_NODDING=5

EAZAMI_COMMIT="c9870e91604ed1514344ddc8f765eae5237dbf27"
SOURCES={
 "nmns_japan38_seed":(
   "data/evidence/japan38_nmns_capitulum_trait_seed_v1.csv",
   "c4da867491e913119d78e4b07df800965a0f021d",22),
 "nmns_japan38_extension":(
   "data/evidence/japan38_nmns_capitulum_trait_seed_extension_v2.csv",
   "346381eeb1326c654cec22888467eb446ab058dd",2),
 "eazami_eastasia20_treatment":(
   "data/evidence/orientation_comp1061_20tip_source_crosswalk_v1.csv",
   "f655544e77f1659323c1e915567dacbb51b01cc3",20),
 "heldout_flora_of_china":(
   "data/evidence/heldout_east_asia_orientation_candidate_panel_v2.csv",
   "2ff3e36c0c195b99a3053446b6e38ba76191ae3e",13),
}
NA_SOURCE_URL=(
  "https://github.com/carol-siniscalchi/cirsium/blob/"
  "c925d4af86160064448f794a9715f70b917a2786/traits/traits.csv"
)


def normalize(s):
    return " ".join(str(s or "").replace("_"," ").strip().split())


def git_blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\x00"+data).hexdigest()


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        records=list(csv.DictReader(f))
    if not records:
        raise ValueError("Original external taxonomy or S1 source has no records")
    return records


def decode_state(raw):
    state=str(raw or "").lower()
    if state in ("0","u","upward_or_erect","upward_or_ascending"):
        return "U"
    if state in ("1","d","downward_or_nodding"):
        return "D"
    if state in ("?","unknown","source_conflict",
                 "source_conflict_index_downward_detail_erect",
                 "lateral_or_slightly_downward_nonbinary",
                 "erect_or_nodding",""):
        return "?"
    raise ValueError("Unrecognized original botanical posture ontology "+str(raw))


def add(index,taxon,state,source,url,claim,raw):
    t=normalize(taxon)
    if not t.startswith("Cirsium "):
        raise ValueError("Botanical source row is not a Cirsium taxon")
    index[t].append({
      "state":decode_state(state),"source":source,
      "authority_url":str(url),"source_claim_boundary":str(claim),
      "source_state":str(raw),
      "taxon_concept_not_sample_voucher":True,
    })


def load_authorities(root,na_states=None):
    if na_states is None:
        na_states=download_published_states()
    if len(na_states)!=NA_BOTANICAL or sum(na_states.values())!=NA_NODDING:
        raise ValueError("Original NA63 botanical authority scope changed")
    index=defaultdict(list)
    for taxon,state in na_states.items():
        add(index,taxon,str(state),"Siniscalchi2023_NA63_HeadPosition",
            NA_SOURCE_URL,"original author state; no identity with Moreyra voucher",
            str(state))
    sha={}
    root=Path(root)
    for category,(rel,expected,nrows) in SOURCES.items():
        path=root/rel
        raw=path.read_bytes()
        sha[category]=git_blob(raw)
        if sha[category]!=expected:
            raise ValueError(f"Pinned EAzami botanical source Git blob mismatch: {category}")
        rows=read_csv(path)
        if len(rows)!=nrows:
            raise ValueError("Original botanical row count drift at "+category)
        for row in rows:
            if category.startswith("nmns"):
                if row["authority_match_status"]!="exact_authority_concept_match":
                    raise ValueError("Unexpected NMNS concept match not exact")
                add(index,row["nmns_taxon_concept"],row["orientation_state"],
                    category,row["source_url"],row["claim_boundary"],row["orientation_state"])
            elif category=="eazami_eastasia20_treatment":
                if str(row["same_primary_voucher_state_verified"]).lower()!="false":
                    raise ValueError("EAzami20 source suddenly claims same voucher without review")
                add(index,row["accepted_taxon"],row["analysis_state"],category,
                    row["source_url"],row["claim_boundary"],row["state_label"])
            else:
                if row["source_lane"]!="Flora_of_China" or row["eligible_for_occurrence_gate"]!="true":
                    raise ValueError("Held-out botanical source not ready or lane drift")
                add(index,row["analysis_taxon"],row["orientation_state"],category,
                    row["source_url"],row["source_basis"],row["orientation_state"])
    return index,sha


def infer(index,taxon):
    """Require unanimous assigned source states; never majority-vote disputes."""
    sources=index.get(normalize(taxon),[])
    states={r["state"] for r in sources if r["state"] in ("U","D")}
    unknown=any(r["state"]=="?" for r in sources)
    if len(states)>1:
        label="?"
        status="conflicting_independent_botanical_sources"
    elif not states:
        label="?"
        status="unscored" if not sources else "indeterminate_or_ambiguous_authority"
    else:
        label=next(iter(states))
        status="authority_taxon_concept_only_ambiguity_also_flagged" if unknown else "authority_taxon_concept_only"
    return {
        "botanical_U_D_state":label,
        "botanical_state_status":status,
        "n_source_documents":len(sources),
        "any_explicit_unresolved_source":unknown,
        "botanical_source_names":"|".join(sorted(set(x["source"] for x in sources))),
        "botanical_source_urls":"|".join(sorted(set(x["authority_url"] for x in sources))),
        "botanical_U_D_same_public_voucher_verified":False,
    }


def join(original,nuclear,index):
    if (len(original)!=ORIGINAL_UNITS or
        len({x["accepted_key"] for x in original})!=ORIGINAL_KEYS or
        sum(int(x["original_S1_photo_count"]) for x in original)!=ORIGINAL_IMAGES or
        len({(x["accepted_key"],normalize(x["S1_original_source_name"])) for x in original})!=ORIGINAL_UNITS):
        raise ValueError("S1 full251/244/170066 photo source/name inventory changed")
    bykey=Counter(x["accepted_key"] for x in original)
    if sum(n==2 for n in bykey.values())!=MIXED_KEYS or any(n>2 for n in bykey.values()):
        raise ValueError("Seven original WCVP collapsed names lost")
    if (len(nuclear)!=MEASURED_NUCLEAR or
        len({x["run"] for x in nuclear})!=MEASURED_NUCLEAR or
        len({x["biosample"] for x in nuclear})!=MEASURED_NUCLEAR or
        len({x["accepted_key"] for x in nuclear})!=MEASURED_NUCLEAR):
        raise ValueError("EAzami preflight 111 real nuclear BioSamples not distinct")
    byname={(x["accepted_key"],normalize(x["taxon_name"])):x for x in nuclear}
    if len(byname)!=MEASURED_NUCLEAR:
        raise ValueError("Duplicated accepted-name/source nuclear specimen")
    out=[]
    seen=set()
    matched=0
    for x in original:
        key=x["accepted_key"]
        taxon=normalize(x["S1_original_source_name"])
        source=byname.get((key,taxon))
        if source: matched+=1;seen.add((key,taxon))
        a=infer(index,taxon)
        out.append({
            "accepted_key":key,
            "original_source_name":taxon,
            "original_S1_photo_count":int(x["original_S1_photo_count"]),
            "nuclear_real120k_sample_in_111":bool(source),
            "source_SRA_run":source["run"] if source else "",
            "source_public_BioSample":source["biosample"] if source else "",
            "source_phylogenetic_tip":source["phylogenetic_tip"] if source else "",
            **a,
            "camera_relative_photo_angle_not_botanical_U_D":True,
            "original_accepted_taxon_taxonomic_revision_verified":False,
        })
    if matched!=MEASURED_NUCLEAR or seen!=set(byname):
        raise ValueError("One or more of 111 source nuclear taxa is not exact original S1 name+key")
    if len(out)!=ORIGINAL_UNITS:
        raise ValueError("251 original S1 source rows lost")
    scored=[x for x in out if x["botanical_U_D_state"] in ("U","D")]
    nuc=[x for x in out if x["nuclear_real120k_sample_in_111"]]
    nucscored=[x for x in nuc if x["botanical_U_D_state"] in ("U","D")]
    conflicts=[x["original_source_name"] for x in out
               if x["botanical_state_status"]=="conflicting_independent_botanical_sources"]
    summary={
      "schema":"azami_V4_T2_original251_111actualSRA_botanical_source_crosswalk_v1",
      "source_original_S1_photos":ORIGINAL_IMAGES,
      "source_original_S1_accepted_keys":ORIGINAL_KEYS,
      "source_original_S1_source_named_units":ORIGINAL_UNITS,
      "source_accepted_name_collisions":MIXED_KEYS,
      "source_real_120k_nuclear_sample_tips":MEASURED_NUCLEAR,
      "original_taxa_with_independently_authored_U_D":len(scored),
      "original_U_D_nodding":sum(x["botanical_U_D_state"]=="D" for x in scored),
      "original_U_D_erect":sum(x["botanical_U_D_state"]=="U" for x in scored),
      "nuclear111_with_independently_authored_U_D":len(nucscored),
      "nuclear111_with_nodding":sum(x["botanical_U_D_state"]=="D" for x in nucscored),
      "nuclear111_with_erect":sum(x["botanical_U_D_state"]=="U" for x in nucscored),
      "nuclear111_without_verified_taxon_concept_state":len(nuc)-len(nucscored),
      "cross_authority_state_disagreements":sorted(set(conflicts)),
      "molecular_BioSample_photophenotype_same_voucher_assessed":False,
      "can_use_camera_image_angles_as_head_nodding":False,
      "U_D_missing_taxa_not_imputed_from_rainfall":True,
      "global111_tip_Mk_with_most_tips_unscored_not_confirmed":True,
      "real_source111_118_ortho_certified_nuclear_tree_completed":False,
      "precipitation_causal_adaptation_inferred":False,
    }
    return out,summary


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source251",type=Path,required=True)
    p.add_argument("--real111-sample-source",type=Path,required=True)
    p.add_argument("--eazami-checkout",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    authorities,blobs=load_authorities(a.eazami_checkout)
    original=read_csv(a.source251)
    source=read_csv(a.real111_sample_source)
    data,summary=join(original,source,authorities)
    summary["EAzami_source_commit"]=EAZAMI_COMMIT
    summary["authority_git_blob_sha1"]=blobs
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"S1_ALL251_ORIGINAL_NAMES_BOTANICAL_UD_SOURCE_SCOPE.csv",data)
    (a.out_dir/"S1_ALL251_REAL111_BOTANICAL_UD_SCOPE_SUMMARY.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("T2_REAL111_ALL251_BOTANICAL_STATE_SCOPE="+json.dumps(summary,sort_keys=True,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()

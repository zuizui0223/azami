"""T3 source-backed photographic eligibility for original 131 Senecioneae rDNA tips.

Pilot only. The user-selected V4 endpoint is S1 continuous EXIF-image-vertical
head axis; a GBIF record with a photo is NOT an angle measurement, and
plant-flora U/D text is neither imported nor substituted.

Conservatively select 15 photo-candidate species from original GenBank ITS
nuclear genealogy by fixed genus/order rule, blind to rainfall, posture and
photographic orientation. Query GBIF species/match and geo-referenced still
image human observations, keep licenses and occurrence links, not image bytes.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from pathlib import Path
from urllib.parse import urlparse
import requests

N_SOURCE=131
GENUS_QUOTAS={"Cremanthodium":3,"Ligularia":3,"Parasenecio":3,"Sinosenecio":3,
              "Senecio":1,"Tussilago":1,"Crassocephalum":1}
API="https://api.gbif.org/v1"
USER_AGENT="azami-v4-S1-continuous-T3-photo-source-feasibility/1.0 scientific source audit"

def roster(path):
    raw=Path(path).read_bytes()
    rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    if len(rows)!=N_SOURCE or len({x["published_source_species"] for x in rows})!=N_SOURCE:
        raise ValueError("131 original ITS source concepts must be present")
    if any(not x["published_ITS_GenBank_accession"] for x in rows):
        raise ValueError("Cannot use taxon without authentic source ITS")
    return rows,hashlib.sha256(raw).hexdigest()

def sample_source_taxa(rows):
    selected=[]
    for genus,quota in GENUS_QUOTAS.items():
        found=sorted(x["published_source_species"] for x in rows
                     if x["published_source_species"].startswith(genus+" "))
        if len(found)<quota:
            raise ValueError("Not enough original taxonomy: "+genus)
        selected+=found[:quota]
    if len(selected)!=15 or len(set(selected))!=15:
        raise ValueError("Original 15 name pilot selection drift")
    return selected

def canonical(name):
    return " ".join(name.split()[:2])

def read_json(session,path,params,max_attempts=3):
    err=None
    for attempt in range(max_attempts):
        try:
            resp=session.get(API+path,params=params,timeout=35)
            if resp.status_code==429 or resp.status_code>=500:
                err=f"HTTP_{resp.status_code}"
                time.sleep(min(2.0*(attempt+1),6))
                continue
            resp.raise_for_status()
            return resp.json()
        except (requests.RequestException,ValueError) as e:
            err=str(e)
            time.sleep(attempt+1)
    raise RuntimeError("GBIF_HTTP_SOURCE_NOT_VERIFIED: "+str(err))

def strict_match(resp,name):
    expected=canonical(name)
    got=canonical(str(resp.get("scientificName","")))
    if not(resp.get("rank")=="SPECIES"
           and resp.get("matchType")=="EXACT"
           and int(resp.get("confidence",0))>=90
           and got==expected):
        return None
    return int(resp.get("usageKey",0) or 0) or None

def valid_media(media):
    return any(str(m.get("type","")).lower()=="stillimage"
               and str(m.get("identifier","")).startswith(("https://","http://"))
               for m in media or [])

def trial(session,names):
    details=[];photos=[]
    for name in names:
        match=read_json(session,"/species/match",{"name":name,"rank":"SPECIES","strict":"true"})
        gbif_id=strict_match(match,name)
        status="NO_EXACT_GBIF_SPECIES_MATCH" if not gbif_id else "GBIF_QUERY_COMPLETED"
        n=0
        if gbif_id:
            # Limit is retrieval ceiling, NOT a global count of all valid photos.
            data=read_json(session,"/occurrence/search",{
                "taxonKey":str(gbif_id),"mediaType":"StillImage",
                "hasCoordinate":"true","basisOfRecord":"HUMAN_OBSERVATION",
                "limit":30,"offset":0})
            for r in data.get("results",[]):
                species=canonical(str(r.get("species","")))
                if species!=canonical(name) or not valid_media(r.get("media",[])):
                    continue
                accepted_media=next(m for m in r["media"]
                    if str(m.get("type","")).lower()=="stillimage"
                    and str(m.get("identifier","")).startswith(("https://","http://")))
                photos.append({
                    "source_taxon":name,"gbif_taxon_key":gbif_id,
                    "gbif_occurrence_key":r.get("key"),
                    "GBIF_occurrence_record_url":"https://www.gbif.org/occurrence/"+str(r.get("key")),
                    "source_image_url":accepted_media.get("identifier"),
                    "license":r.get("license") or accepted_media.get("license") or "",
                    "rightsHolder":r.get("rightsHolder") or accepted_media.get("rightsHolder") or "",
                    "observation_year":r.get("year"),
                    "not_an_S1_continuous_head_angle_measurement":True
                })
                n+=1
        details.append({"source_taxon":name,"strict_original_genbank_ITS_tip":True,
                        "GBIF_exact_species_key":gbif_id,"status":status,
                        "retrieved_usable_stillimage_occurrences_first30":n,
                        "GBIF_media_retrieval_ceiling_per_taxon":30})
        time.sleep(.25)
    return details,photos

def run(source,outdir):
    rows,sha=roster(source)
    chosen=sample_source_taxa(rows)
    with requests.Session() as session:
        session.headers.update({"User-Agent":USER_AGENT})
        details,images=trial(session,chosen)
    root=Path(outdir);root.mkdir(parents=True,exist_ok=True)
    output={
        "schema":"T3_S1_photo_endpoint_external_131_GenBank_true_ITS_blind15_source_pilot_v1",
        "status":"GBIF_REAL_PUBLIC_PHOTO_SOURCE_METADATA_AUDIT_NOT_ANGLE_MEASUREMENT",
        "original_131_GenBank_tip_roster_sha256":sha,
        "source_ITS_taxa":N_SOURCE,"pilot_taxa":len(chosen),
        "selection":"fixed_3_first_source_names_each_major_genus_plus_3_singleton_genus__blind_to_climate_or_U_D",
        "selection_by_flora_posture":False,"original_photo_angular_measurements":0,
        "original_photo_files_downloaded":0,
        "n_strict_GBIF_taxon_matches":sum(x["GBIF_exact_species_key"] is not None for x in details),
        "n_photo_observations_in_first30_per_taxon":len(images),
        "n_species_with_source_photos":sum(x["retrieved_usable_stillimage_occurrences_first30"]>0 for x in details),
        "not_photos_unavailable_beyond_first30":True,
        "not_automatically_transferrable_Cirsium_detector":True,
        "source_licenses_and_attribution_retained_in_csv":True,
        "no_reframing_as_true_gravity_or_botanical_downward":True,
        "has_nuclear_ITS_tip_not_genus_wide_nuclear_species_tree":True,
        "T3_same_S1_angle_comparative_statistical_status":"BLOCKED_UNTIL_VALIDATED_REAL_IMAGE_MEASUREMENTS",
        "taxa":details
    }
    (root/"T3_PUBLIC_IMAGE_SOURCE_ELIGIBILITY.json").write_text(
        json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    if images:
        with (root/"T3_PUBLIC_IMAGES_LICENSED_SOURCE_CANDIDATES.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(images[0]))
            w.writeheader();w.writerows(images)
    print("T3_S1_ENDPOINT_IMAGE_SOURCE_PILOT="+json.dumps({
        "taxa":len(chosen),"matched":output["n_strict_GBIF_taxon_matches"],
        "with_images":output["n_species_with_source_photos"],
        "source_photo_records":len(images),"measured_angle_count":0},sort_keys=True))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    run(a.source,a.out_dir)
if __name__=="__main__":main()

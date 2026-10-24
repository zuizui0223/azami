"""T3 real independent-Senecioneae photo/FOUR-climate S1→S2→T2→T3 bridge.

Census ALL 131 independent originally sequence-linked GenBank ITS taxa, not
posture-selected exemplars, with first-30 geolocated HUMAN_OBSERVATION
StillImages per exact GBIF binomial. Sample ACTUAL original geographic sites
against exactly the S1/S2 full-source CHELSA v2.1 BIO18/BIO05, original
CHELSA shortwave RSDS, and WorldClim v2.1 30 arcsec elevation.

NO T3 photographic orientation measurement or cross-tribe evolutionary
effect is claimed. Original Cirsium-specific YOLO classifier is NOT applied
without validated external-domain detector/keypoint quality. Linked
multicopy ITS/ETS is NOT an independent multilocus species phylogeny.
Never export GBIF individual lat/lon, observer IDs, or downloaded image bytes.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from analysis.v4 import t3_s1_continuous_external_photo_source_pilot as original
from analysis.v4 import s1_s2_real_bio18_bio05_rsds_elevation as s1s2
from analysis.v4.native_cached_climate import sample_cached_pixels,verified_cache

SOURCE_BASE="https://api.gbif.org/v1"
FOUR=s1s2.EXPOSURES
MAX_PHOTO_RECORDS=30
MIN_INDEPENDENT_PHOTO_SUPPORT=10
MIN_INDEPENDENT_CELLS=3
HEADERS={"User-Agent":"azami-S1S2T2T3-independent-Senecioneae-photos-4-climate/1.0"}
ALLOWED={"http://creativecommons.org/licenses/by/4.0/legalcode",
         "http://creativecommons.org/licenses/by-nc/4.0/legalcode",
         "http://creativecommons.org/publicdomain/zero/1.0/"}

def read_json(session,path,params,max_attempts=4):
    err=None
    for attempt in range(max_attempts):
        try:
            response=session.get(SOURCE_BASE+path,params=params,timeout=(8,35))
            if response.status_code==429 or response.status_code>=500:
                err=f"HTTP_{response.status_code}"
                time.sleep(1+2*attempt)
                continue
            response.raise_for_status()
            return response.json()
        except (requests.RequestException,ValueError) as e:
            err=type(e).__name__
            time.sleep(attempt+1)
    raise RuntimeError("GBIF_SOURCE_UNAVAILABLE_"+str(err))

def photo_entry(observation,source_name,taxonkey):
    name=original.canonical(str(observation.get("species","")))
    if name != original.canonical(source_name):return None
    lat=observation.get("decimalLatitude")
    lon=observation.get("decimalLongitude")
    try:
        lat=float(lat);lon=float(lon)
        if not np.isfinite(lat) or not np.isfinite(lon) or abs(lat)>90 or abs(lon)>180:
            return None
    except (TypeError,ValueError):
        return None
    for media in observation.get("media",[]):
        link=str(media.get("identifier",""))
        if str(media.get("type","")).casefold()!="stillimage" or not link.startswith("https://"):
            continue
        license=observation.get("license") or media.get("license") or ""
        if license not in ALLOWED:
            continue
        return {
            "source_taxon":source_name,"gbif_taxon_key":taxonkey,
            "gbif_occurrence_key":str(observation["key"]),
            "source_image_url":link,
            "image_license_url":license,
            "source_rightsHolder":str(observation.get("rightsHolder") or media.get("rightsHolder") or ""),
            "_latitude":lat,"_longitude":lon,
        }
    return None

def collect(rows,session):
    # Full genuine biological source ITS taxa, no U/D or climate selection.
    names=sorted(r["published_source_species"] for r in rows)
    if len(names)!=original.N_SOURCE or len(set(names))!=original.N_SOURCE:
        raise ValueError("Published 131 original DNA species census changed")
    candidates=[];details=[]
    for name in names:
        matched=read_json(session,"/species/match",{
            "name":name,"rank":"SPECIES","strict":"true"})
        taxonkey=original.strict_match(matched,name)
        info={"original_GenBank_ITS_taxon":name,
              "original_GenBank_ITS_exists":True,
              "GBIF_strict_species_key":taxonkey,
              "retrieval_ceiling_per_original_taxon":MAX_PHOTO_RECORDS,
              "n_photo_observation_valid":0,
              "status":"NO_EXACT_GBIF_SPECIES_MATCH" if not taxonkey else "GBIF_SOURCE_QUERIED"}
        if taxonkey:
            got=read_json(session,"/occurrence/search",{
                "taxonKey":taxonkey,
                "mediaType":"StillImage",
                "hasCoordinate":"true",
                "basisOfRecord":"HUMAN_OBSERVATION",
                "limit":MAX_PHOTO_RECORDS,"offset":0})
            seen=set()
            for row in got.get("results",[]):
                item=photo_entry(row,name,taxonkey)
                if item is None or item["gbif_occurrence_key"] in seen:
                    continue
                seen.add(item["gbif_occurrence_key"])
                candidates.append(item)
            info["n_photo_observation_valid"]=len(seen)
            info["GBIF_query_results_total_reported"]=int(got.get("count",0))
            info["more_results_available_beyond_first30"]=bool(got.get("count",0)>MAX_PHOTO_RECORDS)
        details.append(info)
        time.sleep(.08)
    return details,candidates

def exact_four_climates(frame,source_root,rsds_root,meta_path):
    if frame.empty:raise ValueError("No independent source photographs at all")
    lat=frame["_latitude"].to_numpy(float)
    lon=frame["_longitude"].to_numpy(float)
    if len(set(frame["gbif_occurrence_key"]))!=len(frame):
        # Occurrence keys should not represent independent measurements for
        # more than one source biological taxon.
        raise ValueError("Same GBIF observation leaked to multiple source taxa")
    result=frame.copy()
    evidence={}
    for field,member in (("chelsa_bio18","CHELSA_bio18_1981-2010_V.2.1.tif"),
                         ("chelsa_bio05","CHELSA_bio05_1981-2010_V.2.1.tif"),
                         ("elevation_m","wc2.1_30s_elev.tif")):
        data,audit=s1s2.raster_sample(Path(source_root)/member,lat,lon)
        result[field]=data
        evidence[field]=audit
    metadata=json.loads(Path(meta_path).read_text())["chelsa_rsds_mean"]
    root=Path(rsds_root)/"chelsa_rsds_mean"
    grid=verified_cache(root,metadata)
    raw=sample_cached_pixels(root,grid,lat,lon)
    result["chelsa_rsds_mean"]=raw*grid["scale"]+grid["offset"]
    evidence["chelsa_rsds_mean"]={"full_original_raster_sha256":grid["source_sha256"],
      "source":"S1 S2 original SHA-checked and nonresampled 316 tile cache",
      "n_site_complete":int(np.isfinite(result["chelsa_rsds_mean"]).sum())}
    for field in FOUR:
        if np.isinf(result[field]).any():
            raise ValueError("Nonfinite environmental pixel exposure in "+field)
    return result,evidence

def aggregate(rows,census,sha,source_meta):
    table=pd.DataFrame(rows)
    climate=table.dropna(subset=list(FOUR)).copy()
    if climate.empty:raise ValueError("Zero four-environment photo source support")
    climate["_cell"]=(np.floor((climate["_latitude"]+90)/5).astype(int)*72+
                      np.floor((climate["_longitude"]+180)/5).astype(int))
    summary=[]
    grouped=climate.groupby("source_taxon")
    for info in census:
        taxon=info["original_GenBank_ITS_taxon"]
        sample=grouped.get_group(taxon) if taxon in grouped.groups else climate.iloc[:0]
        count=len(sample)
        cells=int(sample["_cell"].nunique()) if count else 0
        med={field:float(sample[field].median()) if count else None for field in FOUR}
        rep={field:[float(v) for v in sample[field].quantile([.05,.95]).to_numpy()] if count else None for field in FOUR}
        summary.append({
            **info,
            "n_four_environment_complete_geo_photo_observations":count,
            "n_independent_5deg_geographic_cells":cells,
            "niche_BIO18_BIO05_RSDS_elev_medians":med,
            "niche_5th_95th_within_first30_observed_photo_sites":rep,
            "candidate_for_manual_S1_angle_protocol_audit":
                bool(count>=MIN_INDEPENDENT_PHOTO_SUPPORT and cells>=MIN_INDEPENDENT_CELLS),
            "has_genuinely_measured_S1_equivalent_photo_angle":False,
        })
    genera={}
    for v in summary:
        genus=v["original_GenBank_ITS_taxon"].split()[0]
        x=genera.setdefault(genus,{"n_original_ITS_source_taxa":0,
                                    "n_strict_GBIF_matched":0,
                                    "n_taxa_with_images":0,
                                    "n_fourclimate_photo_records":0,
                                    "n_taxa_geographically_replicated_photo_pilot_candidate":0})
        x["n_original_ITS_source_taxa"]+=1
        x["n_strict_GBIF_matched"]+=int(v["GBIF_strict_species_key"] is not None)
        x["n_taxa_with_images"]+=int(v["n_photo_observation_valid"]>0)
        x["n_fourclimate_photo_records"]+=v["n_four_environment_complete_geo_photo_observations"]
        x["n_taxa_geographically_replicated_photo_pilot_candidate"]+=int(
            v["candidate_for_manual_S1_angle_protocol_audit"])
    return {
        "schema":"T3_other_tribe_131_authentic_nuclear_ITS_photo_fourclimate_source_census_20261024_v1",
        "status":"REAL_OTHER_TRIBE_GBIF_FOUR_ENVIRONMENT_PHOTO_SOURCE_QC__HEAD_ANGLE_STILL_NOT_MEASURED",
        "source_131_original_GenBank_ITS_roster_sha256":sha,
        "original_genbank_ITS_biological_taxa":131,
        "four_axes_exactly_S1_S2":list(FOUR),
        "original_131_taxon_selection_not_based_on_photography_or_environment_or_U_D":True,
        "max_first30_per_species":"availability_screen_only_NOT_full_gbif_occurrence_census",
        "n_GBIF_exact_taxon_matches":sum(x["GBIF_strict_species_key"] is not None for x in summary),
        "n_source_species_with_photos":sum(x["n_photo_observation_valid"]>0 for x in summary),
        "n_georeferenced_photo_records_before_all4_climate":len(table),
        "n_records_4_environment_complete":len(climate),
        "n_species_4_environment_complete":sum(x["n_four_environment_complete_geo_photo_observations"]>0 for x in summary),
        "n_species_qualify_for_MANUAL_photo_angle_audit":sum(
            x["candidate_for_manual_S1_angle_protocol_audit"] for x in summary),
        "n_genera_with_photo_records":sum(x["n_taxa_with_images"]>0 for x in genera.values()),
        "n_actual_qualified_T3_S1_equivalent_continuous_head_angle_measurements":0,
        "no_cross_tribe_angle_environment_effect_fitted":True,
        "original_ITS_ETS_one_linked_rDNA_gene_tree_not_independent_multilocus_species_tree":True,
        "matched_photo_images_cannot_be_treated_as_S1_angle_or_gravity":True,
        "group_specific_head_and_pedicel_detection_quality_not_yet_validated":True,
        "sampling_observation_lat_lon_and_image_bytes_never_exported":True,
        "validity":"SOURCE_AND_NICHE_FEASIBILITY_ONLY",
        "climate_provenance":source_meta,"by_genus":genera,"by_species":summary,
    }

def execute(ITS_roster,source_root,rsds_cache_root,rsds_metadata,outdir):
    roster,sha=original.roster(ITS_roster)
    with requests.Session() as session:
        session.headers.update(HEADERS)
        details,rows=collect(roster,session)
    if len(details)!=original.N_SOURCE:
        raise ValueError("Real published source 131 taxon census incomplete")
    if len(rows)<50:
        raise ValueError("Too few actual georeferenced other-tribe photo candidate sources")
    frame,evidence=exact_four_climates(pd.DataFrame(rows),source_root,rsds_cache_root,rsds_metadata)
    output=aggregate(frame.to_dict("records"),details,sha,evidence)
    root=Path(outdir);root.mkdir(parents=True,exist_ok=True)
    (root/"T3_OTHER_TRIBE_131_TRUE_DNA_TAXA_PUBLIC_IMAGES_FOUR_CLIMATES.json").write_text(
        json.dumps(output,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    # Public original photo URLs are only references. No locations/observer IDs,
    # raw photos, original image frames, or invented head angles are released.
    public=frame.drop(columns=["_latitude","_longitude",*FOUR]).copy()
    public["not_a_measured_S1_head_angle"]=True
    public.to_csv(root/"T3_SOURCE_PUBLIC_PHOTO_REFERENCE_LICENSES.csv",index=False)
    print("T3_TRUE_SOURCE_OTHER_TRIBE_PHOTO_FOUR_CLIMATE="+json.dumps({
        "source_ITS_taxa":131,
        "exact_GBIF_taxa":output["n_GBIF_exact_taxon_matches"],
        "photo_taxa":output["n_source_species_with_photos"],
        "photo_records":output["n_georeferenced_photo_records_before_all4_climate"],
        "four_climate_photo_records":output["n_records_4_environment_complete"],
        "qualified_for_manual_photo_angle_review":output["n_species_qualify_for_MANUAL_photo_angle_audit"],
        "n_measured_angles":0},sort_keys=True),flush=True)
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ITS-roster",type=Path,required=True)
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--rsds-cache-root",type=Path,required=True)
    p.add_argument("--rsds-metadata",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    x=p.parse_args()
    execute(x.ITS_roster,x.source_root,x.rsds_cache_root,x.rsds_metadata,x.outdir)
if __name__=="__main__":main()

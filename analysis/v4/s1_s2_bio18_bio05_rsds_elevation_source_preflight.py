"""Preflight actual CHELSA BIO5/BIO18 and elevation before S1/S2 trait inspection.

A 206 partial-response with an exact Content-Range and TIFF header is required
to claim range sampling ready. Never silently download an entire global TIFF
on an ignored Range header. The original five S1/S2 axes remain immutable.
This is environmental-source feasibility, NOT biological inference.
"""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import requests

HEADERS={"Range":"bytes=0-16383","Accept-Encoding":"identity",
         "User-Agent":"azami-V4-science-CHELSA-public-range-preflight/1.0"}
SOURCE={
 "chelsa_bio05":[
  "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio05/1981-2010/CHELSA_bio05_1981-2010_V.2.1.tif",
  "https://os.zhdk.cloud.switch.ch/envicloud/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio5_1981-2010_V.2.1.tif",
  "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio5_1981-2010_V.2.1.tif"],
 "chelsa_bio18":[
  "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio18/1981-2010/CHELSA_bio18_1981-2010_V.2.1.tif",
  "https://os.zhdk.cloud.switch.ch/envicloud/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio18_1981-2010_V.2.1.tif",
  "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio18_1981-2010_V.2.1.tif"],
 "elevation_worldclim":[
  "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip",
  "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_elev.zip"],
}
def get_range(session,url):
    begin=time.monotonic()
    try:
        with session.get(url,headers=HEADERS,timeout=(6,14),
                         stream=True,allow_redirects=True) as r:
            status=r.status_code
            ctype=r.headers.get("Content-Type","")
            cr=r.headers.get("Content-Range","")
            cl=r.headers.get("Content-Length")
            evidence=r.raw.read(32) if status in (200,206) else b""
            tiff=evidence[:4] in (b"II*\0",b"MM\0*",b"II+\0",b"MM\0+")
            zip_file=evidence.startswith(b"PK\x03\x04")
            partial=status==206 and cr.startswith("bytes 0-")
            return {"url":url,"status_code":status,"content_type":ctype,
                "content_range":cr,"content_length":cl,
                "range_206_verified":partial,
                "first_bytes_TIFF_header":tiff,
                "first_eight_bytes_hex":evidence[:8].hex(),
                "first_bytes_ZIP_header":zip_file,
                "seconds":round(time.monotonic()-begin,2),
                "ready_to_consider_GDAL_remote_range":
                    bool(partial and (tiff or zip_file))}
    except requests.RequestException as exc:
        return {"url":url,"status":"REMOTE_UNAVAILABLE",
                "reason":type(exc).__name__,"seconds":round(time.monotonic()-begin,2),
                "ready_to_consider_GDAL_remote_range":False}
def main():
    a=argparse.ArgumentParser();a.add_argument("--out",required=True,type=Path);args=a.parse_args()
    out={"schema":"azami_V4_S1_S2_BIO18_BIO5_RSDS_ELEVATION_SOURCE_PREFLIGHT_V1",
      "climate_terms":{"BIO18":"Precipitation of Warmest Quarter (CHELSA v2.1, 1981-2010)",
        "BIO5":"Mean Daily Maximum Temperature of Warmest Month (CHELSA v2.1, 1981-2010)",
        "RSDS":"original S1/S2 CHELSA v2.1 cached radiation not newly acquired",
        "elevation":"WorldClim 2.1 elevation, source distinct from CHELSA"},
      "not_biological_result":True,
      "no_S1_S2_private_photo_locations_fetched":True,
      "not_old_BIO1_BIO15_relabeled":True,"sources":{}}
    with requests.Session() as s:
        for variable,urls in SOURCE.items():
            candidates=[get_range(s,url) for url in urls]
            out["sources"][variable]={"candidates":candidates,
                "first_available_206_url":next((r["url"] for r in candidates
                    if r["ready_to_consider_GDAL_remote_range"]),None)}
    out["raster_provenance_status"]=("BOUNDED_BYTE_RANGE_AVAILABLE_FOR_ALL_NEW_AXES"
       if all(x["first_available_206_url"] for x in out["sources"].values())
       else "BLOCKED_ONE_OR_MORE_NEW_RASTER_SOURCES_NO_VERIFIED_BYTE_RANGE")
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("NEW_S1_S2_4CLIMATE_SOURCE_RANGE="+json.dumps({
      k:v["first_available_206_url"] for k,v in out["sources"].items()},sort_keys=True))
    print("SOURCE_STATUS="+out["raster_provenance_status"])
if __name__=="__main__":main()

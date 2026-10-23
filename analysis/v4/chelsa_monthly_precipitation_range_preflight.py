"""Bounded live provenance probe of CHELSA v2.1 monthly precipitation normals.

A range-aware, authenticated metadata+byte receipt is required before
remote GDAL raster sampling. Never request all 12 multi-GB global TIFFs in
GitHub Actions without a bounded, source-validated plan. HTTP 200 on Range
is NOT real byte-range support and should NOT trigger full TIFF downloads.

Provider-supported 1981-2010 monthly precipitation is water-equivalent
(liquid + snow), NOT rainy-day frequency or a rain-on-anthesis measurement.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import requests

BASES=(
    "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/pr/",
    "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/pr/",
    "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/climatologies/1981-2010/pr/",
)
MONTHS=(1,6,12)
MAX_PROBES=len(BASES)*len(MONTHS)
RANGE_N=2048


def file_url(base,month):
    if base not in BASES or not 1<=month<=12:
        raise ValueError("Unverified CHELSA source base or month")
    return base+f"CHELSA_pr_{month:02d}_1981-2010_V.2.1.tif"


def inspect_source(session=None,months=MONTHS):
    if tuple(months)!=MONTHS:
        raise ValueError("Stable monthly-raster preflight months required")
    s=session or requests.Session()
    evidence=[]
    for base in BASES:
        for month in MONTHS:
            url=file_url(base,month)
            entry={"url":url,"month":month,"range_requested":f"bytes=0-{RANGE_N-1}",
                   "valid_raster_source_with_true_http_206":False}
            try:
                response=s.get(
                    url,headers={"Range":f"bytes=0-{RANGE_N-1}",
                                 "Accept-Encoding":"identity"},
                    timeout=(8,18),stream=True,allow_redirects=True)
                try:
                    entry["http_status"]=int(response.status_code)
                    entry["content_range"]=response.headers.get("Content-Range","")
                    entry["content_type"]=response.headers.get("Content-Type","")
                    entry["archive_total_bytes"]=None
                    match=re.fullmatch(rf"bytes 0-{RANGE_N-1}/(\d+)",
                                       entry["content_range"])
                    if response.status_code==206 and match:
                        entry["archive_total_bytes"]=int(match.group(1))
                        entry["valid_raster_source_with_true_http_206"]=True
                        # TIFF little/big endian magic only, never read remainder.
                        sample=b""
                        for ch in response.iter_content(262144):
                            sample+=ch
                            if len(sample)>=4:break
                        entry["valid_tiff_magic"]=sample[:4] in (
                            b"II*\x00",b"MM\x00*",b"II+\x00",b"MM\x00+")
                        entry["valid_raster_source_with_true_http_206"]=(
                            entry["valid_raster_source_with_true_http_206"]
                            and entry["valid_tiff_magic"])
                finally:response.close()
            except requests.RequestException as e:
                entry["status"]="REMOTE_CHELSA_UNAVAILABLE_NOT_MONTHLY_PRECIP_ABSENCE"
                entry["transport_error_class"]=type(e).__name__
            evidence.append(entry)
    return {
        "status":("AUTHENTIC_CHELSA_MONTHLY_RANGE_SOURCE_READY" if
                  any(x["valid_raster_source_with_true_http_206"] for x in evidence)
                  else "CHELSA_MONTHLY_REMOTE_RANGE_NOT_YET_VERIFIED"),
        "source":"CHELSA V2.1 climatologies / monthly pr / 1981-2010",
        "provider":"https://www.chelsa-climate.org/datasets/chelsa_climatologies",
        "n_pinned_probe_months":len(MONTHS),
        "n_tested_endpoints":len(evidence),
        "source_1981_2010_water_equivalent_includes_snow":True,
        "site_actual_daily_precipitation_observed":False,
        "rainy_day_frequency_from_monthly_normal_proven":False,
        "head_photo_orientation_or_bio12_used_to_select_source":False,
        "probes":evidence,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result=inspect_source()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,ensure_ascii=False,
                                indent=2)+"\n",encoding="utf-8")
    print("CHELSA_MONTHLY_RAIN_RANGES="+json.dumps(result,ensure_ascii=False))
    if result["status"]!="AUTHENTIC_CHELSA_MONTHLY_RANGE_SOURCE_READY":
        print("NOT_READY: no 12-layer cloud access; do not derive rainfall-frequency claims")


if __name__=="__main__":
    main()

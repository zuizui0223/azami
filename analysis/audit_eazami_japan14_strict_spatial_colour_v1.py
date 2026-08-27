#!/usr/bin/env python3
"""Audit all exact Japan38 colour concepts inside the frozen Azami strict-spatial cohort.

This is a coverage/provenance audit, not a trait-history test. It asks whether
already-generated detector-crop colour measurements can supply additional
Japan-local concepts before any new media are acquired.

The admission rule is frozen before execution:
- exact concept/taxon identity from eazami_japan14_exact_concepts_v1.csv;
- finite CIELAB L* with >=1 usable detected head in the strict-spatial row;
- coordinate lies inside the same operational Japan window used by EAzami v1;
- >=2 Japan-local observations AND >=2 distinct frozen ``analysis_cell`` values.

The 2-cell rule was predeclared in EAzami readiness v14 before this audit was
executed. It does not replace or delete earlier n>=5 high-depth analyses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as h:
        return list(csv.DictReader(h))


def sha256(path: Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def finite(value):
    try: x=float(value)
    except (TypeError,ValueError): return None
    return x if math.isfinite(x) else None


def in_japan_window(lat, lon):
    return (
        (24.0 <= lat < 30.8 and 122.0 <= lon <= 131.5)
        or (30.8 <= lat < 41.5 and 129.0 <= lon <= 142.5)
        or (41.5 <= lat <= 45.7 and 139.0 <= lon <= 145.9)
    )


def write_csv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:
        path.write_text("",encoding="utf-8"); return
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cohort",type=Path,required=True)
    p.add_argument("--concepts",type=Path,required=True)
    p.add_argument("--summary-output",type=Path,required=True)
    p.add_argument("--observations-output",type=Path,required=True)
    args=p.parse_args()

    concepts=read_csv(args.concepts)
    if len(concepts)!=14 or len({r['paper_japan_member_id'] for r in concepts})!=14:
        raise ValueError("expected exactly 14 unique exact concepts")
    cohort=read_csv(args.cohort)
    if len(cohort)!=46276:
        raise ValueError(f"frozen strict-spatial row count changed: {len(cohort)}")
    required={'obs_id','taxon_name','latitude','longitude','observed_on','analysis_cell','corolla_lab_lightness_median','corolla_lab_lightness_n_usable_heads'}
    if not required.issubset(cohort[0]):
        raise ValueError(f"cohort columns missing {sorted(required-set(cohort[0]))}")

    concept_rows=[]; admitted=[]; obs_out=[]
    for c in concepts:
        mid,taxon=c['paper_japan_member_id'],c['taxon_name']
        tax_rows=[r for r in cohort if r['taxon_name']==taxon]
        colour=[]; japan=[]
        for r in tax_rows:
            l=finite(r.get('corolla_lab_lightness_median'))
            n=finite(r.get('corolla_lab_lightness_n_usable_heads')) or 0
            lat=finite(r.get('latitude')); lon=finite(r.get('longitude'))
            if l is None or n < 1: continue
            colour.append(r)
            if lat is None or lon is None or not in_japan_window(lat,lon): continue
            rr={
                'paper_japan_member_id':mid,'taxon_name':taxon,'obs_id':r['obs_id'],
                'latitude':lat,'longitude':lon,'observed_on':r.get('observed_on',''),
                'analysis_cell':r.get('analysis_cell',''),'corolla_lab_lightness_median':l,
                'corolla_lab_lightness_n_usable_heads':int(n),
            }
            japan.append(rr); obs_out.append(rr)
        cells=sorted({r['analysis_cell'] for r in japan if r['analysis_cell']})
        blocks=len({(round(float(r['latitude'])),round(float(r['longitude']))) for r in japan})
        passed=len(japan)>=2 and len(cells)>=2
        if passed: admitted.append(mid)
        values=[float(r['corolla_lab_lightness_median']) for r in japan]
        concept_rows.append({
            'paper_japan_member_id':mid,'taxon_name':taxon,
            'strict_spatial_rows_all':len(tax_rows),'colour_usable_rows_global':len(colour),
            'japan_colour_usable_rows':len(japan),'japan_distinct_analysis_cells':len(cells),
            'japan_approx_1deg_blocks':blocks,
            'japan_lightness_median':statistics.median(values) if values else None,
            'admission_rule_pass':passed,
        })

    result={
        'contract_version':'eazami_japan14_strict_spatial_colour_audit_v1',
        'status_date':'2026-08-27',
        'source':{
            'repo':'zuizui0223/azami','workflow_run':29306454759,
            'artifact_id':8301295025,'artifact_name':'ch1-strict-spatial-cohort-only-20260713',
            'artifact_digest':'sha256:f453e36b1867bbb0a239253c2cdfa4a4a6b8f6117a9e814db7224e8b060a440a',
            'source_file':'strict_spatial_thinned_observations.csv',
            'source_file_sha256':sha256(args.cohort),'source_rows':len(cohort),
        },
        'operational_japan_window':'same three frozen coordinate boxes as EAzami geographic-provenance gate v1',
        'admission_rule':'>=2 Japan-local colour-usable strict-spatial observations AND >=2 distinct frozen analysis_cell values; exact concept identity required',
        'admission_threshold_frozen_before_audit':True,
        'concepts':concept_rows,
        'admitted_concepts':admitted,
        'admitted_concept_count':len(admitted),
        'claim_boundary':'Existing-data coverage audit only. Passing admits a continuous local L* proxy to a later predeclared sensitivity panel; it does not define fixed colour, W/C state, transition direction, ancestry, convergence, adaptation or rate.'
    }
    args.summary_output.parent.mkdir(parents=True,exist_ok=True)
    args.summary_output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    write_csv(args.observations_output,obs_out)
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=='__main__': main()

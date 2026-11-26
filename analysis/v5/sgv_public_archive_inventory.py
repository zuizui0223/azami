"""Inventory public ENA read-run metadata, not genetic selection evidence."""
import io
import json
import re
from pathlib import Path
import pandas as pd
import requests

OUT=Path("work/public_sgv_inventory")
API="https://www.ebi.ac.uk/ena/portal/api"
TAXA={"Cirsium arvense":41550,"Cirsium vulgare":92907}
REQUIRED=["run_accession","study_accession","sample_accession","scientific_name","library_strategy","library_source"]
WANTED=REQUIRED+["secondary_study_accession","secondary_sample_accession","library_selection",
    "collection_date","country","geo_loc_name","tissue_type","sample_title","read_count","base_count"]
# Do not export raw free-text sample titles, locations, identifiers of collectors,
# or coordinates. Accessions are public archival identifiers.
PUBLIC=[c for c in WANTED if c!="sample_title"]

def get_tsv(session,endpoint,params):
    r=session.get(API+"/"+endpoint,params=params,timeout=(20,120))
    r.raise_for_status()
    if not r.text.strip():raise RuntimeError("Empty ENA response")
    return pd.read_csv(io.StringIO(r.text),sep="\t",dtype=str,keep_default_na=False)

def pick_fields(columns):
    missing=set(REQUIRED)-set(columns)
    if missing:raise ValueError("Required ENA fields unavailable: "+",".join(sorted(missing)))
    return [k for k in WANTED if k in columns]

def normalize(d,sp):
    d=d.copy()
    assert set(REQUIRED)<=set(d)
    if d.run_accession.eq("").any() or d.run_accession.duplicated().any():
        raise ValueError("Missing or duplicate run accession")
    # Exact binomial is the requested taxon; report exclusion rather than silently pool synonyms/hybrids.
    d["focal_taxon"]=sp
    d["exact_name_match"]=d.scientific_name.str.strip().eq(sp)
    for c in ["secondary_study_accession","secondary_sample_accession","collection_date","country","geo_loc_name","tissue_type","library_selection"]:
        if c not in d:d[c]=""
    d["reported_country"]=d.country.where(d.country.ne(""),d.geo_loc_name).str.split(":").str[0].str.strip()
    d["reported_year"]=d.collection_date.str.extract(r"^((?:18|19|20)\d{2})(?:$|[-/])",expand=False).fillna("")
    d["sample_key"]=d.sample_accession.where(d.sample_accession.ne(""),d.secondary_sample_accession)
    d["project_key"]=d.secondary_study_accession.where(d.secondary_study_accession.ne(""),d.study_accession)
    return d

def summarize(d):
    f=d[d.exact_name_match]
    rows=[]
    for (sp,project,strategy,source),g in f.groupby(["focal_taxon","project_key","library_strategy","library_source"],dropna=False):
        s=g[g.sample_key.ne("")].drop_duplicates("sample_key")
        known=s[s.reported_country.ne("")]
        rows.append(dict(species=sp,project=project,library_strategy=strategy,library_source=source,
            runs=len(g),unique_samples=len(s),runs_without_sample=int(g.sample_key.eq("").sum()),
            samples_with_country=len(known),countries=";".join(sorted(known.reported_country.unique())),
            samples_with_year=int(s.reported_year.ne("").sum()),
            possible_nuclear_population_screen=(source.upper()=="GENOMIC" and strategy.upper() in {"RAD-SEQ","WGS","GBS","TARGETED-CAPTURE","AMPLICON"} and len(s)>=20 and known.reported_country.nunique()>=2)))
    return pd.DataFrame(rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers["User-Agent"]="azami-public-archive-feasibility-audit/1.0"
    fields=get_tsv(session,"returnFields",{"result":"read_run"})
    if "columnId" not in fields:raise ValueError("Unexpected ENA field catalog")
    chosen=pick_fields(fields.columnId)
    data=[];meta=[]
    for sp,taxid in TAXA.items():
        d=get_tsv(session,"search",{"result":"read_run","query":f"tax_eq({taxid})","fields":",".join(chosen),"limit":0,"format":"tsv"})
        q=normalize(d,sp)
        data.append(q)
        meta.append(dict(species=sp,taxid=taxid,runs_returned=len(q),exact_name_runs=int(q.exact_name_match.sum()),
                         distinct_samples=q.loc[q.exact_name_match & q.sample_key.ne(""),"sample_key"].nunique()))
    full=pd.concat(data,ignore_index=True)
    summary=summarize(full)
    summary.to_csv(OUT/"project_inventory.csv",index=False)
    safe=[c for c in PUBLIC if c in full and c not in ["country","geo_loc_name","collection_date","tissue_type"]]
    full[[*safe,"focal_taxon","exact_name_match","reported_country","reported_year","sample_key","project_key"]].to_csv(OUT/"public_run_inventory.csv",index=False)
    (OUT/"receipt.json").write_text(json.dumps({"api":API,"requested_fields":chosen,"taxa":meta,
        "gate":"metadata screen only, not evidence of ancestry, fitness or selection; confirm study design, independent populations and trait linkage manually",
        "candidate_project_rows":int(summary.possible_nuclear_population_screen.sum()) if len(summary) else 0},indent=2)+"\n")
    print(json.dumps(meta),flush=True)

if __name__=="__main__":
    main()

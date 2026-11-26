import pandas as pd
from analysis.v5.sgv_public_archive_inventory import normalize,summarize
def test_independent_samples_not_runs():
    f=pd.DataFrame(dict(run_accession=["r1","r2","r3"],study_accession=["p"]*3,sample_accession=["a","a","b"],
       scientific_name=["Cirsium vulgare"]*3,library_strategy=["RAD-Seq"]*3,library_source=["GENOMIC"]*3,
       country=["France: area","France: area","USA"],collection_date=["2019-01-01","2019-01-01","2021"]))
    q=normalize(f,"Cirsium vulgare")
    s=summarize(q).iloc[0]
    assert s.runs==3 and s.unique_samples==2 and s.samples_with_year==2
    assert s.countries=="France;USA" and not s.possible_nuclear_population_screen

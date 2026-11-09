"""Export taxon-level sufficient statistics from the authenticated S1 source."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from analysis.v4.s1_s2_vif5_si10_alltraits import selected_source

ROOT = Path("work/bio18_geo")
ENV = ["chelsa_bio18", "chelsa_bio05", "chelsa_rsds_mean", "topo_elevation"]
COLOUR = ["corolla_lab_lightness", "corolla_lab_chroma", "corolla_hue_sin", "corolla_hue_cos"]

def main():
    native = pd.read_csv(ROOT/"private/decoded/observations.csv.gz", dtype={"obs_id": str}, low_memory=False)
    ext = pd.read_csv(ROOT/"private/extended-source.csv.gz", dtype={"obs_id": str})
    data = native.merge(ext, on="obs_id", validate="one_to_one")
    assert len(data) == 313350
    subset, coverage = selected_source(data, COLOUR)
    assert coverage["status"] == "SOURCE_LINKED"
    value_cols = ENV + COLOUR
    records = []
    for key, group in subset.groupby("accepted_key", sort=True):
        values = group[value_cols].to_numpy(float)
        mean = values.mean(axis=0)
        deviations = values - mean
        # Sufficient statistics for equal-taxon between and within analysis.
        records.append(dict(taxon=str(key), n=len(group), mean=mean.tolist(),
                            covariance=(deviations.T @ deviations / len(group)).tolist()))
    payload = dict(columns=value_cols, coverage=coverage, taxa=records)
    (ROOT/"public/s1_colour_moments.json").write_text(json.dumps(payload))
    print(coverage)

if __name__ == "__main__":
    main()

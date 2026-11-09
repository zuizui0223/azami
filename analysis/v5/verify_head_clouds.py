"""Check the interactive traces against the three published S2 centroids."""
import base64
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd


def array(value):
    if isinstance(value, dict):
        return np.frombuffer(base64.b64decode(value["bdata"]), dtype=value["dtype"])
    return np.asarray(value)


def main(folder):
    html = (folder / "head_clouds.html").read_text(encoding="utf-8")
    start = html.rindex("Plotly.newPlot(") + len("Plotly.newPlot(")
    decoder = json.JSONDecoder()
    args = html[start:].lstrip()
    _, end = decoder.raw_decode(args)
    traces, _ = decoder.raw_decode(args[end:].lstrip().removeprefix(",").lstrip())
    target = pd.read_csv(folder / "matched_centroids.csv")
    overlap = pd.read_csv(folder / "novel_head_overlap.csv")
    assert len(traces) == 6
    for i,(sp, data) in enumerate(target.groupby("species",sort=True)):
        data = data.set_index("trait").loc[["angle","chroma","outline"]]
        counts = overlap.query("species == @sp and native95_mass == .95 and n_traits == 3 and bandwidth == .5 and design == 'geographic_pairs'").iloc[0]
        for j,name in enumerate(["native","analogous","novel"]):
            trace = traces[3*i+j]
            xyz = np.column_stack([array(trace[k]) for k in "xyz"])
            n = counts.n_native if j==0 else counts.pairs
            assert xyz.shape == (n,3) and np.isfinite(xyz).all()
            np.testing.assert_allclose(xyz.mean(0),data[name],atol=1e-7)
            print(sp, name, len(xyz))
    assert (folder / "head_clouds.png").stat().st_size > 40000


if __name__ == "__main__":
    main(Path(sys.argv[1]))

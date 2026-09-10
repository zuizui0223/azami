import importlib.util
from pathlib import Path
import sys
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[1]/'analysis/v3'
sys.path.insert(0,str(root))
spec=importlib.util.spec_from_file_location('pilot_figures',root/'build_distribution_pilot_figures.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_disjoint_baseline_support_and_pair_reference():
    r=np.random.default_rng(4)
    samples={'A':r.normal(size=(120,1)),'B':r.normal(size=(120,1)),'small':r.normal(size=(70,1))}
    pairs=pd.DataFrame([dict(taxon_a='A',taxon_b='B',wasserstein_w2=0.,centred_w2=0.)])
    base,compare,band=m.calibrate(samples,'synthetic',.5,pairs,9)
    assert {x['taxon_name'] for x in base}=={'A','B'}
    assert len(compare)==1 and not compare[0]['exceeds_reference']
    assert compare[0]['within_reference_p95']==max(x['same_taxon_w2_p95'] for x in base)
    assert {x['factor'] for x in band}=={.75,1.,1.25}
    assert all(0<=x['density_overlap']<=1 for x in band)

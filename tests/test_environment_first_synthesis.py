import importlib.util
from pathlib import Path
import pandas as pd

def test_pinned_synthesis_preserves_families(tmp_path):
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('synthesis',root/'analysis/v3/build_environment_first_synthesis.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    source=root/'analysis_outputs/environment_first_20260910/source'
    report=mod.build(source,tmp_path)
    f=pd.read_csv(tmp_path/'trait_environment_map.csv')
    assert report['rows']==162 and len(f)==162
    for scale,file in [('among_taxon','biological_axes_among_min5.csv'),('within_taxon','biological_axes_within.csv')]:
        original=pd.read_csv(source/file)
        joined=f[f.scale==scale].merge(original,on=['construct_id','predictor'],suffixes=('_new','_old'),validate='one_to_one')
        assert (joined.q_bh_new==joined.q_bh_old).all()
    assert not f[f.scale=='within_taxon'].robustness_status.eq('full_declared_chain_pass').any()

"""Reconcile the transfer audit counts and first-stage prediction scores."""
from pathlib import Path
import numpy as np
import pandas as pd

root=Path(__file__).parent/"results"
d=root/"transfer_support"
support=pd.read_csv(d/"transfer_support.csv").set_index("species")
detail=pd.read_csv(d/"transfer_robustness.csv")
separate=pd.read_csv(d/"transfer_endpoint_cohorts.csv")
single=pd.read_csv(root/"three_stages/native_transfer_test.csv").query("model=='head_month' and trait!='outline'")
joint=pd.read_csv(root/"three_stages/native_transfer_alignment.csv").query("model=='head_month'")
original=pd.concat([single.rename(columns={"trait":"endpoint","q2_vs_no_shift":"expected"}),
                    joint.assign(endpoint="joint").rename(columns={"joint_q2":"expected"})])
matched=detail.merge(original,on=["species","endpoint"],validate="one_to_one")
assert len(matched)==6 and len(separate)==4
np.testing.assert_allclose(matched.q_zero,matched.expected)
assert (support.effective_blocks <= support.blocks+1e-12).all()
assert (separate.pairs<=separate.novel_available).all()
assert (matched.blocks_in_within_test <= matched.species.map(support.blocks)).all()
assert ((detail.loo_min<=detail.loo_max)&detail.within_signflip_p.between(0,1)).all()
print("Original scores, regional counts, and endpoint-specific cohorts verified.")


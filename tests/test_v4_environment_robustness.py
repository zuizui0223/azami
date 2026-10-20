import numpy as np
import pandas as pd

from analysis.v4.environment_robustness import abundance_ranking


def test_abundance_ranking_uses_module_eligibility_not_effects_or_other_traits():
    frame = pd.DataFrame({"accepted_key":[1,1,1,2,2,3,3],
                          "rain":[1,1,np.nan,1,1,1,1],
                          "trait":[0,0,0,1000,1000,-1000,-1000],
                          "trait__eligible":[True,True,True,True,True,False,False],
                          "other__eligible":[False]*7})
    rank = abundance_ranking(frame,["trait"],["rain"])
    assert rank.index.tolist() == [1,2]
    assert rank.tolist() == [2,2]
    frame["trait"] *= -100
    assert rank.equals(abundance_ranking(frame,["trait"],["rain"]))

"""No literature U-D and no 161-tip pseudo nuclear phylogeny."""
from analysis.v4 import t2_s1_continuous_sequence_source_coverage as m
def test_counts_are_real_source_cohorts():
    assert m.COUNTS["all251"]==(251,170066,63)
    assert m.COUNTS["its109"]==(109,162146,16)
    assert m.COUNTS["q20_posthoc110"]==(110,132762,19)
    assert m.COUNTS["source_union161"]==(161,165296,27)
    assert m.COUNTS["shared58"]==(58,129612,8)
    assert m.COUNTS["q20_only52"]==(52,3150,11)
    assert m.COUNTS["missing90"]==(90,4770,36)
    assert len(m.FIVE)==5

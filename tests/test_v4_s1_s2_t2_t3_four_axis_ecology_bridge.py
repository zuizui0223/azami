"""S1/S2/T2/T3 are one 4-axis question not independent replication."""
from analysis.v4 import s1_s2_t2_t3_four_axis_ecology_bridge as m
def test_four_axes_names_and_four_original_verified_receipts():
    assert m.AXES==("chelsa_bio18","chelsa_bio05","chelsa_rsds_mean","elevation_m")
    assert m.PANELS==("NA42_genuine_published_nuclear","Japan12_genuine_176_loci_nuclear")
    assert set(m.SHA)=={"S1_S2","T2_true_nuclear","T2_BIO18_gap","T3_extant_public_photos"}
    assert m.SHA["S1_S2"]=="080c3c66228516834f62a32e8a7849199ff6580818030de19cc038b4c73a45bb"
    assert m.SHA["T2_BIO18_gap"]=="894997bf2faa371b9664b55564c8be527aa699bfd7f379a95ceafe4dd451510b"
    assert m.SHA["T3_extant_public_photos"]=="ffb90a009b031522480f3b5618db2c7ad70050b993e6eee70f0a428e05cd20d2"

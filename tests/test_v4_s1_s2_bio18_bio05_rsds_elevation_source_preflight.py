from analysis.v4 import s1_s2_bio18_bio05_rsds_elevation_source_preflight as m
def test_new_climate_terms_are_not_original_invented_S1_data():
    assert set(m.SOURCE)=={"chelsa_bio05","chelsa_bio18","elevation_worldclim"}
    assert len(m.SOURCE["chelsa_bio05"])>=2
    assert len(m.SOURCE["chelsa_bio18"])>=2
    assert all(u.startswith("https://") for urls in m.SOURCE.values() for u in urls)
    assert m.HEADERS["Range"]=="bytes=0-16383"

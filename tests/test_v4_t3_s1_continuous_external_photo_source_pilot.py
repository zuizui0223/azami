"""Ensure T3 GenBank+GBIF real image source is NEVER botanical U/D or S1 angle."""
from analysis.v4 import t3_s1_continuous_external_photo_source_pilot as m

def test_fixed_genus_quota_and_no_fake_angles():
    assert sum(m.GENUS_QUOTAS.values())==15
    assert m.N_SOURCE==131
    assert m.canonical("Cremanthodium campanulatum (Franch.) Diels")=="Cremanthodium campanulatum"
    assert m.strict_match({"rank":"SPECIES","matchType":"EXACT","confidence":99,"scientificName":"Cremanthodium campanulatum","usageKey":100},"Cremanthodium campanulatum")==100
    assert m.strict_match({"rank":"SPECIES","matchType":"FUZZY","confidence":99,"scientificName":"Cremanthodium campanulatum","usageKey":100},"Cremanthodium campanulatum") is None
    assert m.valid_media([{"type":"StillImage","identifier":"https://example.org/a.jpg"}])
    assert not m.valid_media([{"type":"Sound","identifier":"https://example.org/a.wav"}])

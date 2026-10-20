"""Fail-closed tests for full-source WCVP eligibility and environment cache pins."""
import pandas as pd
import pytest
from analysis.v4.native_665115_pipeline import (
    add_explicit_non_species_placeholders, CLIMATE, ORIGINAL_SHA, AUTHORITY_SHA,
)

def test_all_species_names_queried_and_non_species_explicitly_unresolved():
    frame=pd.DataFrame({'taxon_name':['Cirsium vulgare','Cirsium','Cirsium vulgare vulgare','Cirsium × affine'],
                        'taxon_rank':['species','genus','variety','hybrid']})
    src={'resolutions':[{'input_name':'Cirsium vulgare','resolution_status':'resolved_unique_accepted_key',
                         'accepted_key':123,'accepted_name':'Cirsium vulgare'}]}
    updated,number,missing=add_explicit_non_species_placeholders(frame,src)
    assert (number,missing)==(1,3)
    assert len(updated['resolutions'])==4
    assert all(r['accepted_key'] is None and r['resolution_status']=='unresolved_taxon'
               for r in updated['resolutions'][1:])
    assert len(src['resolutions'])==1

def test_missing_species_name_fails_before_native_classification():
    frame=pd.DataFrame({'taxon_name':['Cirsium vulgare','Cirsium discolor'],
                        'taxon_rank':['species','species']})
    src={'resolutions':[{'input_name':'Cirsium vulgare','accepted_key':1}]}
    with pytest.raises(ValueError,match='Missing exact species'):
        add_explicit_non_species_placeholders(frame,src)

def test_nine_climate_archives_and_two_source_hashes_pinned():
    assert len(CLIMATE)==9 and len({v[0] for v in CLIMATE.values()})==9
    assert all(len(v[1])==64 for v in CLIMATE.values())
    assert len(ORIGINAL_SHA)==len(AUTHORITY_SHA)==64

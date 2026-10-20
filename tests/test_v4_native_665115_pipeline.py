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


@pytest.mark.parametrize('case', ['missing', 'network', 'changed_bytes'])
def test_private_source_download_stops_without_exposing_url(tmp_path, monkeypatch, case):
    import io
    import textwrap
    import urllib.request
    from pathlib import Path
    workflow = (Path(__file__).resolve().parents[1] /
                '.github/workflows/v4-native-665115-replay.yml').read_text()
    # Exercise the actual runner's acquisition code, with no external network.
    script = textwrap.dedent(workflow.split("python - <<'PY'\n", 1)[1].split("\n          PY", 1)[0])
    url = 'https://example.invalid/private?token=do-not-print'
    monkeypatch.setenv('AZAMI_ORIGINAL_METADATA_URL', '' if case == 'missing' else url)
    monkeypatch.chdir(tmp_path)
    def fetch(*args, **kwargs):
        if case == 'network':
            raise OSError(url)
        response = io.BytesIO(b'not the frozen source')
        response.url = url
        return response
    monkeypatch.setattr(urllib.request, 'urlopen', fetch)
    with pytest.raises(SystemExit) as exc:
        exec(compile(script, '<private-source-acquisition>', 'exec'), {})
    assert url not in str(exc.value)
    assert 'do-not-print' not in str(exc.value)
    assert not list(tmp_path.rglob('*.zip'))

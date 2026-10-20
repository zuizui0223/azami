"""Fetch authoritative name/distribution evidence for ALL 288 upstream source taxa.

Acquisition only: no photographs, traits, environment outcomes or native selection
are supplied here. This explicit user-directed run creates a dated source snapshot.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import time
from urllib.parse import urlencode
import requests
from analysis.v4.native_scope import resolve_exact

DATASET = 'f382f0ce-323a-4091-bb9f-add557f3a9a2'
API = 'https://api.gbif.org/v1'
GEOMETRY = 'https://raw.githubusercontent.com/tdwg/wgsrpd/52da7828aba9d461dd133c27b3bd7a4407161f54/geojson/level3.geojson'
GEOMETRY_SHA = 'c172bcf6aba20e19477adc60aebf0023068f0175dca8480f0760b090dcf64840'


def download(url):
    last = None
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=45, headers={'User-Agent':'azami-native-first/20261020'})
            r.raise_for_status()
            return r.content
        except requests.RequestException as exc:
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError('Source acquisition failed: ' + url) from last


def pages(url):
    records, payloads, offset = [], [], 0
    while True:
        sep = '&' if '?' in url else '?'
        payload = json.loads(download(url + sep + urlencode({'limit':100, 'offset':offset})))
        payloads.append(payload)
        rows = payload.get('results', [])
        records.extend(rows)
        if payload.get('endOfRecords', len(rows) < 100):
            return records, payloads
        if not rows:
            raise RuntimeError('Pagination stalled')
        offset += len(rows)


def run(outdir):
    outdir.mkdir(exist_ok=True, parents=True)
    namesfile = Path(__file__).with_name('native_upstream_taxa.txt')
    names = namesfile.read_text().splitlines()
    if len(names) != 288 or len(set(names)) != 288:
        raise ValueError('Expected exactly the full 288-taxon upstream source roster')
    started = datetime.now(timezone.utc).isoformat()
    meta = json.loads(download(API + '/dataset/' + DATASET))
    if meta.get('doi','').lower() != '10.15468/6h8ucr':
        raise ValueError('Wrong WCVP dataset')
    def query_name(name):
        rows, raw = pages(API + '/species/search?' + urlencode({'dataset_key':DATASET, 'q':name}))
        return resolve_exact(name, rows), raw
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(query_name, names))
    resolutions = [r[0] for r in results]
    keys = sorted({r['accepted_key'] for r in resolutions if r['accepted_key'] is not None})
    def query_distribution(key):
        rows, raw = pages(API + f'/species/{key}/distributions')
        return str(key), rows, raw
    with ThreadPoolExecutor(max_workers=4) as pool:
        distributions = list(pool.map(query_distribution, keys))
    geometry_bytes = download(GEOMETRY)
    if hashlib.sha256(geometry_bytes).hexdigest() != GEOMETRY_SHA:
        raise ValueError('TDWG geometry checksum differs')
    snapshot = {'started_utc':started, 'completed_utc':datetime.now(timezone.utc).isoformat(),
                'dataset_metadata':meta, 'resolutions':resolutions,
                'distributions':{key:rows for key,rows,raw in distributions},
                'geometry':json.loads(geometry_bytes)}
    raw = {'name_queries':{name:r[1] for name,r in zip(names,results)},
           'distribution_queries':{key:raw for key,rows,raw in distributions}}
    (outdir/'native_sources.json').write_text(json.dumps(snapshot, ensure_ascii=False))
    (outdir/'raw_gbif_responses.json').write_text(json.dumps(raw, ensure_ascii=False))
    for name in ['native_scope.py','native_upstream_taxa.txt','collect_native_sources.py']:
        shutil.copyfile(Path(__file__).with_name(name),outdir/name)
    hashes = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in outdir.iterdir() if p.is_file()}
    receipt = {'source_dataset':DATASET, 'dataset_doi':meta.get('doi'),
               'dataset_modified':meta.get('modified'), 'n_requested_taxa':len(names),
               'n_resolved_source_taxa':sum(r['accepted_key'] is not None for r in resolutions),
               'n_accepted_taxa':len(keys), 'tdwg_geometry_sha256':GEOMETRY_SHA,
               'files_sha256':hashes, 'traits_used_for_selection':False,
               'analysis_not_run':True, 'native_definition':'existing WCVP native TDWG3 record; absent records remain unlisted'}
    (outdir/'source_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out-dir',type=Path,required=True)
    run(p.parse_args().out_dir)

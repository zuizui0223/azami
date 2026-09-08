"""Plan complete native-component chunks and derive location-blind work packets."""
from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
from pathlib import Path
import sqlite3

from .protected_artifacts import new_json, require
from .reconciled_photo_schedule import STATUS, component_score
from .reconciled_stream_input import pilot_input
from .workflow import ROOT, canonical_digest, digest

CONTRACT = ROOT / 'analysis/v3/measurement_chunk_contract.json'


def contract():
    value = json.loads(CONTRACT.read_text(encoding='utf-8'))
    require(value['status'] == 'bounded_raw_measurement_chunks_authorized_no_ecological_admission'
            and value['raw_measurement_chunk_authorized'] is True
            and value['ecological_fitting_authorized'] is False
            and value['maximum_observations_per_chunk'] == 128
            and value['maximum_requests_per_chunk'] == 64, 'Unexpected chunk authority')
    return value


def validate_algorithm():
    from . import stream_original_traits as worker
    from .detect_cached_images import MODEL_SHA
    value = contract()['algorithm']
    require(hashlib.sha256(inspect.getsource(worker.measure_head).encode()).hexdigest() == value['measure_head_source_sha256']
            and canonical_digest(worker.features.specification()) == value['feature_specification_canonical_sha256']
            and MODEL_SHA == value['model_sha256'], 'Chunk measurement algorithm changed')
    decision = json.loads((ROOT/'analysis/v3/measurement_qualification_decision_20260908.json').read_text())
    require(canonical_digest(decision) == value['measurement_decision_canonical_sha256'], 'Measurement routes changed')


def group_chunks(groups, reused_ids, max_obs=128, max_requests=64):
    """Greedy consecutive groups; limits never discard or split a component."""
    require(len({g['component_id'] for g in groups}) == len(groups), 'Duplicate source component')
    require(set(reused_ids) <= {g['component_id'] for g in groups}, 'Unknown reused component')
    chunks, current = [], None
    for group in groups:
        if group['component_id'] in reused_ids:
            continue
        require(0 < group['observations'] <= max_obs and 0 <= group['requests'] <= max_requests,
                'A complete component exceeds the execution bound; preserve it and stop')
        if current is None or current['observations']+group['observations'] > max_obs or current['requests']+group['requests'] > max_requests:
            current = {'chunk_id':f'c{len(chunks):06d}', 'component_ids':[], 'observations':0, 'requests':0, 'photo_jobs':0}
            chunks.append(current)
        current['component_ids'].append(group['component_id'])
        for key in ('observations','requests','photo_jobs'):
            current[key] += group[key]
    covered = [g for chunk in chunks for g in chunk['component_ids']]
    require(len(covered) == len(set(covered)) and set(covered)|set(reused_ids) == {g['component_id'] for g in groups}, 'Source partition lost a component')
    return chunks


def build_plan(schedule: Path, reused_csv: Path, out: Path):
    spec = contract()
    require(digest(schedule) == spec['source_schedule_sha256'], 'Source schedule identity differs')
    require(digest(reused_csv) == spec['reused_pilot']['selected_observations_file_sha256'], 'Reused pilot identity differs')
    with reused_csv.open(newline='', encoding='utf-8') as handle:
        ids = [r['obs_id'] for r in csv.DictReader(handle)]
    require(len(ids) == len(set(ids)) == spec['reused_pilot']['observations'], 'Reused pilot observation count differs')
    schedule_spec = json.loads((ROOT/'analysis/v3/reconciled_photo_schedule_contract.json').read_text())
    with sqlite3.connect(schedule.resolve().as_uri()+'?mode=ro', uri=True) as db:
        execution = {k:json.loads(v) for k,v in db.execute('SELECT key,value_json FROM execution')}
        require(execution['status'] == STATUS and execution['contract_canonical_sha256'] == canonical_digest(schedule_spec), 'Schedule source contract differs')
        rows = db.execute('''WITH n AS (SELECT component_id,COUNT(*) n FROM native_observations GROUP BY component_id),
            p AS (SELECT component_id,COUNT(*) jobs,SUM(state='request_candidate_not_authorized') requests FROM photo_jobs GROUP BY component_id)
            SELECT n.component_id,n.n,p.jobs,p.requests FROM n JOIN p USING(component_id)''').fetchall()
        groups = [dict(component_id=g,observations=n,photo_jobs=p,requests=q) for g,n,p,q in rows]
        groups.sort(key=lambda r:(component_score(r['component_id'],schedule_spec['ordering_salt']),r['component_id']))
        require(sum(g['observations'] for g in groups) == db.execute('SELECT COUNT(*) FROM native_observations').fetchone()[0]
                and sum(g['photo_jobs'] for g in groups) == db.execute('SELECT COUNT(*) FROM photo_jobs').fetchone()[0], 'Unrepresented native components or photo jobs')
        db.execute('CREATE TEMP TABLE reused(obs_id TEXT PRIMARY KEY)')
        db.executemany('INSERT INTO reused VALUES (?)',[(i,) for i in ids])
        reused_ids = {r[0] for r in db.execute('SELECT DISTINCT component_id FROM native_observations JOIN reused USING(obs_id)')}
        require(sum(g['observations'] for g in groups if g['component_id'] in reused_ids) == len(ids)
                and db.execute('SELECT COUNT(*) FROM native_observations JOIN reused USING(obs_id)').fetchone()[0] == len(ids), 'Reused pilot splits a component or contains absent observations')
    chunks = group_chunks(groups,reused_ids,spec['maximum_observations_per_chunk'],spec['maximum_requests_per_chunk'])
    reused = [g for g in groups if g['component_id'] in reused_ids]
    plan = {'schema_version':1,'contract_canonical_sha256':canonical_digest(spec),'schedule_sha256':digest(schedule),
            'reused_components':reused,'chunks':chunks,'source_totals':{k:sum(g[k] for g in groups) for k in ('observations','photo_jobs','requests')},
            'source_components':len(groups),'trait_values_read':0,'ecological_fitting_authorized':False}
    plan['plan_id'] = canonical_digest(plan)
    new_json(out,plan)
    return plan


def validate_plan(plan):
    body={k:v for k,v in plan.items() if k!='plan_id'}
    require(plan['plan_id'] == canonical_digest(body) and plan['contract_canonical_sha256'] == canonical_digest(contract()), 'Plan identity differs')
    return plan


def packet_for_chunk(schedule: Path, plan: dict, chunk_id: str):
    validate_plan(plan)
    matches=[r for r in plan['chunks'] if r['chunk_id']==chunk_id]
    require(len(matches)==1,'Unknown chunk')
    chunk=matches[0]
    packet=pilot_input(schedule,plan['schedule_sha256'],128,component_ids=chunk['component_ids'])
    require(len(packet['queue'])==chunk['requests'] and len(packet['selected'])==chunk['observations']
            and packet['report']['selected_photo_jobs']==chunk['photo_jobs'],'Chunk packet source conservation differs')
    packet.update(plan_id=plan['plan_id'],chunk_id=chunk_id)
    return packet


def photo_unit_packet(packet, item):
    photo=item['photo_id']
    links=[r for r in packet['links'] if r['photo_id']==photo]
    selected=sorted({r['obs_id'] for r in links})
    require(set(selected)==set(item['obs_ids']) and selected,'Photo observation links differ')
    return {'selected':selected,'selection_scores':{k:packet['selection_scores'][k] for k in selected},
            'queue':[item],'links':links,'input_sha256':packet['input_sha256'],'schedule_sha256':packet['schedule_sha256'],
            'report':{'mode':'reconciled_photo_unit','selected_observations':len(selected),'selected_photo_links':len(links),
                      'request_candidates':1,'plan_id':packet['plan_id'],'chunk_id':packet['chunk_id'],
                      'component_id':item['component_id'],'source_reconciliation_verified':True,
                      'production_image_execution_authorized':False,'images_fetched':0,'trait_values_read':0}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--schedule',type=Path,required=True)
    parser.add_argument('--reused-pilot-selected',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    plan=build_plan(args.schedule,args.reused_pilot_selected,args.out)
    print(json.dumps({'plan_id':plan['plan_id'],'chunks':len(plan['chunks']),'source_totals':plan['source_totals'],'reused_components':len(plan['reused_components'])}))


if __name__=='__main__':
    main()

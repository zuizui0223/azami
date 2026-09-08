"""Reversible observation/module views of independently verified raw streams.

All native source records and all 27 endpoint slots stay visible. This is a
measurement inventory, not an ecological fit or a new measurement qualification.
Raw archives remain the authoritative per-head numerical evidence.
"""
from __future__ import annotations

from collections import defaultdict
import argparse
import json
from pathlib import Path
import sqlite3

from .archived_measurement import DECISION, exact_decision_view, verify_chunk
from .build_dependence_groups import Components
from .build_observation_measurements import aggregate_heads, eligible_endpoints, finite, mean
from .image_features import registry
from .measurement_chunks import contract
from .protected_artifacts import new_json, require
from .recover_native_source_authority import private_directory
from .verify_original_stream import CONDITIONS, number, read_rows, verify
from .workflow import canonical_digest, digest, text_digest

SPECIFICATION={
    'version':'v3_verified_raw_stream_observation_views_v1',
    'source_denominator':'Every observation and photo link in the exact pinned native source schedule; never a union of successful photos.',
    'head_to_photo':'Arithmetic mean of eligible heads. Raw finite and joint-QC means remain separate from 14-route operational means.',
    'photo_to_observation':'Equal weight per distinct EXIF-oriented pixel identity, not per photo alias or detected head.',
    'joint_modules':'Within each module, all coordinates use the SAME operationally eligible heads in each image and the SAME images in each observation. Colour retains its raw coordinates; this does not choose a model representation.',
    'bbox':'Orientation and gross-shape modules retain baseline and all four fixed perturbations on the same head/image set. These are technical envelopes, not calibrated error variances.',
    'quality':'Module/endpoint-matched pixel-size and sharpness means are unavailable if any contributing head or image lacks them.',
    'incomplete':'Observed partial values remain descriptive fields. Final operational values require every request-candidate photo linked to the observation to have a verified terminal record. Rights/metadata-blocked photos are retained, not requested.',
    'identity_conflict':'Identical pixels require identical detector and full head numerical summaries. Conflicting processing outcomes stop the build, never choose the best result.',
    'dependence':'Retain original source component IDs and derived joins from exact bytes/pixels. These groups are not biological individuals; no approximate scene matching is claimed.',
    'colour':'Preserve legacy operational colour and same-head floral/context candidate contrasts separately. Uniform-floral colour still needs its own image-only qualification; neither is pigment concentration.',
    'interpretation':'This is not ecological admission, botanical accuracy, independent same-head replication or separation of photographic and biological variance.',
}


def definition():
    decision=json.loads(DECISION.read_text(encoding='utf-8'))
    admitted=set(decision['next_measurement_execution']['original_stream_required_endpoints'])
    modules=defaultdict(list)
    for row in decision['endpoints']:
        if row['endpoint_id'] in admitted:
            modules[row['module']].append(row['endpoint_id'])
    require(len(admitted)==14 and sum(map(len,modules.values()))==14,'Operational route inventory differs')
    return admitted,{k:sorted(v) for k,v in sorted(modules.items())}


def read_photos(directory):
    """Read only after the caller has independently verified the whole stream."""
    photos={r['photo_id']:{'transfer':r,'heads':{},'bbox':defaultdict(dict)}
            for r in read_rows(directory/'transfer_private.csv')}
    for line in (directory/'photo_detection_private.jsonl').read_text(encoding='utf-8').splitlines():
        row=json.loads(line); photos[row['photo_id']]['detection']=row
    for line in (directory/'head_diagnostics_private.jsonl').read_text(encoding='utf-8').splitlines():
        row=json.loads(line)
        photos[row['photo_id']]['heads'][str(row['head_index'])]={'endpoints':{},'diagnostics':row['diagnostics']}
    for row in read_rows(directory/'endpoint_measurements_private.csv'):
        record={k:v for k,v in row.items() if k not in ('photo_id','head_index')}
        for key in ('value','original','mirror','mirror_abs_difference'): record[key]=number(record[key])
        record['primary_measurement_eligible']=record['primary_measurement_eligible']=='True'
        photos[row['photo_id']]['heads'][row['head_index']]['endpoints'][row['endpoint_id']]=record
    for row in read_rows(directory/'bbox_measurements_private.csv'):
        record={k:v for k,v in row.items() if k not in ('photo_id','head_index','endpoint_id','condition')}
        for key in ('value','original','mirror','mirror_abs_difference'): record[key]=number(record[key])
        photos[row['photo_id']]['bbox'][row['head_index']][(row['endpoint_id'],row['condition'])]=record
    return photos


def quality(heads,key):
    values=[h['diagnostics'].get(key) for h in heads]
    return mean(values) if values and all(finite(v) for v in values) else None


def summarize(photo):
    """Numerical definitions independent of taxon, location or environment."""
    admitted,modules=definition()
    indices=sorted(photo['heads'],key=int)
    heads=[photo['heads'][i] for i in indices]
    rows,_=aggregate_heads(heads)
    qc=[eligible_endpoints(h['endpoints']) for h in heads]
    operational=[{k:ok[k] and k in admitted and h['endpoints'][k]['primary_measurement_eligible']
                  for k in ok} for h,ok in zip(heads,qc)]
    for row in rows:
        key=row['endpoint_id']; selected=[h for h,ok in zip(heads,operational) if ok[key]]
        row.update(operational_mean=mean([h['endpoints'][key]['value'] for h in selected]),
                   n_operational_heads=len(selected),
                   operational_size=quality(selected,'head_min_dimension_px'),
                   operational_sharpness=quality(selected,'head_laplacian_variance'))
    joint=[]
    for module,keys in modules.items():
        selected=[(i,h) for i,h,ok in zip(indices,heads,operational) if all(ok[k] for k in keys)]
        conditions=sorted(CONDITIONS) if module in ('orientation','gross_shape') else ['baseline']
        for condition in conditions:
            for key in keys:
                if condition=='baseline':
                    values=[h['endpoints'][key]['value'] for i,h in selected]
                else:
                    records=[photo['bbox'][i].get((key,condition),{}) for i,h in selected]
                    require(all(r.get('status')=='usable' and finite(r.get('value')) for r in records),
                            'An operational module lacks a complete bbox condition')
                    values=[r['value'] for r in records]
                chosen=[h for i,h in selected]
                joint.append({'module':module,'endpoint_id':key,'condition':condition,'n_heads':len(selected),
                              'value':mean(values),'size':quality(chosen,'head_min_dimension_px'),
                              'sharpness':quality(chosen,'head_laplacian_variance')})
    colour=[h for h,ok in zip(heads,operational) if ok['corolla_lab_chroma'] and ok['corolla_lab_lightness']]
    _,pairs=aggregate_heads(colour)
    return rows,joint,pairs


class ObservationViews:
    """Sparse saved photo summaries plus a full-source all27 observation view."""
    def __init__(self,schedule: Path,out: Path,*,expected_schedule_sha=None):
        self.out=private_directory(out)
        require(not self.out.exists(),'Preserve earlier observation views')
        expected_schedule_sha=expected_schedule_sha or contract()['source_schedule_sha256']
        require(digest(schedule)==expected_schedule_sha,'Full native source schedule bytes differ')
        self.out.mkdir(parents=True)
        self.db=sqlite3.connect(self.out/'observation_views_private.sqlite',uri=True)
        self.db.execute('PRAGMA foreign_keys=ON')
        self.db.execute('ATTACH DATABASE ? AS source',(schedule.resolve().as_uri()+'?mode=ro',))
        self.db.executescript('''
            CREATE TABLE source_observations AS SELECT * FROM source.native_observations;
            CREATE UNIQUE INDEX source_obs ON source_observations(obs_id);
            CREATE TABLE source_photos AS SELECT photo_id,component_id,state FROM source.photo_jobs;
            CREATE UNIQUE INDEX source_photo ON source_photos(photo_id);
            CREATE TABLE source_links AS SELECT * FROM source.native_links;
            CREATE UNIQUE INDEX source_link ON source_links(obs_id,photo_id);
            CREATE INDEX source_link_photo ON source_links(photo_id);
            CREATE TABLE endpoints(endpoint_id TEXT PRIMARY KEY,unit TEXT,operational_route INTEGER);
            CREATE TABLE inputs(input_id TEXT PRIMARY KEY,receipt_json TEXT);
            CREATE TABLE photo_results(photo_id TEXT PRIMARY KEY,source_bytes TEXT,pixel_id TEXT,status TEXT,
                n_heads INTEGER,fingerprint TEXT,source_identity_sha TEXT);
            CREATE INDEX result_pixels ON photo_results(pixel_id);
            CREATE TABLE photo_inputs(photo_id TEXT,input_id TEXT,PRIMARY KEY(photo_id,input_id));
            CREATE TABLE pixel_identity(pixel_id TEXT PRIMARY KEY,fingerprint TEXT,representative_photo TEXT);
            CREATE TABLE photo_endpoints(photo_id TEXT,endpoint_id TEXT,raw_mean REAL,qc_mean REAL,op_mean REAL,
                n_finite INTEGER,n_qc INTEGER,n_op INTEGER,size REAL,sharpness REAL,PRIMARY KEY(photo_id,endpoint_id));
            CREATE TABLE photo_modules(photo_id TEXT,module TEXT,endpoint_id TEXT,condition TEXT,value REAL,n_heads INTEGER,
                size REAL,sharpness REAL,PRIMARY KEY(photo_id,module,endpoint_id,condition));
            CREATE TABLE photo_colour_pairs(photo_id TEXT,context_kind TEXT,statistic TEXT,n_pairs INTEGER,
                legacy REAL,floral REAL,context REAL,contrast REAL,PRIMARY KEY(photo_id,context_kind,statistic));
        ''')
        self.db.execute('DETACH DATABASE source')
        admitted,_=definition()
        self.db.executemany('INSERT INTO endpoints VALUES (?,?,?)',[(r['endpoint_id'],r['unit'],int(r['endpoint_id'] in admitted)) for r in registry()])
        self.context={'specification':SPECIFICATION,'source_schedule_sha256':expected_schedule_sha,
                      'decision_canonical_sha256':canonical_digest(json.loads(DECISION.read_text(encoding='utf-8'))),
                      'implementation_sha256_text_lf':text_digest(Path(__file__))}
        new_json(self.out/'execution_contract.json',self.context)
        self.failed=False

    def add_directory(self,directory,input_id):
        """Internal callback: only receive independently verified directory bytes."""
        for photo_id,photo in read_photos(directory).items():
            source=self.db.execute('SELECT state FROM source_photos WHERE photo_id=?',(photo_id,)).fetchone()
            require(source is not None and source[0]=='request_candidate_not_authorized','Measured photo outside requested source universe')
            links={r[0] for r in self.db.execute('SELECT obs_id FROM source_links WHERE photo_id=?',(photo_id,))}
            require(links==set(photo['detection']['obs_ids']),'Measured photo links differ from the full source schedule')
            transfer=photo['transfer']; pixel=transfer['oriented_rgb_pixel_sha256'] or None
            fingerprint=canonical_digest({'heads':photo['heads'],'bbox':{i:sorted((k,v) for k,v in d.items()) for i,d in photo['bbox'].items()},
                                          'detector_results':photo['detection']['detections'],'status':transfer['status']})
            identity=canonical_digest({k:v for k,v in photo['detection'].items() if k not in ('transfer','status','detections')})
            record=(photo_id,transfer['source_byte_sha256'] or None,pixel,transfer['status'],len(photo['heads']),fingerprint,identity)
            prior=self.db.execute('SELECT * FROM photo_results WHERE photo_id=?',(photo_id,)).fetchone()
            if prior is not None:
                require(prior==record,'Conflicting repeated photo outcome; do not choose a version')
                self.db.execute('INSERT OR IGNORE INTO photo_inputs VALUES (?,?)',(photo_id,input_id))
                continue
            if pixel:
                prior=self.db.execute('SELECT fingerprint FROM pixel_identity WHERE pixel_id=?',(pixel,)).fetchone()
                require(prior is None or prior[0]==fingerprint,'Identical pixels have inconsistent processing summaries')
                self.db.execute('INSERT OR IGNORE INTO pixel_identity VALUES (?,?,?)',(pixel,fingerprint,photo_id))
            self.db.execute('INSERT INTO photo_results VALUES (?,?,?,?,?,?,?)',record)
            self.db.execute('INSERT INTO photo_inputs VALUES (?,?)',(photo_id,input_id))
            rows,joint,pairs=summarize(photo)
            self.db.executemany('INSERT INTO photo_endpoints VALUES (?,?,?,?,?,?,?,?,?,?)',[
                (photo_id,r['endpoint_id'],r['raw_mean'],r['eligible_mean'],r['operational_mean'],r['n_finite_heads'],
                 r['n_eligible_heads'],r['n_operational_heads'],r['operational_size'],r['operational_sharpness']) for r in rows])
            self.db.executemany('INSERT INTO photo_modules VALUES (?,?,?,?,?,?,?,?)',[
                (photo_id,r['module'],r['endpoint_id'],r['condition'],r['value'],r['n_heads'],r['size'],r['sharpness']) for r in joint])
            self.db.executemany('INSERT INTO photo_colour_pairs VALUES (?,?,?,?,?,?,?,?)',[
                (photo_id,r['context_kind'],r['statistic'],r['n_pairs'],r['legacy_mean'],r['floral_mean'],r['context_mean'],r['contrast_mean']) for r in pairs])

    def add_pilot(self,directory: Path):
        try:
            report=json.loads((directory/'original_stream_report.json').read_text(encoding='utf-8'))
            decision=exact_decision_view(self.out,report['measurement_decision_sha256'])
            checked=verify(directory,decision)
            execution=json.loads((directory/'execution_contract.json').read_text(encoding='utf-8'))
            require(execution['reconciled_schedule_sha256']==self.context['source_schedule_sha256'],'Pilot source schedule differs')
            input_id='pilot:'+digest(directory/'original_stream_report.json')
            with self.db:
                self.add_directory(directory,input_id)
                self.db.execute('INSERT INTO inputs VALUES (?,?)',(input_id,json.dumps(checked,sort_keys=True)))
        except BaseException:
            self.failed=True
            raise

    def add_archive(self,bundle,asset,batch,chunk_id,packet_path):
        try:
            input_id='archive:'+asset['bundle_sha256']
            require(self.db.execute('SELECT 1 FROM inputs WHERE input_id=?',(input_id,)).fetchone() is None,'Archive already collected')
            with self.db:
                report=verify_chunk(bundle,asset,batch,chunk_id,self.out/f'archive-{chunk_id}',packet_path=packet_path,
                                    on_unit=lambda directory,packet,checked:self.add_directory(directory,input_id))
                self.db.execute('INSERT INTO inputs VALUES (?,?)',(input_id,json.dumps(report,sort_keys=True)))
        except BaseException:
            self.failed=True
            raise

    def finish(self):
        require(not self.failed,'An input failed; do not promote a partial import')
        self.db.executescript('''
            CREATE TABLE observation_support AS
            SELECT o.obs_id,o.component_id,COUNT(p.photo_id) n_source_photos,
                SUM(CASE WHEN p.state='request_candidate_not_authorized' THEN 1 ELSE 0 END) n_requests,
                SUM(CASE WHEN p.state='request_candidate_not_authorized' AND r.photo_id IS NULL THEN 1 ELSE 0 END) n_pending,
                COUNT(r.photo_id) n_terminal,SUM(CASE WHEN r.status='success' AND r.n_heads=0 THEN 1 ELSE 0 END) n_no_detection,
                SUM(CASE WHEN r.n_heads>0 THEN 1 ELSE 0 END) n_detected,
                SUM(CASE WHEN r.status IS NOT NULL AND r.status!='success' THEN 1 ELSE 0 END) n_failed
            FROM source_observations o LEFT JOIN source_links l USING(obs_id)
                LEFT JOIN source_photos p ON l.photo_id=p.photo_id LEFT JOIN photo_results r ON p.photo_id=r.photo_id GROUP BY o.obs_id;
            CREATE UNIQUE INDEX observation_support_id ON observation_support(obs_id);
            CREATE TABLE observation_images AS SELECT DISTINCT l.obs_id,r.pixel_id,i.representative_photo photo_id
                FROM source_links l JOIN photo_results r USING(photo_id) JOIN pixel_identity i USING(pixel_id);
            CREATE UNIQUE INDEX observation_pixels ON observation_images(obs_id,pixel_id);
            CREATE TABLE observation_endpoint_partial AS SELECT i.obs_id,e.endpoint_id,
                AVG(e.raw_mean) partial_raw_mean,AVG(e.qc_mean) partial_qc_mean,AVG(e.op_mean) partial_op_mean,
                COUNT(e.op_mean) n_op_images,COUNT(e.raw_mean) n_finite_images,
                CASE WHEN COUNT(e.size)=COUNT(e.op_mean) THEN AVG(e.size) END size,
                CASE WHEN COUNT(e.sharpness)=COUNT(e.op_mean) THEN AVG(e.sharpness) END sharpness
                FROM observation_images i JOIN photo_endpoints e USING(photo_id) GROUP BY i.obs_id,e.endpoint_id;
            CREATE UNIQUE INDEX observation_endpoint_partial_key ON observation_endpoint_partial(obs_id,endpoint_id);
            CREATE VIEW observation_endpoints AS SELECT s.obs_id,s.component_id,e.endpoint_id,e.unit,e.operational_route,
                s.n_pending,s.n_source_photos,s.n_requests,s.n_terminal,p.partial_raw_mean,p.partial_qc_mean,p.partial_op_mean,
                COALESCE(p.n_op_images,0) n_op_images,COALESCE(p.n_finite_images,0) n_finite_images,
                CASE WHEN s.n_pending=0 THEN p.partial_raw_mean END raw_mean,
                CASE WHEN s.n_pending=0 AND e.operational_route=1 THEN p.partial_op_mean END operational_mean,
                CASE WHEN s.n_pending=0 AND e.operational_route=1 THEN p.size END size,
                CASE WHEN s.n_pending=0 AND e.operational_route=1 THEN p.sharpness END sharpness,
                CASE WHEN s.n_pending>0 THEN 'pending_source_photos'
                     WHEN e.operational_route=0 THEN 'held_measurement_route'
                     WHEN p.n_op_images>0 THEN 'operationally_usable_not_ecologically_admitted'
                     WHEN s.n_requests=0 THEN 'source_rights_or_metadata_blocked'
                     WHEN s.n_detected>0 THEN 'qc_unusable'
                     WHEN s.n_no_detection=s.n_terminal THEN 'no_detection'
                     ELSE 'transfer_or_measurement_failed' END status
                FROM observation_support s CROSS JOIN endpoints e LEFT JOIN observation_endpoint_partial p
                    ON s.obs_id=p.obs_id AND e.endpoint_id=p.endpoint_id;
            CREATE TABLE observation_modules AS SELECT i.obs_id,m.module,m.endpoint_id,m.condition,
                COUNT(m.value) n_images,AVG(m.value) partial_value,
                CASE WHEN s.n_pending=0 THEN AVG(m.value) END value,
                CASE WHEN s.n_pending=0 AND COUNT(m.size)=COUNT(m.value) THEN AVG(m.size) END size,
                CASE WHEN s.n_pending=0 AND COUNT(m.sharpness)=COUNT(m.value) THEN AVG(m.sharpness) END sharpness,
                s.n_pending
                FROM observation_images i JOIN photo_modules m USING(photo_id) JOIN observation_support s USING(obs_id)
                GROUP BY i.obs_id,m.module,m.endpoint_id,m.condition;
            CREATE UNIQUE INDEX observation_module_key ON observation_modules(obs_id,module,endpoint_id,condition);
            CREATE TABLE observation_colour_pairs AS SELECT i.obs_id,p.context_kind,p.statistic,
                COUNT(p.floral) n_images,AVG(p.legacy) partial_legacy,AVG(p.floral) partial_floral,
                AVG(p.context) partial_context,AVG(p.contrast) partial_contrast,s.n_pending
                FROM observation_images i JOIN photo_colour_pairs p USING(photo_id) JOIN observation_support s USING(obs_id)
                GROUP BY i.obs_id,p.context_kind,p.statistic;
        ''')
        groups=Components()
        for component, in self.db.execute('SELECT DISTINCT component_id FROM source_observations'): groups.find(component)
        for column in ('source_bytes','pixel_id'):
            owners={}
            for identity,component in self.db.execute(f'SELECT r.{column},p.component_id FROM photo_results r JOIN source_photos p USING(photo_id) WHERE r.{column} IS NOT NULL'):
                if identity in owners: groups.union(owners[identity],component)
                else: owners[identity]=component
        self.db.execute('CREATE TABLE derived_components(source_component TEXT PRIMARY KEY,derived_component TEXT)')
        self.db.executemany('INSERT INTO derived_components VALUES (?,?)',[(c,groups.find(c)) for c in sorted(groups.parents)])
        self.db.commit()
        scalar=lambda sql:self.db.execute(sql).fetchone()[0]
        report={'status':'PARTIAL_RAW_STREAM_OBSERVATION_VIEWS_BUILT_NO_ECOLOGY',**self.context,
                'source_observations':scalar('SELECT COUNT(*) FROM source_observations'),
                'source_photo_jobs':scalar('SELECT COUNT(*) FROM source_photos'),
                'source_links':scalar('SELECT COUNT(*) FROM source_links'),
                'verified_photo_units':scalar('SELECT COUNT(*) FROM photo_results'),
                'distinct_measured_pixel_identities':scalar('SELECT COUNT(*) FROM pixel_identity'),
                'detected_heads':scalar('SELECT COALESCE(SUM(n_heads),0) FROM photo_results'),
                'complete_request_observations':scalar('SELECT COUNT(*) FROM observation_support WHERE n_pending=0 AND n_requests>0'),
                'pending_observations':scalar('SELECT COUNT(*) FROM observation_support WHERE n_pending>0'),
                'source_only_blocked_observations':scalar('SELECT COUNT(*) FROM observation_support WHERE n_requests=0'),
                'retained_endpoints':27,'operational_routes':14,'held_routes':13,
                'source_images_persisted':0,'local_image_requests':0,'environment_values_read':0,
                'ecological_fitting_authorized':False,'ecological_models_executed':0}
        require(report['source_observations']==report['complete_request_observations']+report['pending_observations']+report['source_only_blocked_observations'],
                'Observation denominator does not reconcile')
        self.db.close()
        report['database_sha256']=digest(self.out/'observation_views_private.sqlite')
        new_json(self.out/'public_report.json',report)
        return report


def build(inputs: Path,out: Path):
    """Replay a private list of saved inputs; no network or image model loaded."""
    spec=json.loads(inputs.read_text(encoding='utf-8'))
    require(spec['source_schedule_sha256']==contract()['source_schedule_sha256'],'Input list changes the full source authority')
    builder=ObservationViews(Path(spec['source_schedule']),out)
    new_json(builder.out/'input_list_identity.json',{'input_list_sha256':digest(inputs)})
    try:
        if spec.get('pilot'):
            builder.add_pilot(Path(spec['pilot']))
        for row in spec['archives']:
            batch=json.loads(Path(row['batch']).read_text(encoding='utf-8'))
            builder.add_archive(Path(row['bundle']),row['asset'],batch,row['chunk_id'],Path(row['packet']))
        return builder.finish()
    except BaseException as error:
        builder.db.close()
        new_json(builder.out/'incomplete_build.json',{'status':'INCOMPLETE_DO_NOT_USE','error_type':type(error).__name__})
        raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    try:
        print(json.dumps(build(args.inputs,args.out)))
    except Exception as error:
        print(json.dumps({'status':'OBSERVATION_VIEWS_INCOMPLETE','error_type':type(error).__name__}))
        raise SystemExit(1) from None


if __name__=='__main__':
    main()

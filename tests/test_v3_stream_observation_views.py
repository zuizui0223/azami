import copy
import json
import sqlite3

import pytest

from analysis.v3 import cloud_measurement_chunks as cloud
from analysis.v3 import stream_observation_views as views
from analysis.v3 import stream_original_traits as worker
from analysis.v3.image_features import registry
from analysis.v3.build_observation_measurements import COMPOSITION
from analysis.v3.protected_artifacts import new_json
from analysis.v3.verify_original_stream import CONDITIONS
from analysis.v3.workflow import digest
from test_v3_cloud_measurement_chunks import prepared
from test_v3_original_stream_traits import offline_run


def source_file(path,extra_pending=False):
    with sqlite3.connect(path) as db:
        db.executescript('''CREATE TABLE native_observations(obs_id TEXT PRIMARY KEY,component_id TEXT,expected_photo_count INTEGER);
            CREATE TABLE photo_jobs(photo_id TEXT PRIMARY KEY,component_id TEXT,state TEXT);
            CREATE TABLE native_links(obs_id TEXT,photo_id TEXT);''')
        db.executemany('INSERT INTO native_observations VALUES (?,?,?)',[('1','component-one',2+int(extra_pending)),('2','component-two',1),('3','component-three',1)])
        db.executemany('INSERT INTO photo_jobs VALUES (?,?,?)',[
            ('11','component-one','request_candidate_not_authorized'),('12','component-one','request_candidate_not_authorized'),
            ('99','component-two','request_candidate_not_authorized'),('88','component-three','license_unavailable')])
        db.executemany('INSERT INTO native_links VALUES (?,?)',[('1','11'),('1','12'),('2','99'),('3','88')])
        if extra_pending:
            db.execute("INSERT INTO photo_jobs VALUES ('77','component-one','request_candidate_not_authorized')")
            db.execute("INSERT INTO native_links VALUES ('1','77')")
    return path


def import_chunk(prepared,tmp_path,monkeypatch,extra_pending=False):
    packet,store,fixture=prepared
    # The shared fixture intentionally increments values between calls. Reset
    # it per image here so identical pixel aliases have identical measurements.
    original=worker.run_photo_unit
    def repeatable(*a,**k):
        fixture.calls['calls']=0
        return original(*a,**k)
    monkeypatch.setattr(worker,'run_photo_unit',repeatable)
    cloud_report=cloud.execute_packet(packet,store,tmp_path/'cloud',fixture.weights)
    pinned=tmp_path/'packet.json'; new_json(pinned,packet)
    batch={'plan_id':packet['plan_id'],'chunks':{packet['chunk_id']:{'packet_sha256':digest(pinned)}}}
    source=source_file(tmp_path/'source.sqlite',extra_pending)
    builder=views.ObservationViews(source,tmp_path/'views',expected_schedule_sha=digest(source))
    builder.add_archive(tmp_path/'cloud/completed/result.zip',cloud_report['output_asset'],batch,packet['chunk_id'],pinned)
    return builder


def test_full_source_and_all27_survive_sparse_verified_import(prepared,tmp_path,monkeypatch):
    builder=import_chunk(prepared,tmp_path,monkeypatch)
    report=builder.finish()
    assert report['source_observations']==3 and report['verified_photo_units']==2
    assert report['complete_request_observations']==report['pending_observations']==report['source_only_blocked_observations']==1
    assert report['distinct_measured_pixel_identities']==1 # aliases of the same decoded test image
    with sqlite3.connect(tmp_path/'views/observation_views_private.sqlite') as db:
        assert db.execute('SELECT COUNT(*) FROM observation_endpoints').fetchone()[0]==3*27
        assert db.execute("SELECT COUNT(*) FROM observation_endpoints WHERE obs_id='1' AND status='held_measurement_route'").fetchone()[0]==13
        assert db.execute("SELECT MAX(n_op_images) FROM observation_endpoints WHERE obs_id='1'").fetchone()[0]==1
        assert db.execute("SELECT COUNT(*) FROM observation_endpoints WHERE obs_id='2' AND status='pending_source_photos'").fetchone()[0]==27
    assert report['ecological_fitting_authorized'] is False


def test_partial_observation_is_not_a_failure_or_a_final_mean(prepared,tmp_path,monkeypatch):
    builder=import_chunk(prepared,tmp_path,monkeypatch,extra_pending=True)
    builder.finish()
    with sqlite3.connect(tmp_path/'views/observation_views_private.sqlite') as db:
        row=db.execute("SELECT n_pending,partial_raw_mean,operational_mean,status FROM observation_endpoints WHERE obs_id='1' AND endpoint_id='corolla_lab_chroma'").fetchone()
        assert row[0]==1 and row[1] is not None and row[2] is None and row[3]=='pending_source_photos'
        assert db.execute("SELECT COUNT(*) FROM observation_modules WHERE obs_id='1' AND value IS NOT NULL").fetchone()[0]==0


def synthetic_photo(values):
    admitted,_=views.definition()
    heads={}; bbox={}
    for i,value in enumerate(values):
        endpoints={r['endpoint_id']:{'value':value,'status':'usable','primary_measurement_eligible':r['endpoint_id'] in admitted}
                   for r in registry()}
        for key in COMPOSITION: endpoints[key]['value']=0.25
        heads[str(i)]={'endpoints':endpoints,'diagnostics':{'head_min_dimension_px':100+i,'head_laplacian_variance':200+i}}
        bbox[str(i)]={(key,condition):{'value':value+offset,'status':'usable'}
                     for key in views.definition()[1]['gross_shape']+views.definition()[1]['orientation']
                     for offset,condition in enumerate(sorted(CONDITIONS))}
    return {'heads':heads,'bbox':bbox}


def test_same_head_joint_colour_and_shape_support_and_bbox_values():
    photo=synthetic_photo([2,4,9])
    photo['heads']['1']['endpoints']['corolla_lab_lightness']['primary_measurement_eligible']=False
    photo['heads']['2']['endpoints']['capitulum_outline_circularity']['primary_measurement_eligible']=False
    rows,joint,pairs=views.summarize(photo)
    colour=[r for r in joint if r['module']=='visible_colour']
    assert {r['n_heads'] for r in colour}=={2}
    chroma=next(r for r in colour if r['endpoint_id']=='corolla_lab_chroma')
    assert chroma['value']==5.5 # same heads as the full module, not the endpoint's three heads
    shape=[r for r in joint if r['module']=='gross_shape']
    assert len(shape)==4*5 and {r['n_heads'] for r in shape}=={2}
    for condition in CONDITIONS:
        values={r['value'] for r in shape if r['condition']==condition}
        offset=sorted(CONDITIONS).index(condition)
        assert values=={3 if condition=='baseline' else 3+offset}
    assert next(r for r in rows if r['endpoint_id']=='bract_projection_maximum')['operational_mean'] is None


def test_missing_matched_quality_is_not_averaged_over_a_different_subset():
    photo=synthetic_photo([2,4])
    photo['heads']['1']['diagnostics']['head_min_dimension_px']=None
    rows,joint,_=views.summarize(photo)
    assert all(r['size'] is None for r in joint)
    assert all(r['operational_size'] is None for r in rows)
    assert next(r for r in rows if r['endpoint_id']=='corolla_lab_chroma')['operational_sharpness']==200.5


def test_missing_bbox_is_rejected_not_relabelled_as_low_uncertainty():
    photo=synthetic_photo([2])
    photo['bbox']['0'].pop(('orientation_image_vertical_angle','bbox_left_5pct'))
    with pytest.raises(ValueError,match='complete bbox'):
        views.summarize(photo)


def test_same_pixel_conflict_stops_instead_of_choosing_better(prepared,tmp_path,monkeypatch):
    builder=import_chunk(prepared,tmp_path,monkeypatch)
    directory=tmp_path/'cloud/new_units/u0000'
    photos=views.read_photos(directory)
    photo=next(iter(photos.values()))
    photo['transfer']['photo_id']='99'
    photo['detection']['obs_ids']=['2']
    key=next(iter(photo['heads']))
    photo['heads'][key]['endpoints']['corolla_lab_chroma']['value']=999
    monkeypatch.setattr(views,'read_photos',lambda _: {'99':photo})
    with pytest.raises(ValueError,match='inconsistent processing'):
        builder.add_directory(directory,'conflict')
    builder.db.close()


def test_equal_image_weight_not_equal_head_weight(tmp_path):
    source=source_file(tmp_path/'source.sqlite')
    builder=views.ObservationViews(source,tmp_path/'views',expected_schedule_sha=digest(source))
    # Exercise the real SQL aggregation seam: two images, with 1 and 9 heads.
    for photo,pixel,value,n in [('11','pixels-a',10.0,1),('12','pixels-b',30.0,9)]:
        builder.db.execute('INSERT INTO photo_results VALUES (?,?,?,?,?,?,?)',(photo,'bytes-'+photo,pixel,'success',n,'f','s'))
        builder.db.execute('INSERT INTO pixel_identity VALUES (?,?,?)',(pixel,'f',photo))
        builder.db.execute('INSERT INTO photo_endpoints VALUES (?,?,?,?,?,?,?,?,?,?)',(photo,'corolla_lab_chroma',value,value,value,n,n,n,100,200))
        builder.db.execute('INSERT INTO photo_modules VALUES (?,?,?,?,?,?,?,?)',(photo,'visible_colour','corolla_lab_chroma','baseline',value,n,100,200))
    builder.finish()
    with sqlite3.connect(tmp_path/'views/observation_views_private.sqlite') as db:
        assert db.execute("SELECT operational_mean FROM observation_endpoints WHERE obs_id='1' AND endpoint_id='corolla_lab_chroma'").fetchone()[0]==20
        assert db.execute("SELECT value FROM observation_modules WHERE obs_id='1'").fetchone()[0]==20

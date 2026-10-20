import json
import numpy as np
import pytest
from PIL import Image
from analysis.v4.upstream import digest
from analysis.v4.run_upstream_acquisition import decode
from analysis.v4.measure_upstream_anatomy import run, point, polygon_mask

def fixture(tmp_path, mode='RGB'):
    source=tmp_path/'source';(source/'originals').mkdir(parents=True)
    path=source/'originals'/'1.png'
    Image.new(mode,(40,30), 200 if mode=='L' else (220,140,210)).save(path)
    record=dict(photo_id=1,obs_id=1,status='downloaded_original',original_path='originals/1.png',sha256=digest(path),current_taxon_name='test species')
    (source/'acquisition_ledger.json').write_text(json.dumps([record]))
    annotation=dict(photo_id=1,original_sha256=digest(path),review_status='assistant_visual_proposal',flowering_stage='anthesis',floral_polygons=[[[.1,.1],[.9,.1],[.9,.9],[.1,.9]]])
    annotations=tmp_path/'annotations.json';annotations.write_text(json.dumps([annotation]))
    return source,annotations,path

def test_monochrome_is_not_white_flower(tmp_path):
    source,annotations,path=fixture(tmp_path,'L')
    image,profile=decode(path)
    assert image.mode=='RGB'
    assert profile=='monochrome_source_no_colour_inference'
    result=run(source,annotations,tmp_path/'result')
    assert result['results'][0]['colour']['status']=='source_colour_unassessable'
    assert 'lab_L_median' not in result['results'][0]['colour']

def test_invalid_profile_preserves_geometry_not_calibrated_colour(tmp_path):
    _,_,path=fixture(tmp_path)
    Image.new('RGB',(30,20),(100,140,200)).save(path,icc_profile=b'not_an_ICC_profile')
    image,profile=decode(path)
    assert image.size==(30,20)
    assert profile=='invalid_ICC_colour_unassessable'

def test_actual_pixels_and_missing_landmarks_stay_separate(tmp_path):
    source,annotations,_=fixture(tmp_path)
    result=run(source,annotations,tmp_path/'result')
    a=result['results'][0]
    assert a['colour']['status']=='measured_image_colour'
    assert a['presentation']['status']=='missing_landmarks'
    assert a['review_status']=='assistant_visual_proposal'
    assert not a['eligible_for_primary_ecology']
    assert not result['v3_phenotype_csv_used']

def test_annotation_must_match_exact_image_bytes(tmp_path):
    source,annotations,_=fixture(tmp_path)
    a=json.loads(annotations.read_text());a[0]['original_sha256']='0'*64;annotations.write_text(json.dumps(a))
    with pytest.raises(ValueError,match='byte mismatch'):
        run(source,annotations,tmp_path/'result')

def test_no_floral_annotation_is_not_zero_colour(tmp_path):
    source,annotations,_=fixture(tmp_path)
    a=json.loads(annotations.read_text());a[0]['floral_polygons']=[];annotations.write_text(json.dumps(a))
    result=run(source,annotations,tmp_path/'result')
    assert result['results'][0]['colour']=={'status':'missing_anatomical_mask'}

def test_coordinate_and_mask_validation():
    with pytest.raises(ValueError):point([1.1,0.5],(10,20))
    with pytest.raises(ValueError):point([float('nan'),0.5],(10,20))
    with pytest.raises(ValueError):polygon_mask([[[0,0],[1,1]]],(10,20))
    assert polygon_mask([], (10,20)) is None

def test_replay_checks_source_archive_before_any_extraction(tmp_path):
    from analysis.v4.prepare_upstream_replay import prepare
    import zipfile
    archive=tmp_path/'fake.zip'
    with zipfile.ZipFile(archive,'w') as z:z.writestr('../escape.txt','no')
    with pytest.raises(ValueError,match='checksum'):
        prepare(archive,tmp_path/'output')
    assert not (tmp_path/'escape.txt').exists()

def test_floral_exclusions_remove_pixels_without_reselecting_colour(tmp_path):
    source,annotations,_=fixture(tmp_path)
    a=json.loads(annotations.read_text());a[0]['floral_exclusions']=[[[.1,.1],[.5,.1],[.5,.9],[.1,.9]]]
    annotations.write_text(json.dumps(a))
    result=run(source,annotations,tmp_path/'result')
    assert 0<result['results'][0]['colour']['n_pixels']<600

def test_anatomical_body_ratio_is_not_spine_count():
    from analysis.v4.anatomical_geometry import involucre
    mask=np.zeros((100,100),bool);mask[10:90,30:70]=True
    x=involucre(mask,(50,90),(50,10),[])
    assert 1.9<x['body_length_width_ratio']<2.2
    assert x['visible_annotated_phyllaries']==0
    assert x['phyllary_length_relative_to_body_width']==[]
    assert 'spine_count' not in x

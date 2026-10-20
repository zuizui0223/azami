"""Synthetic and original-byte schema tests for functional v4 head coordinates."""
import json
import numpy as np
import pytest
from PIL import Image
from analysis.v4.upstream import colour,srgb_pixels_to_lab,digest
from analysis.v4.functional_coordinates import presentation_geometry,phyllary_geometry
from analysis.v4.measure_upstream_anatomy import run

def body():
    m=np.zeros((200,200),bool)
    m[40:160,75:125]=True
    return m,(100.,160.),(100.,40.)

def test_floret_lab_uses_all_pixels_and_exposes_continuous_colour_dispersion():
    rgb=np.empty((20,20,3),dtype=np.uint8)
    rgb[:10]=[245,245,245]
    rgb[10:]=[160,40,170]
    result=colour(rgb,np.ones((20,20),bool))
    assert result['status']=='measured_image_colour'
    assert result['n_pixels']==400
    assert result['low_chroma_fraction']==pytest.approx(.5)
    assert result['chromatic_fraction']==pytest.approx(.5)
    assert result['lab_dispersion_median']>0
    assert result['lab_dispersion_q90']>=result['lab_dispersion_median']
    assert np.allclose(result['lab_joint_median'],
        np.median(srgb_pixels_to_lab(rgb.reshape(-1,3)),axis=0))
    assert result['colour_distribution_calibrated'] is False
    assert set(('lab_a_median','lab_b_median','lab_L_median'))<=set(result)

def test_all_gray_has_no_biological_hue_without_discarding_pixels():
    rgb=np.full((14,12,3),245,np.uint8)
    x=colour(rgb,np.ones((14,12),bool))
    assert x['low_chroma_fraction']==1
    assert x['hue_degrees'] is None
    assert x['lab_L_median']>95 and x['n_pixels']==168

def test_presentation_is_rotation_invariant_to_camera_roll_for_stalk_bend():
    a=presentation_geometry((50,60),(70,30),(50,100))
    rot=lambda p:(-p[1],p[0])
    b=presentation_geometry(*map(rot,[(50,60),(70,30),(50,100)]))
    assert a['head_vs_peduncle_bend_unsigned_deg']==pytest.approx(
        b['head_vs_peduncle_bend_unsigned_deg'])
    assert a['head_vs_peduncle_bend_signed_image_deg']==pytest.approx(
        b['head_vs_peduncle_bend_signed_image_deg'])
    # But the angle versus image vertical varies with camera orientation.
    assert a['image_axis_degrees']!=pytest.approx(b['image_axis_degrees'])
    assert a['gravity_degrees'] is None
    assert a['absolute_gravity_or_nodding_state_inferred'] is False

def test_mirror_flips_image_chirality_but_not_unsigned_bend():
    base=(50,60);top=(70,30);stalk=(50,100)
    x=presentation_geometry(base,top,stalk)
    flip=lambda p:(-p[0],p[1])
    y=presentation_geometry(*map(flip,[base,top,stalk]))
    assert x['head_vs_peduncle_bend_unsigned_deg']==pytest.approx(
        y['head_vs_peduncle_bend_unsigned_deg'])
    assert x['head_vs_peduncle_bend_signed_image_deg']==pytest.approx(
        -y['head_vs_peduncle_bend_signed_image_deg'])

def test_no_phyllary_annotation_is_missing_evidence_not_absent_phyllaries():
    m,base,rim=body()
    q=phyllary_geometry(m,base,rim,())
    assert q['status']=='not_individually_annotated'
    assert q['annotated_visible_phyllaries']==0
    assert q['phyllary_count_estimated'] is False
    assert q['per_whorl']=={}
    assert q['sticky_exudate_measured'] is False

def test_relative_lengths_angles_and_whorls_are_anatomical_only():
    m,base,rim=body()
    phyl=[
       {'base':(80.,80.),'tip':(60.,50.),'whorl':'outer'},
       {'base':(120.,80.),'tip':(140.,50.),'whorl':'outer'},
       {'base':(85.,60.),'tip':(80.,30.),'whorl':'middle'}]
    x=phyllary_geometry(m,base,rim,phyl)
    assert x['status']=='measured_2d_visible_phyllary_axes'
    assert x['annotated_visible_phyllaries']==3
    assert x['per_whorl']['outer']['n_annotated_visible']==2
    assert x['per_whorl']['middle']['n_annotated_visible']==1
    assert x['per_whorl']['outer']['median_outward_displacement_over_body_width']>0
    assert all(z['projected_axis_angle_to_involucre_deg']>0 for z in x['per_phyllary'])
    assert all(not z['botanical_posture_class_inferred'] for z in x['per_phyllary'])
    assert x['true_spine_length_measured'] is False
    assert x['sticky_exudate_measured'] is False

def test_phyllary_axis_angle_is_rotation_and_mirror_invariant_in_2d():
    m,base,rim=body()
    ann=[{'base':(80.,80.),'tip':(60.,50.),'whorl':'outer'}]
    a=phyllary_geometry(m,base,rim,ann)
    mirror=lambda p:(199-p[0],p[1])
    am=phyllary_geometry(np.fliplr(m),mirror(base),mirror(rim),
         [{'base':mirror(ann[0]['base']),'tip':mirror(ann[0]['tip']),'whorl':'outer'}])
    assert a['per_phyllary'][0]['projected_axis_angle_to_involucre_deg']==pytest.approx(
        am['per_phyllary'][0]['projected_axis_angle_to_involucre_deg'])
    assert a['per_phyllary'][0]['projected_outward_displacement_over_body_width']==pytest.approx(
        am['per_phyllary'][0]['projected_outward_displacement_over_body_width'])

def test_unknown_whorl_is_retained_but_fabricated_whorl_is_rejected():
    m,base,rim=body()
    sample={'base':(80.,80.),'tip':(60.,50.)}
    assert 'unknown' in phyllary_geometry(m,base,rim,[sample])['per_whorl']
    with pytest.raises(ValueError,match='whorl'):
        phyllary_geometry(m,base,rim,[{**sample,'whorl':'appressed'}])

def test_record_runner_adds_coordinates_without_promoting_unreviewed_photo(tmp_path):
    source=tmp_path/'source';(source/'originals').mkdir(parents=True)
    photo=source/'originals/1.png'
    Image.new('RGB',(100,100),(210,110,180)).save(photo)
    record={'photo_id':1,'obs_id':1,'status':'downloaded_original',
       'original_path':'originals/1.png','sha256':digest(photo),
       'current_taxon_name':'Cirsium test'}
    (source/'acquisition_ledger.json').write_text(json.dumps([record]))
    a={'photo_id':1,'original_sha256':digest(photo),'review_status':'assistant_visual_proposal',
       'flowering_stage':'anthesis','floral_polygons':[[[.1,.1],[.5,.1],[.5,.4],[.1,.4]]],
       'involucre_polygon':[[.4,.4],[.6,.4],[.6,.9],[.4,.9]],
       'head_base':[.5,.9],'floral_end':[.5,.1],'involucre_rim':[.5,.4],
       'peduncle_point':[.5,.99],
       'anatomical_phyllaries':[{'base':[.45,.6],'tip':[.3,.45],'whorl':'outer'}]}
    annotations=tmp_path/'annotations.json';annotations.write_text(json.dumps([a]))
    result=run(source,annotations,tmp_path/'out')['results'][0]
    q=result['functional_coordinates']
    assert q['visible_floret_colour']['lab_joint_median'] is not None
    assert q['presentation_2d']['status']=='image_projected_presentation'
    assert q['visible_phyllary_axes_2d']['annotated_visible_phyllaries']==1
    assert not result['eligible_for_primary_ecology']
    assert not q['genetic_or_adaptive_function_inferred']
    assert not q['ready_for_primary_ecology']

"""Photo-measurable functional *coordinates* for three capitulum domains.

These are 2-D observed phenotypes. Never infer gravity, real armature density,
mucilage, visitor access efficacy or retained genetic LEGO modules from pixels.
"""
from __future__ import annotations
import numpy as np
from analysis.v4.upstream import angles
from analysis.v4.anatomical_geometry import involucre

WHORLS=frozenset({'outer','middle','inner','unknown'})

def presentation_geometry(head_base,floral_end,peduncle_point=None):
    """Report image vertical and head-to-peduncle bend as DISTINCT quantities.

    Signed bend is screen chirality (flips with mirrored photos), not biologically
    upward/downward orientation. All angles are 2D projections.
    """
    if head_base is None or floral_end is None:
        return {'status':'missing_head_axis','gravity_degrees':None}
    base=np.asarray(head_base,float);top=np.asarray(floral_end,float)
    if base.shape!=(2,) or top.shape!=(2,) or not np.isfinite(base).all() or not np.isfinite(top).all():
        raise ValueError('Head landmark coordinates must be finite xy points')
    h=top-base
    if np.linalg.norm(h)<1e-8:
        return {'status':'degenerate_head_axis','gravity_degrees':None}
    result=angles(base,top,peduncle_point)
    result.update(status='image_projected_presentation',
                  head_vs_peduncle_bend_unsigned_deg=result['peduncle_axis_degrees'],
                  head_vs_peduncle_bend_signed_image_deg=None,
                  signed_image_angle_is_not_gravity=True,
                  absolute_gravity_or_nodding_state_inferred=False)
    if peduncle_point is not None:
        ped=np.asarray(peduncle_point,float)
        if ped.shape!=(2,) or not np.isfinite(ped).all():
            raise ValueError('Peduncle landmark must be a finite xy point')
        stalk=base-ped
        if np.linalg.norm(stalk)>1e-8:
            cross=stalk[0]*h[1]-stalk[1]*h[0]
            result['head_vs_peduncle_bend_signed_image_deg']=float(
                np.degrees(np.arctan2(cross,np.dot(stalk,h))))
    return result


def phyllary_geometry(body_mask,body_base,body_rim,phyllaries=()):
    """Measure annotated individual phyllary axes in an image-aligned frame.

    whorl labels must be explicitly annotated; 'unknown' remains unknown.
    An axis angle is NOT a botanical appressed/ascending/spreading state. It
    is not an estimate of true botanical spine counts or material stickiness.
    """
    attrs={'phyllary_count_estimated':False,'true_spine_length_measured':False,
           'sticky_exudate_measured':False,'access_effect_demonstrated':False}
    if not phyllaries:
        return dict(attrs,status='not_individually_annotated',
                    annotated_visible_phyllaries=0,per_phyllary=[],
                    per_whorl={})
    if body_mask is None or body_base is None or body_rim is None:
        return dict(attrs,status='missing_involucre_body_reference',
                    annotated_visible_phyllaries=len(phyllaries),
                    per_phyllary=[],per_whorl={})
    m=np.asarray(body_mask,dtype=bool)
    if m.ndim!=2:raise ValueError('Expected 2-D involucre body mask')
    body=involucre(m,body_base,body_rim,())
    if body['status']!='measured_2d_anatomical_body':
        return dict(attrs,status=body['status'],annotated_visible_phyllaries=len(phyllaries),
                    per_phyllary=[],per_whorl={})
    base=np.asarray(body_base,float);rim=np.asarray(body_rim,float)
    axis=(rim-base)/np.linalg.norm(rim-base)
    perp=np.array([-axis[1],axis[0]])
    width=body['body_width_pixels'];items=[]
    for a in phyllaries:
        if not isinstance(a,dict):
            raise ValueError('Each phyllary needs explicit base/tip/whorl fields')
        whorl=a.get('whorl','unknown')
        if whorl not in WHORLS:raise ValueError('Unsupported whorl label')
        vbase=np.asarray(a['base'],float);vtip=np.asarray(a['tip'],float)
        if vbase.shape!=(2,) or vtip.shape!=(2,) or not np.isfinite(vbase).all() or not np.isfinite(vtip).all():
            raise ValueError('Phyllary base/tip must be finite xy coordinates')
        v=vtip-vbase;length=np.linalg.norm(v)
        if length<1e-8:raise ValueError('Degenerate annotated phyllary')
        axial_cos=np.clip(np.dot(v,axis)/length,-1,1)
        angle=float(np.degrees(np.arccos(axial_cos)))
        side=float(np.dot(vbase-base,perp))
        # A near-midline insertion provides no reliable image-side label.
        outward=None if abs(side)<.025*width else float(
            np.sign(side)*np.dot(v,perp)/width)
        items.append({'whorl':whorl,
            'projected_free_length_over_body_width':float(length/width),
            'projected_axis_angle_to_involucre_deg':angle,
            'projected_outward_displacement_over_body_width':outward,
            'image_side_identifiable':outward is not None,
            'botanical_posture_class_inferred':False,
            'true_apex_spine_measured':False,
            'local_surface_tangent_measured':False})
    by={}
    for whorl in ('outer','middle','inner','unknown'):
        these=[x for x in items if x['whorl']==whorl]
        if not these:continue
        by[whorl]={'n_annotated_visible':len(these),
                   'median_free_length_over_body_width':float(np.median(
                       [x['projected_free_length_over_body_width'] for x in these])),
                   'median_axis_angle_deg':float(np.median(
                       [x['projected_axis_angle_to_involucre_deg'] for x in these])),
                   'median_outward_displacement_over_body_width':
                       (float(np.median(q)) if (q:=[
                           x['projected_outward_displacement_over_body_width']
                           for x in these
                           if x['projected_outward_displacement_over_body_width'] is not None]) else None)}
    return dict(attrs,status='measured_2d_visible_phyllary_axes',
                annotated_visible_phyllaries=len(items),
                involucre_body_width_pixels=width,per_phyllary=items,per_whorl=by,
                coverage_denominator_total_phyllaries_known=False)

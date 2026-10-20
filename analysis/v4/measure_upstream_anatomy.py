"""Measure direct source-image pixels within explicit anatomical annotations."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw
from analysis.v4.upstream import digest,colour,angles
from analysis.v4.anatomical_geometry import involucre
from analysis.v4.functional_coordinates import presentation_geometry,phyllary_geometry
from analysis.v4.run_upstream_acquisition import decode

def point(p,size):
    if p is None:return None
    if len(p)!=2 or not np.isfinite(p).all() or not all(0<=v<=1 for v in p):raise ValueError('Normalized image coordinates must lie in [0,1]')
    return tuple(float(v)*(n-1) for v,n in zip(p,size))

def polygon_mask(polygons,size):
    if not polygons:return None
    out=Image.new('1',size);draw=ImageDraw.Draw(out)
    for poly in polygons:
        if len(poly)<3:raise ValueError('Polygon needs at least three vertices')
        draw.polygon([point(p,size) for p in poly],fill=1)
    return np.asarray(out,dtype=bool)

def run(source,annotations,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    ledger={r['photo_id']:r for r in json.loads((source/'acquisition_ledger.json').read_text()) if r['status']=='downloaded_original'}
    ann=json.loads(Path(annotations).read_text());results=[];seen=set()
    for a in ann:
        key=(a['photo_id'],a.get('unit_id','head_1'))
        if key in seen:raise ValueError('Duplicate annotated head')
        seen.add(key);r=ledger[a['photo_id']];path=source/r['original_path']
        if digest(path)!=r['sha256'] or a['original_sha256']!=r['sha256']:raise ValueError('Annotation/image byte mismatch')
        im,profile=decode(path);rgb=np.asarray(im)
        floral=polygon_mask(a.get('floral_polygons'),im.size)
        excluded=polygon_mask(a.get('floral_exclusions'),im.size)
        if floral is not None and excluded is not None:floral = floral & ~excluded
        body=polygon_mask([a['involucre_polygon']] if a.get('involucre_polygon') else [],im.size)
        base=point(a.get('head_base'),im.size);top=point(a.get('floral_end'),im.size);stem=point(a.get('peduncle_point'),im.size)
        pairs=[[point(p,im.size) for p in pair] for pair in a.get('anatomical_phyllary_pairs',[])]
        if any(len(pair)!=2 for pair in pairs):raise ValueError('A phyllary annotation needs base and tip')
        new_phyllaries=a.get('anatomical_phyllaries',[])
        if pairs and new_phyllaries:
            raise ValueError('Do not silently duplicate legacy pairs and new phyllary annotations')
        candidates=[{'base':pair[0],'tip':pair[1],'whorl':'unknown'} for pair in pairs]
        candidates.extend({'base':point(row.get('base'),im.size),
                           'tip':point(row.get('tip'),im.size),
                           'whorl':row.get('whorl','unknown')} for row in new_phyllaries)
        body_rim=point(a.get('involucre_rim'),im.size)
        if any(row['base'] is None or row['tip'] is None for row in candidates):
            raise ValueError('Annotated phyllary must have both base and tip')
        colour_result=({'status':'source_colour_unassessable','reason':profile}
                       if profile not in {'ICC_to_sRGB','assumed_sRGB_no_embedded_profile'}
                       else colour(rgb,floral)) if a.get('flowering_stage')=='anthesis' else {
                           'status':'anthesis_not_confirmed'}
        module_colour={'status':colour_result['status'],'lab_joint_median':colour_result.get('lab_joint_median'),
                       'lab_within_head_dispersion_median':colour_result.get('lab_dispersion_median'),
                       'low_chroma_fraction':colour_result.get('low_chroma_fraction'),
                       'camera_colour_not_absolute_reflectance':True}
        result={'photo_id':a['photo_id'],'unit_id':key[1],'obs_id':r['obs_id'],'taxon_name':r.get('current_taxon_name'),
                'source_sha256':r['sha256'],'annotation_sha256':digest(annotations),'review_status':a.get('review_status','unreviewed'),
                'flowering_stage':a.get('flowering_stage','unknown'),'colour_profile':profile,
                'sampling_role':a.get('sampling_role','unknown'),'tissue_scope':a.get('tissue_scope','unspecified'),
                'source_attribution':r.get('attribution'),'source_photo_license':r.get('license'),'source_url':r.get('url'),
                'eligible_for_primary_ecology':False,
                'colour':colour_result,
                'presentation':angles(base,top,stem) if base is not None and top is not None else {'status':'missing_landmarks'},
                'involucre':involucre(body,base,body_rim,pairs),
                'functional_coordinates':{
                    'visible_floret_colour':module_colour,
                    'presentation_2d':presentation_geometry(base,top,stem),
                    'visible_phyllary_axes_2d':phyllary_geometry(body,base,body_rim,candidates),
                    'same_annotated_head':True,
                    'source_image_semantics':'pixel_and_annotated_2d',
                    'genetic_or_adaptive_function_inferred':False,
                    'ready_for_primary_ecology':False}}
        results.append(result)
        overlay=im.copy().convert('RGBA');tint=np.zeros((im.height,im.width,4),dtype=np.uint8)
        if body is not None:tint[body]=[0,150,255,75]
        if floral is not None:tint[floral]=[255,0,180,85]
        overlay=Image.alpha_composite(overlay,Image.fromarray(tint)).convert('RGB');d=ImageDraw.Draw(overlay)
        if base is not None and top is not None:d.line([base,top],fill='yellow',width=max(2,im.width//300))
        if stem is not None and base is not None:d.line([stem,base],fill='cyan',width=max(2,im.width//300))
        overlay.thumbnail((900,900));overlay.save(out/f"{a['photo_id']}_{key[1]}_annotation.png")
    receipt={'measurement_source':'original_image_bytes_and_explicit_anatomical_annotations',
             'n_head_records':len(results),'review_status_counts':{s:sum(r['review_status']==s for r in results) for s in {r['review_status'] for r in results}},
             'v3_phenotype_csv_used':False,'results':results,'primary_ecology_ready':False}
    (out/'anatomical_measurements.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-dir',required=True);p.add_argument('--annotations',required=True);p.add_argument('--out-dir',required=True)
    a=p.parse_args();r=run(a.source_dir,a.annotations,a.out_dir);print(json.dumps({k:v for k,v in r.items() if k!='results'},indent=2))

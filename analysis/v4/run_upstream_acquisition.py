"""Fetch original source bytes from a pre-detector queue; never import v3 traits."""
from __future__ import annotations
import argparse, io, json, time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import requests
from PIL import Image, ImageOps, ImageCms, ImageDraw
from analysis.v4.upstream import LICENSES, original_url, digest

MAX_IMAGE_BYTES=25*1024*1024
MAX_BATCH_BYTES=256*1024*1024

def decode(path):
    with Image.open(path) as src:
        if src.width*src.height>50_000_000:raise ValueError('Image exceeds 50MP safety limit')
        icc=src.info.get('icc_profile'); im=ImageOps.exif_transpose(src).convert('RGB')
        if icc:
            im=ImageCms.profileToProfile(im,ImageCms.ImageCmsProfile(io.BytesIO(icc)),ImageCms.createProfile('sRGB'),outputMode='RGB')
        return im.copy(), ('ICC_to_sRGB' if icc else 'assumed_sRGB_no_embedded_profile')

def run(manifest,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for d in ['originals','previews','metadata']:(out/d).mkdir(exist_ok=True)
    spec=json.loads(Path(manifest).read_text());records=spec['records']
    if not 1<=len(records)<=100:raise ValueError('Diagnostic acquisition cap is 100 images; split larger queues explicitly')
    sess=requests.Session();sess.headers['User-Agent']='azami-v4-upstream-research/0.1 (zuizui0223/azami)'
    def get(url,**kwargs):
        time.sleep(1.1)
        response=sess.get(url,timeout=(15,60),allow_redirects=False,**kwargs)
        response.raise_for_status()
        if response.is_redirect:raise ValueError('Unexpected redirect; no alternate-size fallback')
        return response
    observations={}
    ids=sorted({int(r[0]) for r in records})
    for start in range(0,len(ids),40):
        payload=get('https://api.inaturalist.org/v1/observations',params={'id':','.join(map(str,ids[start:start+40])),'per_page':200}).json()
        for obj in payload.get('results',[]):observations[int(obj['id'])]=obj
    ledger=[];thumbs=[];total=0
    for obs_id,photo_id,taxon_id in records:
        row={'obs_id':obs_id,'photo_id':photo_id,'source_taxon_id':taxon_id,'status':'not_attempted'}
        try:
            obj=observations.get(obs_id)
            if not obj:raise ValueError('observation_unavailable')
            photo=next((p for p in obj.get('photos',[]) if int(p['id'])==photo_id),None)
            if not photo:raise ValueError('photo_no_longer_attached')
            license_code=photo.get('license_code')
            if license_code not in LICENSES:raise ValueError('current_license_outside_public_diagnostic_scope')
            if obj.get('captive'):raise ValueError('currently_captive')
            taxon=obj.get('taxon') or {}
            if taxon.get('rank')!='species':raise ValueError('current_rank_not_species')
            url=original_url(photo['url'])
            row.update(current_taxon_name=taxon.get('name'),current_taxon_id=taxon.get('id'),taxon_changed=taxon.get('id')!=taxon_id,
                       url=url,license=license_code,attribution=photo.get('attribution'),observed_on=obj.get('observed_on'),
                       environment_eligible=not (obj.get('obscured') or obj.get('geoprivacy') in ['private','obscured']) and bool(obj.get('geojson')))
            if row['environment_eligible']:row['public_geojson']=obj.get('geojson')
            (out/'metadata'/f'{photo_id}.json').write_text(json.dumps(row,indent=2)+'\n')
            path=out/'originals'/f'{photo_id}{Path(url).suffix.lower()}'
            downloaded=0
            with get(url,stream=True) as response, path.open('wb') as target:
                for chunk in response.iter_content(128*1024):
                    downloaded+=len(chunk);total+=len(chunk)
                    if downloaded>MAX_IMAGE_BYTES or total>MAX_BATCH_BYTES:raise ValueError('byte_budget_exceeded')
                    target.write(chunk)
            im,profile=decode(path)
            row.update(status='downloaded_original',sha256=digest(path),bytes=downloaded,width=im.width,height=im.height,colour_profile=profile,original_path=str(path.relative_to(out)))
            im.thumbnail((500,500));im.save(out/'previews'/f'{photo_id}.png')
            thumbs.append((photo_id,row['current_taxon_name'],im.copy()))
        except Exception as exc:
            row.update(status='failed_or_ineligible',error=str(exc))
        ledger.append(row)
        (out/'acquisition_ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
        print(json.dumps({k:row.get(k) for k in ['photo_id','status','error']},ensure_ascii=False),flush=True)
    template=[{'photo_id':r['photo_id'],'original_sha256':r['sha256'],'review_status':'unreviewed',
               'flowering_stage':'unknown','head_box':None,'floral_polygons':[], 'involucre_polygon':[],
               'head_base':None,'floral_end':None,'peduncle_point':None,'anatomical_phyllary_pairs':[]}
              for r in ledger if r['status']=='downloaded_original']
    (out/'anatomical_annotations.json').write_text(json.dumps(template,indent=2)+'\n')
    for start in range(0,len(thumbs),16):
        items=thumbs[start:start+16];page=Image.new('RGB',(1200,330*int(np.ceil(len(items)/4))),'white');draw=ImageDraw.Draw(page)
        for k,(pid,name,im) in enumerate(items):
            im.thumbnail((288,285));x=(k%4)*300;y=(k//4)*330;page.paste(im,(x,y));draw.text((x+3,y+290),f'{pid} {name}',fill='black')
        page.save(out/f'review_sheet_{start//16+1}.jpg',quality=90)
    counts=Counter(r['status'] for r in ledger)
    report={'input_manifest_sha256':digest(manifest),'source_metadata_sha256':spec['source_archive_sha256'],
            'queried_at_utc':datetime.now(timezone.utc).isoformat(),'requested_photos':len(records),'requested_taxa':len({r[2] for r in records}),
            'status_counts':dict(counts),'downloaded_bytes':total,'measurement_source':'original_photographs',
            'v3_measurement_inputs_used':False,'anatomical_measurements_completed':0,'anatomical_validation':'pending',
            'bulk_scope':'bounded diagnostic batch, not a replacement population sample','detector_status':'not_run_in_this_acquisition_stage'}
    (out/'source_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    if not counts.get('downloaded_original'):raise RuntimeError('No original image bytes obtained')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--out-dir',required=True)
    a=p.parse_args();print(json.dumps(run(a.manifest,a.out_dir),indent=2))

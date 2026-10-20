"""Verify acquired original bytes and recover geometry for a GRAY-profile photo."""
from __future__ import annotations
import argparse,json,zipfile,shutil
from collections import Counter
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit
from analysis.v4.upstream import digest, original_url
from analysis.v4.run_upstream_acquisition import decode

SOURCE_SHA='5600dbe9db39dddbf0f74a5ba9534d2abfd865be0c3071f781c1c0a4efc7aece'

def prepare(archive,out):
    if digest(archive)!=SOURCE_SHA:raise ValueError('Original-source archive checksum mismatch')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            p=PurePosixPath(info.filename)
            if p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe archive member')
            if not (p.parts[0]=='originals' or info.filename in {'source_receipt.json','acquisition_ledger.json'}):continue
            path=out/p
            if info.is_dir():path.mkdir(parents=True,exist_ok=True);continue
            path.parent.mkdir(parents=True,exist_ok=True)
            with z.open(info) as src,path.open('wb') as dst:shutil.copyfileobj(src,dst)
    rows=json.loads((out/'acquisition_ledger.json').read_text());recovered=[]
    for r in rows:
        is_recovery=r['status']=='failed_or_ineligible' and r.get('error')=='cannot build transform'
        if r['status']!='downloaded_original' and not is_recovery:continue
        if is_recovery:
            r['original_path']='originals/'+str(r['photo_id'])+Path(urlsplit(original_url(r['url'])).path).suffix.lower()
        path=out/r['original_path']
        actual=digest(path)
        if not is_recovery and actual!=r['sha256']:raise ValueError('Individual photo checksum mismatch')
        im,profile=decode(path)
        if is_recovery:
            r['prior_decode_error']=r.pop('error');recovered.append(r['photo_id'])
            r.update(status='downloaded_original',sha256=actual,bytes=path.stat().st_size,width=im.width,height=im.height)
        r['colour_profile']=profile;r['colour_eligible_source']=profile in {'ICC_to_sRGB','assumed_sRGB_no_embedded_profile'}
    (out/'acquisition_ledger.json').write_text(json.dumps(rows,indent=2)+'\n')
    receipt={'original_source_artifact_id':12241906850,'source_archive_sha256':SOURCE_SHA,
             'original_photos':sum(r['status']=='downloaded_original' for r in rows),
             'source_colour_eligible':sum(r.get('colour_eligible_source',False) for r in rows),
             'recovered_geometry_photo_ids':recovered,'profile_counts':dict(Counter(r.get('colour_profile') for r in rows)),
             'live_api_requeried':False,'image_bytes_changed':False,'v3_traits_used':False}
    (out/'replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source-zip',required=True);p.add_argument('--out-dir',required=True)
    a=p.parse_args();print(json.dumps(prepare(a.source_zip,a.out_dir),indent=2))

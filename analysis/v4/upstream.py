"""V4 pre-detector sampling and anatomical-ROI measurements; no v3 trait inputs."""
from __future__ import annotations
import argparse, hashlib, json, re, zipfile
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit
import numpy as np
import pandas as pd

META_SHA = '2f216652602a0b3c42b64780d9b31a3ecb2e7fc1300b28b303e8c369b31884ae'
LICENSES = {'cc0', 'cc-by', 'cc-by-sa'}
HOSTS = {'inaturalist-open-data.s3.amazonaws.com', 'static.inaturalist.org'}
FIELDS = ['obs_id','photo_id','taxon_id','taxon_name','taxon_rank','observed_on',
          'photo_license_code','photo_attribution','captive','obscured','geoprivacy',
          'latitude','longitude','positional_accuracy','raw_image_url']

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def original_url(url):
    u=urlsplit(str(url))
    if u.scheme!='https' or u.hostname not in HOSTS or u.username or u.port:
        raise ValueError('Unsupported photograph host')
    path,n=re.subn(r'/(square|small|medium|large|original)\.(jpg|jpeg|png)$',r'/original.\2',u.path,flags=re.I)
    if n!=1 or not re.fullmatch(r'/photos/\d+/original\.(jpg|jpeg|png)',path,re.I):
        raise ValueError('Unsupported photograph path')
    return urlunsplit(('https',u.netloc,path,'',''))

def truth(s):return s.fillna(False).astype(str).str.lower().isin(['true','1','yes'])

def eligible(f):
    """Legal/species/wild-image frame; lack of coordinates does NOT erase a phenotype."""
    x=f.copy()
    for k in ['obs_id','photo_id','taxon_id']:x[k]=x[k].astype(str)
    keep=x.taxon_rank.eq('species') & x.photo_license_code.fillna('').str.lower().isin(LICENSES) & ~truth(x.captive)
    x=x.loc[keep].copy()
    x['latitude']=pd.to_numeric(x.latitude,errors='coerce');x['longitude']=pd.to_numeric(x.longitude,errors='coerce')
    private=truth(x.obscured)|x.geoprivacy.fillna('').str.lower().isin(['obscured','private'])
    x.loc[private,['latitude','longitude']]=np.nan
    x['environment_eligible']=x.latitude.between(-90,90)&x.longitude.between(-180,180)&~private
    def safe_url(u):
        try:return original_url(u)
        except (ValueError,TypeError):return ''
    x['original_url']=x.raw_image_url.map(safe_url)
    x['download_url_supported']=x.original_url.ne('')
    x['cell']=np.where(x.environment_eligible,
                       np.floor((x.latitude+90)/10).fillna(-1).astype(int).astype(str)+'_'+np.floor((x.longitude+180)/10).fillna(-1).astype(int).astype(str), 'unknown')
    x['selection_hash']=[hashlib.sha256(('20261020|'+p).encode()).hexdigest() for p in x.photo_id]
    return x

def balanced_queue(x,per_taxon=2):
    """Outcome-blind equal taxon cap, cells interleaved, one photo per observation."""
    if per_taxon<1:raise ValueError('per_taxon must be positive')
    a=x.loc[x.download_url_supported].sort_values(['selection_hash','photo_id']).drop_duplicates('obs_id').copy()
    a['cell_round']=a.groupby(['taxon_id','cell']).cumcount()
    a=a.sort_values(['taxon_id','cell_round','selection_hash','photo_id'])
    return a.groupby('taxon_id',sort=True).head(per_taxon).reset_index(drop=True)

def colour(rgb,mask):
    """All pixels in an explicit floral mask. Never select the mask by colour."""
    if mask is None:return {'status':'missing_anatomical_mask'}
    rgb=np.asarray(rgb);mask=np.asarray(mask,dtype=bool)
    if rgb.ndim!=3 or rgb.shape[2]!=3 or mask.shape!=rgb.shape[:2]:raise ValueError('RGB/mask shape mismatch')
    if mask.sum()<10:return {'status':'insufficient_mask_pixels'}
    if rgb.dtype!=np.uint8:raise ValueError('Expected explicitly sRGB-encoded uint8 RGB')
    s=rgb[mask].astype(float)/255.
    linear=np.where(s<=.04045,s/12.92,((s+.055)/1.055)**2.4)
    xyz=linear@np.array([[.4124564,.3575761,.1804375],[.2126729,.7151522,.0721750],[.0193339,.1191920,.9503041]]).T
    xyz=xyz/np.array([.95047,1.,1.08883]);delta=6/29
    f=np.where(xyz>delta**3,np.cbrt(xyz),xyz/(3*delta**2)+4/29)
    lab=np.column_stack([116*f[:,1]-16,500*(f[:,0]-f[:,1]),200*(f[:,1]-f[:,2])])
    pix=rgb[mask];C=np.hypot(lab[:,1],lab[:,2]);chromatic=C>=5
    out={'status':'measured_image_colour','n_pixels':int(len(pix)),
         'rgb_median':np.median(pix,axis=0).tolist(),'low_chroma_fraction':float((~chromatic).mean()),
         'clipped_pixel_fraction':float(((pix<=1)|(pix>=254)).any(axis=1).mean()),
         'hue_degrees':None,'hue_resultant':None,'calibrated_reflectance':False}
    for key,val in [('L',lab[:,0]),('a',lab[:,1]),('b',lab[:,2]),('C',C)]:
        out['lab_'+key+'_median']=float(np.median(val));out['lab_'+key+'_q10_q90']=np.quantile(val,[.1,.9]).tolist()
    if chromatic.sum()>=10:
        h=np.arctan2(lab[chromatic,2],lab[chromatic,1]);v=np.array([np.cos(h).mean(),np.sin(h).mean()])
        out['hue_resultant']=float(np.linalg.norm(v))
        if out['hue_resultant']>=.1:out['hue_degrees']=float(np.degrees(np.arctan2(v[1],v[0]))%360)
    inds=np.argsort(C,kind='stable')[np.linspace(0,len(C)-1,5).astype(int)]
    out['observed_rgb_swatches']=pix[inds].tolist()
    return out

def angles(base,top,peduncle=None):
    def angle(a,b):
        den=np.linalg.norm(a)*np.linalg.norm(b)
        if den<=1e-8:return None
        return float(np.degrees(np.arctan2(abs(a[0]*b[1]-a[1]*b[0]),np.dot(a,b))))
    v=np.asarray(top,float)-np.asarray(base,float)
    return {'image_axis_degrees':angle(v,np.array([0.,-1.])),
            'peduncle_axis_degrees':None if peduncle is None else angle(v,np.asarray(base,float)-np.asarray(peduncle,float)),
            'gravity_degrees':None}

def plan(archive,out,per_taxon=10):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if digest(archive)!=META_SHA:raise ValueError('Metadata archive SHA mismatch')
    pieces=[];source_obs=set();source_rows=0
    with zipfile.ZipFile(archive) as z, z.open('photo_metadata_merged.csv') as stream:
        for f in pd.read_csv(stream,usecols=FIELDS,dtype={'obs_id':str,'photo_id':str,'taxon_id':str},chunksize=100000,low_memory=False):
            source_rows+=len(f);source_obs.update(f.obs_id);pieces.append(eligible(f))
    f=pd.concat(pieces,ignore_index=True).drop_duplicates('photo_id')
    f.to_csv(out/'eligible_photo_frame.csv.gz',index=False,compression='gzip')
    queue=balanced_queue(f,per_taxon);queue.to_csv(out/'balanced_queue.csv',index=False)
    counts=f.groupby(['taxon_id','taxon_name']).agg(photos=('photo_id','size'),observations=('obs_id','nunique')).reset_index()
    counts.to_csv(out/'source_taxon_counts.csv',index=False)
    report={'source_archive_sha256':META_SHA,'source_photo_rows':source_rows,'source_observations':len(source_obs),
            'eligible_photo_rows':len(f),'eligible_observations':int(f.obs_id.nunique()),'eligible_taxa':int(f.taxon_id.nunique()),
            'queue_photos':len(queue),'queue_taxa':int(queue.taxon_id.nunique()),'per_taxon_cap':per_taxon,
            'unsupported_photo_urls':int((~f.download_url_supported).sum()),'environment_eligible_photos':int(f.environment_eligible.sum()),'source_snapshot':'2026-07-03',
            'detector_positive_filter':False,'v3_cohort_filter':False,'flower_colour_filter':False,
            'sampling_rule':'hashed one-photo-per-observation; equal taxon cap; interleave 10-degree cells; retain unknown geography',
            'license_scope':sorted(LICENSES),'current_live_metadata_refresh_required':True}
    (out/'upstream_plan.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--metadata-zip',required=True);p.add_argument('--out-dir',required=True);p.add_argument('--per-taxon',type=int,default=10)
    a=p.parse_args();print(json.dumps(plan(a.metadata_zip,a.out_dir,a.per_taxon),indent=2))

import numpy as np
import pandas as pd
import pytest
from analysis.v4.upstream import original_url, eligible, balanced_queue, colour, angles

def test_original_not_square_or_large_and_host_safe():
    assert original_url('https://inaturalist-open-data.s3.amazonaws.com/photos/1/square.JPG') == 'https://inaturalist-open-data.s3.amazonaws.com/photos/1/original.JPG'
    with pytest.raises(ValueError): original_url('https://example.org/photos/1/square.jpg')

def frame():
    return pd.DataFrame([dict(obs_id=str(i),photo_id=str(i),taxon_id=t,taxon_name=t,taxon_rank='species',photo_license_code='cc-by',captive=False,obscured=False,geoprivacy='',latitude=20+(i%2)*20,longitude=10,raw_image_url=f'https://inaturalist-open-data.s3.amazonaws.com/photos/{i}/square.jpg',photo_attribution='test') for i,t in enumerate(['A']*12+['B']*3+['C']*3)])

def test_balance_independent_of_row_order_and_no_phenotype_selection():
    f=eligible(frame())
    a=balanced_queue(f,per_taxon=2)
    b=balanced_queue(f.sample(frac=1,random_state=9),per_taxon=2)
    assert a.photo_id.tolist()==b.photo_id.tolist()
    assert a.groupby('taxon_id').size().to_dict()=={'A':2,'B':2,'C':2}
    assert a.obs_id.is_unique

def test_eligibility_and_environment_flags_are_separate():
    f=frame();f.loc[0,'photo_license_code']='all-rights-reserved';f.loc[1,'captive']=True
    f.loc[2,'obscured']=True;f.loc[3,'latitude']=np.nan
    x=eligible(f)
    assert not {'0','1'} & set(x.photo_id)
    assert not x.set_index('photo_id').loc[['2','3'],'environment_eligible'].any()
    assert x.set_index('photo_id').loc['2','latitude']!=x.set_index('photo_id').loc['2','latitude']

def test_white_not_dropped_and_direct_pixel_lab():
    img=np.full((20,20,3),[245,245,245],dtype=np.uint8)
    out=colour(img,np.ones((20,20),dtype=bool))
    assert out['n_pixels']==400 and out['lab_L_median']>95
    assert out['lab_C_median']<1 and out['hue_degrees'] is None
    assert out['low_chroma_fraction']==1
    assert colour(img,None)['status']=='missing_anatomical_mask'

def test_colour_mask_not_selected_by_chroma():
    img=np.zeros((20,20,3),np.uint8); img[:10]=[245,245,245];img[10:]=[160,40,170]
    out=colour(img,np.ones((20,20),dtype=bool))
    assert out['n_pixels']==400 and np.isclose(out['low_chroma_fraction'],.5)

def test_image_vs_stem_not_gravity():
    a=angles((10,10),(20,0),(0,20))
    assert a['image_axis_degrees']==pytest.approx(45)
    assert a['peduncle_axis_degrees']==pytest.approx(0)
    assert a['gravity_degrees'] is None
    assert angles((10,10),(20,0),None)['peduncle_axis_degrees'] is None

def test_malformed_legacy_url_is_retained_but_not_downloaded():
    f=frame();f.loc[0,'raw_image_url']='https://inaturalist-open-data.s3.amazonaws.com/photos/0/square.'
    x=eligible(f)
    assert x.set_index('photo_id').loc['0','download_url_supported']==False
    assert '0' not in set(balanced_queue(x,20).photo_id)

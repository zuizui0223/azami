"""Whole-range niche occupancy and trait responses in introduced climates."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.ndimage import gaussian_filter
from statsmodels.stats.multitest import multipletests
from analysis.v4.native_introduced_opportunity import authority_source,angular_cells
from analysis.v4.native_cached_climate import VARIABLES
from analysis.v4.s1_s2_extreme_environment import retrieve,open_raster,sample_grid
from analysis.v5.native_pca_env import (background_mesh,_background_species_audit,_pca_env_2d,
                                         _probability_content,climatic_shift_summary)
from analysis.v5.native_pca_env_fullcoverage import apply_original_five_CHELSA
from analysis.v5.run import cluster_fit,stats_row

BASE=Path('work/submission')
PRIVATE=BASE/'private'
OUT=BASE/'novelty'
TRAITS={'angle':'orientation_angle_degrees_median','chroma':'corolla_lab_chroma_median','aspect':'shape_aspect_ratio_median'}
AXES={'angle':['bio18','bio05','rsds_mean'],'chroma':['rsds_mean','bio05','bio18'],
      'aspect':['bio05','rsds_mean','bio18']}
TEST={'angle':['bio18'],'chroma':['rsds_mean','bio05'],'aspect':['bio05','rsds_mean']}

def wmean(a,w):
    return np.average(a,axis=0,weights=w)

def equalcell(d):
    return (1/d.groupby('_cell').obs_id.transform('size')).to_numpy()

def date_terms(g):
    d=pd.to_datetime(g.observed_on,format='mixed',errors='coerce')
    return pd.DataFrame({'sin':np.sin(d.dt.month*2*np.pi/12),'cos':np.cos(d.dt.month*2*np.pi/12),
                         'year':d.dt.year-2020},index=g.index)

def groups(g,r):
    return g.loc[g.range_state==r].copy()

def traits_comparison(species,intro,coverage,contrasts,slopes,pairs):
    for t,col in TRAITS.items():
        axes=AXES[t]
        d=intro.dropna(subset=[col,'novel',*axes]).copy()
        if t=='angle':d=d[d[col].between(0,180)]
        dt=date_terms(d)
        d=d.loc[dt.notna().all(axis=1)]
        dt=dt.loc[d.index]
        info=dict(species=species,trait=t,coverage=coverage)
        raw=[]
        for k in [0,1]:
            g=d[d.novel==k]; raw.append((len(g),g._cell.nunique(),wmean(g[col],equalcell(g)) if len(g) else np.nan))
        support=dict(n_sim=raw[0][0],cells_sim=raw[0][1],mean_sim=raw[0][2],n_new=raw[1][0],cells_new=raw[1][1],mean_new=raw[1][2])
        if min(q[1] for q in raw)<8:
            contrasts.append(dict(**info,**support,status='insufficient_cells'));continue
        w=1/d.groupby(['novel','_cell']).obs_id.transform('size')
        z=d.novel.to_numpy(float)
        X=np.column_stack([np.ones(len(d)),z,dt])
        fit=cluster_fit(d[col],X,w,d._cell)
        contrasts.append(dict(**info,**support,status='estimated',**stats_row(fit,1)))
        vals=d[axes].to_numpy()
        means=wmean(vals,w);sd=np.sqrt(wmean((vals-means)**2,w))
        clim=(vals-means)/sd
        X=np.column_stack([np.ones(len(d)),z,clim,z[:,None]*clim,dt])
        fit=cluster_fit(d[col],X,w,d._cell)
        for a in TEST[t]:
            j=axes.index(a); idx=2+j; contrast=np.zeros(X.shape[1]);contrast[idx]=contrast[idx+len(axes)]=1
            test=fit.t_test(contrast);lo,hi=np.asarray(test.conf_int()).ravel()
            slopes.append(dict(**info,axis=a,climate_sd=sd[j],slope_sim=fit.params.iloc[idx],sim_lower=fit.conf_int().iloc[idx,0],sim_upper=fit.conf_int().iloc[idx,1],p_sim=fit.pvalues.iloc[idx],slope_new=float(np.asarray(test.effect).item()),
                               new_lower=lo,new_upper=hi,p_new=float(test.pvalue),**{'slope_difference':fit.params.iloc[idx+len(axes)],
                               'p_difference':fit.pvalues.iloc[idx+len(axes)],'difference_lower':fit.conf_int().iloc[idx+len(axes),0],'difference_upper':fit.conf_int().iloc[idx+len(axes),1]}))
    if coverage!=.95:return
    d=intro.dropna(subset=[*TRAITS.values(),'novel','bio18','bio05','rsds_mean']).copy()
    d=d[d[TRAITS['angle']].between(0,180)]
    dt=date_terms(d);d=d[dt.notna().all(axis=1)]
    for k,g in d.groupby('novel'):
        if g._cell.nunique()<8:continue
        w=equalcell(g)
        def z(col):
            arr=g[col].to_numpy();m=wmean(arr,w);return (arr-m)/np.sqrt(wmean((arr-m)**2,w))
        for src,dst in [('angle','aspect'),('angle','chroma'),('aspect','chroma')]:
            for adjust in [False,True]:
                x=np.column_stack([np.ones(len(g)),z(TRAITS[src])])
                if adjust:x=np.column_stack([x,g[['bio18','bio05','rsds_mean']],dt.loc[g.index]])
                f=cluster_fit(z(TRAITS[dst]),x,w,g._cell)
                pairs.append(dict(species=species,novel=int(k),predictor=src,outcome=dst,climate_date_adjusted=adjust,
                                  n=len(g),cells=g._cell.nunique(),**stats_row(f,1)))

def main():
    OUT.mkdir(exist_ok=True)
    auth=authority_source((PRIVATE/'authority350.zip').read_bytes())
    wild=pd.read_csv(PRIVATE/'focal_wild.csv.gz',dtype={'obs_id':str},low_memory=False)
    assert len(wild)==297622
    wild['_cell']=angular_cells(wild.latitude,wild.longitude)
    wild=wild.sort_values('obs_id').drop_duplicates(['taxon_name','range_state','_cell','user_id']).copy()
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path('analysis/v4/environment_source_metadata.json'),work=PRIVATE/'rasters')
    for var in ['chelsa_bio18','chelsa_bio05']:
        p,_=retrieve(var,PRIVATE/'rasters');t=open_raster(p,var,PRIVATE/'rasters');wild[var]=sample_grid(t,var,wild)[0]
        t.unlink(missing_ok=True);p.unlink(missing_ok=True)
    measured=pd.read_csv(PRIVATE/'original-measured.csv.gz',usecols=['obs_id',*TRAITS.values()],dtype={'obs_id':str})
    wild=wild.drop(columns=list(TRAITS.values()),errors='ignore').merge(measured,on='obs_id',how='left',validate='one_to_one')
    frozen=json.loads(Path('analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json').read_text())
    occ=[];breadth=[];contrasts=[];slopes=[];pairs=[];head_overlap=[];configuration=[];centroids=[]
    for sp,rows in wild.groupby('taxon_name'):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        complete=rows.dropna(subset=list(VARIABLES));sn=groups(complete,'native');si=groups(complete,'introduced')
        args=[env[env._range=='native'],env[env._range=='introduced'],sn,si]
        transform,edges,nat,inv=_pca_env_2d(*args,return_grid=True)
        D=1-abs(nat-inv).sum()/2
        assert abs(D-frozen['species'][sp]['PCA_env_background_corrected_Schoener_D'])<1e-5
        for a,stat in climatic_shift_summary(sn,si).items():
            for group in ['native','introduced']:
                lo,mid,hi=stat[group+'_cell_equal_q05_q50_q95'];breadth.append(dict(species=sp,axis=a,group=group,q05=lo,median=mid,q95=hi,width=hi-lo))
        coords=transform(si)
        bins=[np.searchsorted(e,coords[:,i],side='right')-1 for i,e in enumerate(edges)]
        inside=np.ones(len(si),bool)
        for i,e in enumerate(edges):
            bins[i][coords[:,i]==e[-1]]=len(e)-2
            inside &= (bins[i]>=0)&(bins[i]<len(e)-1)
        d=si.loc[inside].rename(columns=lambda c:c.removeprefix('chelsa_')).copy()
        for level in [.9,.95,.99]:
            mask=_probability_content(nat,coverage=level)
            novel=~mask[bins[0][inside],bins[1][inside]]
            frac=inv[~mask].sum()
            if level==.95:assert abs(frac-frozen['species'][sp]['introduced_density_outside_native_95percent_probability_content'])<1e-5
            occ.append(dict(species=sp,coverage=level,D=D,outside_density=frac,n_native=len(sn),n_intro=len(si),outside_grid=int((~inside).sum()),n_novel=int(novel.sum())))
            d['novel']=novel.astype(int)
            from analysis.v5.novel_head_overlap import analysis
            head_overlap.extend(analysis(sp,sn,d,level))
            if level==.95:
                from analysis.v5.config_overlap import run as config_run
                config_rows, vectors = config_run(sp,sn,d)
                configuration.extend(config_rows)
                centroids.extend(vectors)
                from analysis.v5.rain_transfer import run_transfer
                run_transfer(sp,sn,d,OUT)
                from analysis.v5.paired_prediction import run as predict_pairs
                predict_pairs(sp,sn,d,OUT)
            traits_comparison(sp,d,level,contrasts,slopes,pairs)
    vector_table = pd.DataFrame(centroids)
    vector_table["holm_six"] = multipletests(vector_table.p_block_flip, method="holm")[1]
    vector_table.to_csv(OUT/"matched_centroids.csv", index=False)
    config=pd.DataFrame(configuration)
    primary=config.within_block & config.dimensions.eq(3)
    config.loc[primary,"holm_two_species"]=multipletests(config.loc[primary,"p"],method="holm")[1]
    config.to_csv(OUT/"configuration_overlap.csv",index=False)
    head=pd.DataFrame(head_overlap)
    primary=head.design.eq('geographic_pairs') & head.native95_mass.eq(.95) & head.n_traits.eq(3) & head.bandwidth.eq(.5)
    assert primary.sum()==2
    head.loc[primary,'holm_primary']=multipletests(head.loc[primary,'p_swap'],method='holm')[1]
    head.to_csv(OUT/'novel_head_overlap.csv',index=False)
    for name,rs in [('niche',occ),('breadth',breadth),('trait_contrasts',contrasts),('trait_gradients',slopes),('trait_covariation',pairs)]:
        df=pd.DataFrame(rs)
        if name in ['trait_contrasts','trait_gradients','trait_covariation']:
            key='climate_date_adjusted' if name=='trait_covariation' else 'coverage'
            for pv in ['p','p_sim','p_new','p_difference']:
                if pv in df:
                    df['holm_'+pv]=df.groupby(key)[pv].transform(lambda s:s.where(s.isna(),pd.Series(multipletests(s.dropna(),method='holm')[1],index=s.dropna().index)) if s.notna().any() else s)
        df.to_csv(OUT/(name+'.csv'),index=False)
if __name__=='__main__':main()

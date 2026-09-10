"""Show frozen shared-space KDE regions, not just scalar volume summaries."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from run_distribution_breadth_pilot import ENV, rng_for, log_kde


def region(sample, h, seed):
    r=rng_for(seed,'region')
    draws=sample[r.integers(len(sample),size=12000)]+r.normal(size=(12000,3))*h
    cutoff=float(np.quantile(log_kde(draws,sample,h),.1))
    check=sample[r.integers(len(sample),size=12000)]+r.normal(size=(12000,3))*h
    mass=float(np.mean(log_kde(check,sample,h)>=cutoff))
    return cutoff,mass


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--traits',type=Path,required=True);p.add_argument('--environment',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    base=Path('analysis_outputs/pr92_distribution_breadth_pilot_20260910')
    contract=json.loads((base/'contract.json').read_text())
    for k,path in [('traits',a.traits),('environment',a.environment)]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==contract['sha256'][k]
    spec={'question':'Where do taxon distributions occupy and overlap shared trait space?',
          'renderer':'matplotlib 3D voxel boundary and 2D projected silhouettes',
          'selection':'two taxa with most eligible observations per module, before inspecting geometry',
          'sample':'one deterministic 50-observation illustration per taxon, not median of 30 volumes',
          'region':'90 percent KDE highest-density region in frozen shared 3D coordinates; 12000 Monte Carlo calibration and independent check draws',
          'grid':'40 cells per coordinate; common bounds within each module; approximation, not convex hull',
          'palette':['#2878A0','#B57C18'],'distinction':'circle/triangle; solid/dashed projected boundary',
          'limits':'image phenotype, not anatomical or functional support; involucre projection retains only 63.6 percent; orientation stays 1D in previous density figure'}
    (a.out/'contract.json').write_text(json.dumps(spec,indent=2)+'\n')
    t=pd.read_csv(a.traits,dtype={'obs_id':str},usecols=['obs_id','taxon_name','endpoint_id','value','measurement_available'])
    t=t[t.measurement_available.astype(str).str.lower().isin(['true','1'])&np.isfinite(t.value)]
    w=t.pivot(index=['obs_id','taxon_name'],columns='endpoint_id',values='value').reset_index()
    e=pd.read_csv(a.environment,dtype={'obs_id':str}).replace([np.inf,-np.inf],np.nan).dropna(subset=[*ENV,'latitude','longitude'])
    bases={b['group']:b for b in json.loads((base/'shared_bases.json').read_text())}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(16,10));records=[]
    overview,oa=plt.subplots(2,3,figsize=(17,10),gridspec_kw={'width_ratios':[1,1.15,1]})
    support=[];loading_rows=[]
    for col,g in enumerate(['colour','outline','involucre']):
        b=bases[g];m=b['members'];h=b['bandwidth']
        f=w[['obs_id','taxon_name',*m]].dropna().merge(e[['obs_id','taxon_name']],on=['obs_id','taxon_name'],validate='one_to_one')
        counts=f.groupby('taxon_name').size().sort_values(ascending=False,kind='stable');names=list(counts[counts>=50].index[:2])
        samples=[]
        for name in names:
            raw=f[f.taxon_name==name][m].to_numpy();x=((raw-b['centre'])/b['scale'])@np.array(b['loadings'])
            r=rng_for(20260910,g,name,'space_illustration');samples.append(x[r.choice(len(x),50,replace=False)])
        if col<2:
            scatter=oa[col,0]
            for taxon in sorted(counts[counts>=50].index):
                raw=f[f.taxon_name==taxon][m].to_numpy()
                x=((raw-b['centre'])/b['scale'])@np.array(b['loadings'])
                r=rng_for(20260910,g,taxon,'space_illustration');selected=x[r.choice(len(x),50,replace=False)]
                support.append(dict(group=g,taxon=taxon,eligible_observations=int(counts[taxon]),displayed_observations=50))
                if taxon not in names:scatter.scatter(selected[:,0],selected[:,1],s=5,c='0.6',alpha=.3,rasterized=True)
            for j,s in enumerate(samples):
                scatter.scatter(s[:,0],s[:,1],s=16,c=['#2878A0','#B57C18'][j],marker=['o','^'][j],alpha=.85,label=names[j].replace('Cirsium','C.'))
            scatter.set_title(f'{g.capitalize()}: continuous observation scores\n{int((counts>=50).sum())} taxa x 50 observations',fontsize=12)
            scatter.set(xlabel='Shared PC1',ylabel='Shared PC2');scatter.grid(alpha=.12);scatter.legend(fontsize=9)
            labels={'corolla_lab_lightness':'Lightness','corolla_lab_chroma':'Chroma','corolla_hue_sin':'Hue sine','corolla_hue_cos':'Hue cosine','capitulum_outline_aspect_ratio':'Aspect ratio','capitulum_outline_circularity':'Circularity','capitulum_outline_solidity':'Solidity','capitulum_width_profile_cv':'Width-profile variation'}
            la=oa[col,1];v=np.array(b['loadings'])
            from matplotlib.colors import LinearSegmentedColormap
            cmap=LinearSegmentedColormap.from_list('signed',['#B57C18','white','#2878A0'])
            la.imshow(v,vmin=-1,vmax=1,cmap=cmap,aspect='auto')
            la.set_xticks(range(3),['PC1','PC2','PC3']);la.set_yticks(range(len(m)),[labels[z] for z in m])
            la.set_title(f'Which measurements define the axes?\n3 axes retain {b["retained_variance"]:.1%} of variance',fontsize=12)
            for k,member in enumerate(m):
                for q in range(3):
                    la.text(q,k,f'{v[k,q]:+.2f}',ha='center',va='center',fontsize=11)
                    loading_rows.append(dict(group=g,endpoint=member,axis=q+1,coefficient=float(v[k,q])))
            la.set_xlabel('PCA coefficients: gold negative / blue positive\nAxis sign is arbitrary; hue uses sine and cosine',fontsize=9)
        combined=np.concatenate(samples);lo=combined.min(axis=0)-4*h;hi=combined.max(axis=0)+4*h
        edges=[np.linspace(lo[k],hi[k],41) for k in range(3)];mid=[(v[:-1]+v[1:])/2 for v in edges]
        mesh=np.meshgrid(*mid,indexing='ij');points=np.column_stack([v.ravel() for v in mesh])
        corners=np.meshgrid(*edges,indexing='ij')
        ax=fig.add_subplot(2,3,col+1,projection='3d');proj=fig.add_subplot(2,3,col+4);handles=[]
        for j,(name,s) in enumerate(zip(names,samples)):
            colour=['#2878A0','#B57C18'][j];cut,mass=region(s,h,f'{g}|{name}')
            mask=(log_kde(points,s,h)>=cut).reshape((40,40,40))
            assert mask.any() and .87<mass<.93
            assert not any(np.take(mask,side,axis=k).any() for k in range(3) for side in [0,-1])
            ax.voxels(*corners,mask,facecolors=matplotlib.colors.to_rgba(colour,.16),edgecolors=None,shade=False)
            ax.scatter(*s.T,color=colour,s=7,marker=['o','^'][j],alpha=.65)
            silhouette=mask.any(axis=2)
            proj.contourf(mid[0],mid[1],silhouette.T,levels=[.5,1.5],colors=[colour],alpha=.18)
            proj.contour(mid[0],mid[1],silhouette.T,levels=[.5],colors=[colour],linestyles=['solid','dashed'][j],linewidths=1.4)
            proj.scatter(s[:,0],s[:,1],s=10,c=colour,marker=['o','^'][j],alpha=.65)
            if col<2:
                hv=oa[col,2]
                hv.contourf(mid[0],mid[1],silhouette.T,levels=[.5,1.5],colors=[colour],alpha=.18)
                hv.contour(mid[0],mid[1],silhouette.T,levels=[.5],colors=[colour],linestyles=['solid','dashed'][j],linewidths=1.4)
                hv.scatter(s[:,0],s[:,1],s=10,c=colour,marker=['o','^'][j],alpha=.65)
            handles.append(Patch(facecolor=matplotlib.colors.to_rgba(colour,.3),edgecolor=colour,label=name.replace('Cirsium','C.')))
            records.append(dict(group=g,taxon=name,eligible_observations=int(counts[name]),sample_size=50,log_density_cutoff=cut,independent_mass=mass,retained_variance=b['retained_variance']))
        ax.set_title(f'{g.capitalize()} | retained variance {b["retained_variance"]:.1%}',pad=15)
        ax.set_xlabel('Shared PC1');ax.set_ylabel('Shared PC2');ax.set_zlabel('Shared PC3')
        ax.set_box_aspect(hi-lo);ax.view_init(elev=22,azim=-55)
        ax.legend(handles=handles,loc='upper left',fontsize=9)
        proj.set(xlabel='Shared PC1',ylabel='Shared PC2',xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]))
        proj.set_aspect('equal',adjustable='box');proj.set_title('Projection of the same 3D regions',fontsize=11);proj.grid(alpha=.15)
        print(g,names,flush=True)
        if col<2:
            hv.set(xlabel='Shared PC1',ylabel='Shared PC2')
            hv.set_title('Taxon-specific hypervolumes\nProjection of 90% regions in 3D',fontsize=12)
            hv.legend(handles=handles,fontsize=9);hv.grid(alpha=.12)
            # Same displayed coordinates in the first and third columns, including all plotted points.
            xl=[scatter.get_xlim(),hv.get_xlim()];yl=[scatter.get_ylim(),hv.get_ylim()]
            for target in [scatter,hv]:
                target.set_xlim(min(z[0] for z in xl),max(z[1] for z in xl));target.set_ylim(min(z[0] for z in yl),max(z[1] for z in yl))
    overview.suptitle('Continuous image measurements -> shared PCA -> trait-space regions',fontsize=19,y=.98)
    overview.text(.5,.035,'All eligible taxa contribute equally to the displayed point cloud. Two most-sampled taxa illustrate the regions.\nPCA uses the frozen equal-taxon basis; regions use one 50-observation sample. No categorical benchmark or biological-accuracy claim.',ha='center',fontsize=10)
    overview.subplots_adjust(left=.055,right=.98,top=.86,bottom=.15,wspace=.60,hspace=.55)
    overview.savefig(a.out/'continuous_trait_workflow.png',dpi=170);overview.savefig(a.out/'continuous_trait_workflow.svg');plt.close(overview)
    pd.DataFrame(support).to_csv(a.out/'display_support.csv',index=False)
    pd.DataFrame(loading_rows).to_csv(a.out/'axis_coefficients.csv',index=False)
    fig.suptitle('Capitulum trait-space hypervolumes\n90% KDE regions; 50 observations per taxon; shared axes within each module',fontsize=18,y=.98)
    fig.text(.5,.025,'Illustrative single resample, not the 30-resample median. Lower panels are silhouettes, not 2D 90% regions.\nInvolucre is an incomplete 3D projection. Smoothed image phenotypes include photographic and measurement variation.',ha='center',fontsize=10)
    fig.subplots_adjust(left=.045,right=.96,bottom=.12,top=.86,wspace=.25,hspace=.25)
    fig.savefig(a.out/'hypervolume_space.png',dpi=170);fig.savefig(a.out/'hypervolume_space.svg');plt.close(fig)
    pd.DataFrame(records).to_csv(a.out/'region_checks.csv',index=False)


if __name__=='__main__':main()

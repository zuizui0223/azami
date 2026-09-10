"""Bounded calibration and static figures for PR92 distribution pilots."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_distribution_breadth_pilot import ENV, log_kde, rng_for
from run_distribution_overlap_pilot import overlap, wasserstein_parts

CORE=['orientation','colour','outline']
BLUE='#2878A0'; GOLD='#B57C18'; INK='#252525'


def calibrate(samples,group,bandwidth,pairs,seed):
    baselines=[]
    for taxon,x in samples.items():
        if len(x)<100:continue
        r=rng_for(seed,group,taxon,'disjoint_baseline')
        values=[]
        for rep in range(30):
            take=r.choice(len(x),100,replace=False)
            a,b=x[take[:50]],x[take[50:]]
            values.append(wasserstein_parts(a,b))
        v=np.array(values)
        baselines.append(dict(group=group,taxon_name=taxon,observations=len(x),
                              median_same_taxon_w2=float(np.median(v[:,0])),
                              same_taxon_w2_p95=float(np.quantile(v[:,0],.95)),
                              centred_same_taxon_w2_p95=float(np.quantile(v[:,2],.95))))
    lookup={row['taxon_name']:row for row in baselines}
    comparisons=[]
    for row in pairs.itertuples():
        if row.taxon_a not in lookup or row.taxon_b not in lookup:continue
        reference=max(lookup[row.taxon_a]['same_taxon_w2_p95'],lookup[row.taxon_b]['same_taxon_w2_p95'])
        cref=max(lookup[row.taxon_a]['centred_same_taxon_w2_p95'],lookup[row.taxon_b]['centred_same_taxon_w2_p95'])
        comparisons.append(dict(group=group,taxon_a=row.taxon_a,taxon_b=row.taxon_b,
                                between_w2=row.wasserstein_w2,within_reference_p95=reference,
                                exceeds_reference=bool(row.wasserstein_w2>reference),
                                centred_between_w2=row.centred_w2,centred_within_reference_p95=cref,
                                centred_exceeds_reference=bool(row.centred_w2>cref)))
    # Deterministic outcome-blind subset, not hand-picked stable pairs.
    r=rng_for(seed,group,'bandwidth_pair_subset')
    chosen=pairs.iloc[np.sort(r.choice(len(pairs),min(100,len(pairs)),replace=False))]
    bandwidth_rows=[]
    for row in chosen.itertuples():
        values={m:[] for m in [.75,1.,1.25]}
        r=rng_for(seed,group,row.taxon_a,row.taxon_b,'bandwidth')
        for rep in range(10):
            x=samples[row.taxon_a];y=samples[row.taxon_b]
            a=x[r.choice(len(x),50,replace=False)];b=y[r.choice(len(y),50,replace=False)]
            ia=r.integers(50,size=128);ib=r.integers(50,size=128)
            za=r.normal(size=(128,a.shape[1]));zb=r.normal(size=(128,b.shape[1]))
            for factor in values:
                h=bandwidth*factor;da=a[ia]+za*h;db=b[ib]+zb*h
                values[factor].append(overlap(a,b,h,da,db,log_kde(da,a,h),log_kde(db,b,h)))
        bandwidth_rows.extend(dict(group=group,taxon_a=row.taxon_a,taxon_b=row.taxon_b,factor=f,
                                   density_overlap=float(np.median(v))) for f,v in values.items())
    return baselines,comparisons,bandwidth_rows


def finish_figure(fig,path):
    fig.savefig(path.with_suffix('.png'),dpi=180,facecolor='white')
    fig.savefig(path.with_suffix('.svg'),facecolor='white')
    plt.close(fig)


def scatter_figure(breadths,associations,out):
    fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained')
    for j,group in enumerate(CORE):
        f=breadths[breadths.group==group]
        for i,predictor in enumerate(['environment_breadth','geographic_breadth_km']):
            ax=axes[i,j]
            ax.scatter(f[predictor],f.kde_volume,s=26,facecolor=BLUE,edgecolor='white',linewidth=.4,alpha=.8)
            for taxon,offset in [('Cirsium arvense',(-8,16)),('Cirsium vulgare',(-8,-18))]:
                one=f[f.taxon_name==taxon]
                if len(one):
                    ax.annotate(taxon.replace('Cirsium','C.'),(one[predictor].iloc[0],one.kde_volume.iloc[0]),xytext=offset,textcoords='offset points',fontsize=8,ha='right',arrowprops={'arrowstyle':'-','color':'0.4','lw':.6})
            r=associations[(associations.group==group)&(associations.predictor==predictor)].iloc[0]
            ax.set_title(f'{group.capitalize()} | rho = {r.rho:+.2f}',loc='left',fontsize=11)
            ax.set_xlabel('Sampled environmental breadth (standardized RMS)' if i==0 else 'Sampled geographic breadth (RMS km)',fontsize=9)
            ax.set_ylabel('90% KDE volume in shared coordinates',fontsize=9)
            ax.grid(alpha=.15);ax.set_axisbelow(True);ax.margins(.12)
    fig.suptitle('Head-phenotype breadth and sampled range breadth\n54 taxa per module; 50 observations per taxon, median of 30 rarefactions',fontsize=13)
    fig.supxlabel('Exploratory: no spatial/phylogenetic or multiple-comparison correction. Volume scales differ among modules.',fontsize=9)
    finish_figure(fig,out/'breadth_scatter')


def distribution_figure(sample_sets,bases,pairs,out):
    fig,axes=plt.subplots(3,2,figsize=(13,13),layout='constrained')
    references=pd.read_csv(out/'between_vs_same_taxon.csv')
    selections=[]
    for i,group in enumerate(CORE):
        f=pairs[pairs.group==group].copy()
        # Common support rule before choosing diagnostic examples.
        eligible={t for t,x in sample_sets[group].items() if len(x)>=100}
        f=f[f.taxon_a.isin(eligible)&f.taxon_b.isin(eligible)]
        near=f[f.centroid_distance<=f.centroid_distance.quantile(.25)]
        examples=[near.sort_values(['centred_w2','taxon_a','taxon_b'],ascending=[False,True,True]).iloc[0],
                  f.sort_values(['density_overlap','taxon_a','taxon_b'],ascending=[False,True,True]).iloc[0]]
        for j,row in enumerate(examples):
            ax=axes[i,j];rule='Nearby means; largest centred distance' if j==0 else 'Largest estimated density overlap'
            ref=references[(references.group==group)&(references.taxon_a==row.taxon_a)&(references.taxon_b==row.taxon_b)].iloc[0]
            above=bool(ref.centred_exceeds_reference)
            selections.append(dict(group=group,selection_rule=rule,taxon_a=row.taxon_a,taxon_b=row.taxon_b,
                                   density_overlap=row.density_overlap,centroid_distance=row.centroid_distance,centred_w2=row.centred_w2,
                                   centred_within_reference_p95=ref.centred_within_reference_p95,centred_exceeds_reference=above))
            for taxon,colour,marker,style in [(row.taxon_a,BLUE,'o','-'),(row.taxon_b,GOLD,'^','--')]:
                raw=sample_sets[group][taxon]
                r=rng_for(20260910,group,taxon,'figure');x=raw[r.choice(len(raw),50,replace=False)]
                label=taxon.replace('Cirsium','C.')
                if group=='orientation':
                    combined=np.concatenate([sample_sets[group][row.taxon_a],sample_sets[group][row.taxon_b]])
                    h=bases[group]['bandwidth'];grid=np.linspace(combined.min()-2*h,combined.max()+2*h,300)[:,None]
                    ax.plot(grid[:,0],np.exp(log_kde(grid,x,h)),color=colour,ls=style,label=label,lw=1.8)
                    ax.set_xlabel('Shared standardized angle coordinate');ax.set_ylabel('KDE density')
                else:
                    ax.scatter(x[:,0],x[:,1],s=23,color=colour,marker=marker,alpha=.55,label=label)
                    ax.scatter(x[:,0].mean(),x[:,1].mean(),s=140,color=colour,marker='X',edgecolor=INK,linewidth=.6)
                    ax.set_xlabel('Shared PC1');ax.set_ylabel('Shared PC2')
            ax.set_title(f'{group.capitalize()} | {rule}\nFull-coordinate overlap = {row.density_overlap:.2f}\nCentred distance above same-taxon reference: {"yes" if above else "no"}',loc='left',fontsize=10)
            ax.legend(loc='best',fontsize=8,framealpha=.9);ax.grid(alpha=.15);ax.set_axisbelow(True)
        xlimits=[ax.get_xlim() for ax in axes[i]];ylimits=[ax.get_ylim() for ax in axes[i]]
        for ax in axes[i]:
            ax.set_xlim(min(x[0] for x in xlimits),max(x[1] for x in xlimits))
            ax.set_ylim(min(y[0] for y in ylimits),max(y[1] for y in ylimits))
    fig.suptitle('Examples of taxon phenotype distributions\nRule-selected illustrations; 50 observations per taxon; not independent tests',fontsize=14)
    fig.supxlabel('Colour/outline panels show 2 of 3 coordinates; overlap uses all coordinates and 30 rarefactions. X marks plotted-sample means.',fontsize=9)
    finish_figure(fig,out/'distribution_examples')
    pd.DataFrame(selections).to_csv(out/'example_selection.csv',index=False)


def main():
    p=argparse.ArgumentParser();p.add_argument('--traits',type=Path,required=True);p.add_argument('--environment',type=Path,required=True)
    p.add_argument('--breadth-dir',type=Path,required=True);p.add_argument('--overlap-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    contract=json.loads((args.breadth_dir/'contract.json').read_text())
    for key,path in [('traits',args.traits),('environment',args.environment)]:assert hashlib.sha256(path.read_bytes()).hexdigest()==contract['sha256'][key]
    chart_contract=dict(renderer='matplotlib static PNG/SVG scientific figures',
                        scatter=dict(grain='taxon per module',rows=54,panels=6,takeaway='show every eligible taxon without fitted causal trend',family='scatter'),
                        distributions=dict(family='density lines and projected scatter',panels=6,selection='>=100 observations; bottom-quartile centroid distance then maximum centred W2, and maximum overlap',
                                           caveat='diagnostic outcome-selected illustrations, not validation'),
                        palette=dict(blue=BLUE,gold=GOLD,neutral=INK,policy='two-root cap',noncolour='circle versus triangle, solid versus dashed'),
                        baseline='30 disjoint 50+50 samples per taxon; compare prior between-taxon median W2 to larger taxon-specific same-taxon p95; not a hypothesis test',
                        bandwidth='100 outcome-blind seeded random pairs per core module; 10 matched repetitions at 0.75/1/1.25 bandwidth; 128 draws per taxon')
    (args.out/'chart_and_analysis_contract.json').write_text(json.dumps(chart_contract,indent=2)+'\n')
    t=pd.read_csv(args.traits,usecols=['obs_id','taxon_name','endpoint_id','value','measurement_available'],dtype={'obs_id':str})
    t=t[t.measurement_available.astype(str).str.lower().isin(['true','1'])&np.isfinite(t.value)]
    wide=t.pivot(index=['obs_id','taxon_name'],columns='endpoint_id',values='value').reset_index()
    e=pd.read_csv(args.environment,dtype={'obs_id':str}).replace([np.inf,-np.inf],np.nan).dropna(subset=[*ENV,'latitude','longitude'])
    bases={b['group']:b for b in json.loads((args.breadth_dir/'shared_bases.json').read_text())}
    pairs=pd.read_csv(args.overlap_dir/'pairwise_distributions.csv')
    allbase=[];allcompare=[];allband=[];sample_sets={}
    for group in CORE:
        b=bases[group];members=b['members']
        f=wide[['obs_id','taxon_name',*members]].dropna().merge(e[['obs_id','taxon_name']],on=['obs_id','taxon_name'],validate='one_to_one').reset_index(drop=True)
        x=((f[members].to_numpy()-b['centre'])/b['scale'])@np.array(b['loadings'])
        samples={taxon:x[np.array(list(idx))] for taxon,idx in f.groupby('taxon_name').groups.items() if len(idx)>=50}
        sample_sets[group]=samples
        a,c,d=calibrate(samples,group,b['bandwidth'],pairs[pairs.group==group],20260910)
        allbase+=a;allcompare+=c;allband+=d
        print(f'{group}: {len(a)} same-taxon baselines; {len(c)} calibrated pairs',flush=True)
    baseline=pd.DataFrame(allbase);compare=pd.DataFrame(allcompare);band=pd.DataFrame(allband)
    baseline.to_csv(args.out/'same_taxon_baselines.csv',index=False);compare.to_csv(args.out/'between_vs_same_taxon.csv',index=False);band.to_csv(args.out/'bandwidth_check.csv',index=False)
    summary=[]
    for group in CORE:
        c=compare[compare.group==group];b=band[band.group==group].pivot(index=['taxon_a','taxon_b'],columns='factor',values='density_overlap')
        summary.append(dict(group=group,baseline_taxa=int((baseline.group==group).sum()),pairs=len(c),
                            between_above_within_p95=int(c.exceeds_reference.sum()),centred_above_within_p95=int(c.centred_exceeds_reference.sum()),
                            bandwidth_rank_rho_075=float(spearmanr(b[1.],b[.75]).statistic),bandwidth_rank_rho_125=float(spearmanr(b[1.],b[1.25]).statistic),
                            median_overlap_by_bandwidth={str(k):float(b[k].median()) for k in b.columns}))
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'text.color':INK,'axes.labelcolor':INK,'svg.fonttype':'none'})
    scatter_figure(pd.read_csv(args.breadth_dir/'taxon_breadths.csv'),pd.read_csv(args.breadth_dir/'breadth_associations.csv'),args.out)
    distribution_figure(sample_sets,bases,pairs,args.out)


if __name__=='__main__':main()

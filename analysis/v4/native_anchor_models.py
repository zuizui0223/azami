"""Exploratory native-only scalar checks; not frozen v3 confirmatory effects."""
import argparse,json
from pathlib import Path
import pandas as pd
import statsmodels.api as sm
from analysis.v4.native_environment import VARS,spherical,weights,taxon_summary,normalize


def run(native_measurements,native_environment,outdir):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    f=pd.read_csv(native_measurements,dtype={'obs_id':str},low_memory=False)
    e=pd.read_csv(native_environment,dtype={'obs_id':str},low_memory=False)
    f=f.merge(e[['obs_id',*VARS]],on='obs_id',validate='one_to_one')
    top=['wcvp:'+str(x) for x in e.groupby('accepted_key').size().sort_values(ascending=False).index]
    out=[]
    for response,predictor in [('corolla_lab_chroma','chelsa_rsds_mean'),('orientation_image_vertical_angle','chelsa_bio12')]:
        for drop in [0,2,5]:
            x=f[~f.taxon_name.isin(top[:drop])].dropna(subset=[response]+VARS).copy()
            n=x.groupby('taxon_name').size();x=x[x.taxon_name.isin(n[n>=5].index)].reset_index(drop=True)
            w=weights(x);s=spherical(x.latitude,x.longitude);cols=['space_'+str(i) for i in range(s.shape[1])];x[cols]=s
            a=taxon_summary(x,[response]+VARS+cols,w);y=normalize(a[[response]])[0][:,0]
            for kind,terms in [('marginal',[predictor]),('joint_environment',VARS),('joint_environment_and_space',VARS+cols)]:
                fit=sm.OLS(y,sm.add_constant(normalize(a[terms])[0])).fit(cov_type='HC3');idx=1+terms.index(predictor)
                out.append({'response':response,'predictor':predictor,'exclude_top':drop,'n_observations':len(x),'n_taxa':len(a),'model':kind,
                            'beta':float(fit.params[idx]),'ci95':fit.conf_int()[idx].tolist(),'p_unadjusted':float(fit.pvalues[idx]),'r2':float(fit.rsquared)})
    (outdir/'native_anchor_models.json').write_text(json.dumps(out,indent=2))
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--measurements',type=Path,required=True);p.add_argument('--environment',type=Path,required=True);p.add_argument('--out-dir',type=Path,required=True)
    a=p.parse_args();run(a.measurements,a.environment,a.out_dir)

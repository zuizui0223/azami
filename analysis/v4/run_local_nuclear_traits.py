"""Execute continuous S1 ancestral/branch analyses on the explicit local panel."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo
from analysis.v4.native_integration import FEATURES
from analysis.v4 import s1s2_multitrait_branch_core as core
from analysis.v4.local_nuclear_trait_analysis import join_moments,rooted_ingroup,evaluate_unit

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True)
    p.add_argument('--moments-dir',type=Path,required=True);p.add_argument('--permutations',type=int,default=9999)
    a=p.parse_args();w=a.workspace;r=w/'work/s1-local-nuclear';out=w/'outputs'
    if a.permutations<9999:raise ValueError('Prespecified analysis requires >=9999 shuffles')
    source=a.moments_dir/'S1_FROZEN_PRIMARY_TAXON_MOMENTS.csv'
    receipt=json.loads((a.moments_dir/'S1_FROZEN_PRIMARY_TAXON_COHORTS_RECEIPT.json').read_text())
    assert digest(source)==receipt['aggregate_moments']['sha256']
    assert digest(a.moments_dir/'S1_FROZEN_PRIMARY_TAXON_COHORTS.csv')==digest(out/'S1_FROZEN_PRIMARY_TAXON_COHORTS.csv')
    elevations=[x for x in receipt['raster_sources'] if x['variable']=='topo_elevation']
    assert len(elevations)==1 and elevations[0]['variable_unit']=='m'
    moments=pd.read_csv(source,dtype={'accepted_key':str})
    moments['coordinate']=moments.coordinate.replace({'topo_elevation':'elevation_m'})
    membership=pd.read_csv(a.moments_dir/'S1_FROZEN_PRIMARY_TAXON_COHORTS.csv',dtype={'accepted_key':str})
    joined=moments.merge(membership[['cohort','accepted_key','n_observations']],on=['cohort','accepted_key'],validate='many_to_one',suffixes=('','_frozen'))
    assert len(joined)==len(moments) and (joined.n_observations==joined.n_observations_frozen).all()
    manifest=pd.read_csv(r/'expanded24_run_manifest.csv',dtype={'accepted_key':str})
    assert len(manifest)==24 and manifest.run.is_unique and manifest.accepted_key.is_unique
    assert set(manifest.accepted_key)<=set(membership.loc[membership.cohort=='common8','accepted_key'])
    trees={'ML':out/'S1_twentyfour_taxon_concatenated_ML.nwk',
           'ASTRAL_ML_lengths':out/'S1_twentyfour_taxon_ASTRAL_topology_ML_lengths.nwk'}
    expected={'ML':'af50447d535c2a6136648fe50e23d5f3012b7be86dbe6381fbbce4ae6b55ab0d',
              'ASTRAL_ML_lengths':'af8e92f166050d178ea08c11f6ef7aca48560b9e9412bdcc0ef2180d30438e92'}
    records=[];tests=[];inputs=[]
    for label,path in trees.items():
        assert digest(path)==expected[label]
        raw=Phylo.read(path,'newick')
        for exclude in [False,True]:
            panel=label+('_23' if exclude else '_24')
            selected=manifest.loc[manifest.run.ne('SRR30887185')] if exclude else manifest
            tree=rooted_ingroup(raw,manifest.run,selected.run)
            for cohort_mode in ['common8','endpoint_specific']:
                for unit,columns in FEATURES.items():
                    cohort='common8' if cohort_mode=='common8' else unit
                    means=join_moments(moments,selected,cohort,[*columns,*core.FOUR])
                    seed=core.stable_seed('|'.join([cohort_mode,unit,*sorted(selected.run)]))
                    result=evaluate_unit(tree,means,columns,nperm=a.permutations,seed=seed)
                    base=dict(panel=panel,cohort_mode=cohort_mode,unit=unit,n_tips=len(means),
                        n_original_observations=int(means.n_original_observations.sum()))
                    records.append({**base,**result})
                    for scope,s in result['scopes'].items():
                        for axis,v in s.get('axes',{}).items():
                            tests.append({**base,'scope':scope,'environment':axis,'edge_count':s['edge_count'],**v})
                    for run,row in means.iterrows():inputs.append({**base,'run':run,**row.to_dict()})
                print(panel,cohort_mode,'complete',flush=True)
                (out/'T2_local24_execution_status.json').write_text(json.dumps(dict(state='RUNNING',last_completed_panel=panel,last_completed_cohort=cohort_mode)))
    table=pd.DataFrame(tests)
    for _,g in table.groupby(['panel','cohort_mode','scope']):
        assert len(g)==32 and len(g[['unit','environment']].drop_duplicates())==32
        adjusted=core.holm_family([(str(i),v) for i,v in zip(g.index,g.conditional_tip_shuffle_fraction)])
        for i in g.index:table.loc[i,'holm_conditional_fraction']=adjusted[str(i)]
    result=dict(status='COMPLETED_CONDITIONAL_CONTINUOUS_ANALYSIS',n_permutations=a.permutations,
        constructs=FEATURES,environments=list(core.FOUR),source_sha256={str(source):digest(source),**{str(p):digest(p) for p in trees.values()}},
        tree_units='substitutions/site; ASTRAL topology independently fitted with sequence likelihood',
        cohort_primary='common8, same photographs for all constructs and environments',
        family='8 constructs x 4 environments, Holm separately per tree/tip-set/cohort/edge scope; alternatives are sensitivities, not independent replication',
        uncertainty_limits=['Fixed observed tips; sampling SD is exported but not propagated','Intervals condition on one tree and plug-in diffusion rate','No historical climate, discrete U/D transitions, adaptation or genetic causation inferred','Conditional unrestricted tip shuffles are not a fitted phylogenetic-process null'],
        records=records)
    (out/'T2_LOCAL24_CONTINUOUS_ANCESTRY_ENVIRONMENT.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    table.to_csv(out/'T2_local24_environment_branch_tests.csv',index=False)
    pd.DataFrame(inputs).to_csv(out/'T2_local24_matched_tip_means.csv',index=False)
    (out/'T2_local24_execution_status.json').write_text(json.dumps(dict(state='COMPLETED',models=len(records),environment_tests=len(table),inference_accepted=False)))
    print('Completed',len(records),'models and',len(table),'environment tests',flush=True)

if __name__=='__main__':main()

"""Observation-level within/between models with joint species-cluster inference.

Each response has a Gaussian random species intercept and the same eight climate
terms. Response-specific missingness is retained. Stacked cluster scores provide
joint colour/projection tests; this is marginal working-covariance inference,
not an unrestricted multivariate random-effects likelihood.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import f, t
from statsmodels.stats.multitest import multipletests

ENV = ['chelsa_bio01','chelsa_rsds_mean','chelsa_vpd_mean','chelsa_bio12']
TERMS = ['intercept']+[v+s for v in ENV for s in ('__between','__within')]
MODULES = {
    'colour': ['corolla_lab_lightness','corolla_lab_chroma','corolla_hue_sin','corolla_hue_cos'],
    'orientation': ['orientation_image_vertical_angle'],
    'bract_projection': ['bract_projection_p95','bract_spread_fraction'],
}

def design_diagnostics(frame):
    """Describe actual endpoint-specific eight-term designs, without changing fits."""
    terms=TERMS[1:]
    def inspect(d, columns):
        a=d[columns].to_numpy(float)
        sd=a.std(axis=0)
        if len(a)<=len(columns) or np.any(sd<=1e-12):
            return {'status':'not_identifiable','n':len(a),'terms':columns}
        z=(a-a.mean(axis=0))/sd
        c=z.T@z/len(z)
        eigen=np.linalg.eigvalsh(c)
        if eigen[0]<=1e-12:
            return {'status':'singular','n':len(a),'terms':columns,
                    'correlation':c.tolist()}
        vif=np.diag(np.linalg.inv(c))
        return {'status':'estimated','n':len(a),'terms':columns,
                'correlation':c.tolist(),'vif':vif.tolist(),
                'remaining_sd_fraction_after_other_terms':(1/np.sqrt(vif)).tolist(),
                'standardized_condition_index':float(np.sqrt(eigen[-1]/eigen[0]))}
    result={}
    for outcomes in MODULES.values():
        for outcome in outcomes:
            mask=(frame[outcome+'__eligible'].eq(True)&np.isfinite(frame[outcome])
                  &np.isfinite(frame[ENV+terms]).all(axis=1))
            d=frame.loc[mask]
            result[outcome]={level:inspect(d,columns) for level,columns in [
                ('full_eight_terms',terms),
                ('between_observation_weighted',[v+'__between' for v in ENV]),
                ('within',[v+'__within' for v in ENV])]}
    return {'scope':'unweighted design diagnostics on exact endpoint model rows; '
                    'not REML precision, causal identification, or species-equal summaries',
            'endpoints':result}


def fit_response(x, y, groups):
    """Exact random-intercept REML using sums of crossproducts, not mean responses."""
    x,y,groups = np.asarray(x,float),np.asarray(y,float),np.asarray(groups)
    n,p=x.shape
    labels,inverse=np.unique(groups,return_inverse=True)
    g=len(labels)
    if g<=p+2 or np.linalg.matrix_rank(x)<p or np.std(y)<=1e-12:
        return {'status':'not_identifiable','n':n,'n_species':g}
    sizes=np.bincount(inverse).astype(float)
    sx=np.zeros((g,p));sy=np.zeros(g)
    np.add.at(sx,inverse,x);np.add.at(sy,inverse,y)
    xx=np.zeros((g,p,p));xy=np.zeros((g,p))
    for j in range(p):
        np.add.at(xy[:,j],inverse,x[:,j]*y)
        for k in range(p):
            np.add.at(xx[:,j,k],inverse,x[:,j]*x[:,k])
    yy=float(y@y)
    def evaluate(ratio, details=False):
        a=ratio/(1+ratio*sizes)
        info=xx.sum(axis=0)-np.einsum('g,gi,gj->ij',a,sx,sx)
        rhs=xy.sum(axis=0)-sx.T@(a*sy)
        beta=np.linalg.solve(info,rhs)
        sse=yy-np.sum(a*sy*sy)-beta@rhs
        sign,logdet=np.linalg.slogdet(info)
        if sign<=0 or sse<=0:
            return np.inf
        objective=(n-p)*np.log(sse/(n-p))+np.log1p(ratio*sizes).sum()+logdet
        if not details:return objective
        inv=np.linalg.inv(info)
        score=xy-np.einsum('gij,j->gi',xx,beta)-sx*(a*(sy-sx@beta))[:,None]
        influence=score@inv
        influence*=np.sqrt(g/(g-1)*(n-1)/(n-p))
        return beta,influence,sse/(n-p)
    opt=minimize_scalar(lambda log_ratio:evaluate(np.exp(log_ratio)),
                        bounds=(-16,16),method='bounded',options={'xatol':1e-8})
    if not opt.success or not np.isfinite(opt.fun):
        return {'status':'not_converged','n':n,'n_species':g}
    ratio=float(np.exp(opt.x)) if opt.fun<evaluate(0.) else 0.
    beta,influence,sigma2=evaluate(ratio,True)
    return {'status':'estimated','n':n,'n_species':g,'labels':labels,
            'beta':beta,'influence':influence,'variance_ratio':ratio,
            'residual_variance':float(sigma2),
            'variance_boundary':bool(ratio==0 or abs(opt.x)>15.99),
            'objective':float(evaluate(ratio))}


def fit_module(frame, outcomes, scales, predictors=None):
    predictors=ENV if predictors is None else predictors
    terms=['intercept']+[v+s for v in predictors for s in ('__between','__within')]
    # Even reduced sensitivity models retain exactly the four-predictor complete rows.
    finite_environment=np.isfinite(frame[ENV]).all(axis=1)
    labels=np.sort(frame.accepted_key.unique())
    label_index={label:i for i,label in enumerate(labels)}
    fits=[];influences=[];components={}
    for outcome in outcomes:
        eligible=finite_environment & frame[outcome+'__eligible'].eq(True) & np.isfinite(frame[outcome])
        d=frame.loc[eligible]
        x=np.column_stack([np.ones(len(d))]+[
            ((d[v+s]-(scales[v]['mean'] if s=='__between' else 0))/scales[v]['sd']).to_numpy()
            for v in predictors for s in ('__between','__within')])
        y=d[outcome].to_numpy(float)
        center=float(y.mean()) if len(y) else 0.
        scale=float(y.std()) if len(y) else 0.
        fit=fit_response(x,(y-center)/scale if scale>0 else y,d.accepted_key.to_numpy())
        components[outcome]={k:v for k,v in fit.items() if k not in ('labels','beta','influence')}
        if fit['status']!='estimated':
            return {'status':'not_identifiable','components':components}
        influence=np.zeros((len(labels),len(terms)))
        influence[[label_index[a] for a in fit['labels']]]=fit['influence']
        influences.append(influence);fits.append(fit)
        se=np.sqrt(np.sum(fit['influence']**2,axis=0))
        critical=t.ppf(.975,fit['n_species']-1)
        components[outcome]['effects']={
            term:{'standardized_beta':float(fit['beta'][j]),
                  'raw_units_per_environment_sd':float(fit['beta'][j]*scale),
                  'raw_ci95':[float((fit['beta'][j]-critical*se[j])*scale),
                              float((fit['beta'][j]+critical*se[j])*scale)]}
            for j,term in enumerate(terms) if j}
        components[outcome]['response_sd']=scale
        components[outcome]['response_mean']=center
        components[outcome]['reference_prediction']=float(center+fit['beta'][0]*scale)
    scores=np.concatenate(influences,axis=1)
    covariance=scores.T@scores
    coefficients=np.concatenate([a['beta'] for a in fits])
    tests={}
    denominator=min(a['n_species'] for a in fits)-1
    def wald(ix):
        b=coefficients[ix];c=covariance[np.ix_(ix,ix)]
        sd=np.sqrt(np.maximum(np.diag(c),0))
        corr=c/np.outer(sd,sd) if np.all(sd>0) else np.zeros_like(c)
        if np.linalg.matrix_rank(corr)<len(ix):
            return {'status':'not_identifiable'}
        statistic=float((b/sd)@np.linalg.solve(corr,b/sd)/len(ix))
        return {'status':'estimated','F':statistic,'df1':len(ix),
                'df2':denominator,'p':float(f.sf(statistic,len(ix),denominator))}
    for j,term in enumerate(terms[1:],1):
        tests[term]=wald([i*len(terms)+j for i in range(len(outcomes))])
    joint_temperature_vpd={}
    if ENV[0] in predictors and ENV[2] in predictors:
        for suffix in ('__between','__within'):
            ix=[i*len(terms)+terms.index(v+suffix)
                for i in range(len(outcomes)) for v in (ENV[0],ENV[2])]
            joint_temperature_vpd[suffix]=wald(ix)
    hue={}
    if 'corolla_hue_sin' in components and 'corolla_hue_cos' in components:
        sine,cosine=components['corolla_hue_sin'],components['corolla_hue_cos']
        s0,c0=sine['reference_prediction'],cosine['reference_prediction']
        if np.hypot(s0,c0)>1e-6:
            for term in terms[1:]:
                s1=s0+sine['effects'][term]['raw_units_per_environment_sd']
                c1=c0+cosine['effects'][term]['raw_units_per_environment_sd']
                delta=np.rad2deg(np.arctan2(s1,c1)-np.arctan2(s0,c0))
                hue[term]={'reference_degrees':float(np.rad2deg(np.arctan2(s0,c0))%360),
                           'change_degrees':float((delta+180)%360-180),
                           'scope':'descriptive fixed-effect +1 native predictor SD contrast; not an independent test'}
    return {'status':'estimated','outcomes':outcomes,'components':components,'tests':tests,
            'predictors':predictors,'joint_temperature_vpd':joint_temperature_vpd,
            'conditional_hue_contrasts':hue,
            'inference':'approximate stacked marginal REML species-cluster Wald F',
            'random_effects':'response-specific species intercept; no phylogenetic correction',
            'pointwise_component_intervals_not_joint_significance_tests':True}


def execute(table, receipt, output):
    report=json.loads(Path(receipt).read_text())
    payload=Path(table).read_bytes()
    if hashlib.sha256(payload).hexdigest()!=report['prepared_table_sha256']:
        raise ValueError('Prepared input SHA mismatch')
    frame=pd.read_csv(table,dtype={'obs_id':str},low_memory=False)
    if len(frame)!=report['native_observations'] or not frame.obs_id.is_unique:
        raise ValueError('Observation population changed')
    scales={v:{'mean':float(frame[v].mean()),'sd':float(frame[v].std(ddof=0))} for v in ENV}
    if any(a['sd']<=0 for a in scales.values()):raise ValueError('Constant climate predictor')
    ranking=frame.accepted_key.value_counts().index
    result={'input':report,'scales':scales,'models':{},
            'actual_design_diagnostics':design_diagnostics(frame),
            'model_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'primary_family':'24 joint tests: 3 modules x 4 predictors x within/between; Holm',
            'scope':'retrospective observational associations; not adaptation or measured defence',
            'limitations':['working random-intercept covariance and conditional cluster inference',
                          'no broad spatial or phylogenetic adjustment in this primary stage',
                          'legacy image colour is not calibrated reflectance or anthocyanin',
                          'bract projection and spread fraction are not anatomical opening angles',
                          'missing measurements are not imputed or inferred to be absent traits']}
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for drop in (0,2,5):
        d=frame.loc[~frame.accepted_key.isin(ranking[:drop])]
        for module,outcomes in MODULES.items():
            print(f'Fitting {module}, exclude {drop} most represented native taxa',flush=True)
            fit=fit_module(d,outcomes,scales)
            fit['excluded_taxa_count']=drop
            fit['excluded_taxa']=[str(frame.loc[frame.accepted_key.eq(taxon),'taxon_name'].mode().iloc[0])
                                  for taxon in ranking[:drop]]
            result['models'][f'{module}/drop{drop}']=fit
    result['attribution_sensitivity']={}
    for omitted in (ENV[0],ENV[2]):
        for module,outcomes in MODULES.items():
            print(f'Attribution sensitivity {module}, omit {omitted}; fixed rows',flush=True)
            fit=fit_module(frame,outcomes,scales,[v for v in ENV if v!=omitted])
            if fit['status']=='estimated':
                for outcome,c in fit['components'].items():
                    if c['n']!=result['models'][module+'/drop0']['components'][outcome]['n']:
                        raise ValueError('Attribution comparison observation set changed')
            result['attribution_sensitivity'][module+'/omit_'+omitted]=fit
    block_family=[test for name,model in result['models'].items() if name.endswith('/drop0')
                  for test in model.get('joint_temperature_vpd',{}).values()]
    # Missing block tests do not shrink the prespecified six-test family.
    block_p=[test.get('p',1.) for test in block_family]
    if len(block_family)==6:
        for test,adjusted in zip(block_family,multipletests(block_p,method='holm')[1]):
            if test['status']=='estimated':test['holm_p_six_block_tests']=float(adjusted)
    result['attribution_scope']='retrospective six joint temperature-VPD block tests; '
    result['attribution_scope']+='omission models are same-row diagnostics, not preferred causal models'
    family=[test for name,model in result['models'].items() if name.endswith('/drop0')
            for test in model.get('tests',{}).values()]
    if len(family)!=24 or any(test['status']!='estimated' for test in family):
        result['primary_status']='incomplete_inference'
    else:
        for test,q in zip(family,multipletests([test['p'] for test in family],method='holm')[1]):
            test['holm_p']=float(q)
        result['primary_status']='estimated'
    (output/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    rows=[]
    for name,model in result['models'].items():
        for term,test in model.get('tests',{}).items():
            rows.append({'model':name,'term':term,**test})
    pd.DataFrame(rows).to_csv(output/'joint_tests.csv',index=False)
    rows=[]
    for name,model in result['models'].items():
        for outcome,component in model.get('components',{}).items():
            for term,effect in component.get('effects',{}).items():
                rows.append({'model':name,'outcome':outcome,'term':term,
                             'n':component['n'],'n_species':component['n_species'],**effect})
    pd.DataFrame(rows).to_csv(output/'component_effects.csv',index=False)
    write_report(result,output)
    return result


def write_report(result, output):
    labels={'colour':'色全体','orientation':'向き','bract_projection':'総苞片の輪郭突出'}
    names=dict(zip(ENV,['温度','日射','VPD','降水']))
    lines=['# v4 観測単位解析：GEBを想定した主解析結果','',
           '温度・日射・VPD・降水を同時に含め、各観測を保持した種ランダム切片モデル。'
           '環境を種平均と種内偏差に分けた。色は明度・彩度・円周色相の共同検定、'
           '総苞片は突出量と広がり割合の共同検定である。','',
           '**この結果は環境との条件付き関連であり、適応・可塑性・防御機構の直接証明ではない。**','',
           '## 主結果','',
           '以下は24共同検定のHolm補正P値。成分別の有意性を拾って色全体の結論にしていない。','',
           '|形質群|環境|種間 P（Holm）|種内 P（Holm）|',
           '|---|---|---:|---:|']
    for module,label in labels.items():
        model=result['models'][module+'/drop0']
        for v in ENV:
            vals=[model.get('tests',{}).get(v+s,{}).get('holm_p') for s in ('__between','__within')]
            text=[f'{x:.4g}' if x is not None else '評価不能' for x in vals]
            lines.append(f'|{label}|{names[v]}|{text[0]}|{text[1]}|')
    lines+=['','## 観測数と方向','',
            '効果量は、全自生域母集団での環境1標準偏差当たり。区間は成分ごとの'
            '点別95%区間で、共同検定の代わりではない。','',
            '|座標|観測数|分類群数|','|---|---:|---:|']
    for module in labels:
        for outcome,c in result['models'][module+'/drop0']['components'].items():
            lines.append(f"|{outcome}|{c['n']:,}|{c['n_species']}|")
    lines+=['','### 色の全体像','',
            'ΔL*、ΔC*と円周色相の変化を併記する。色相の度数差は、種平均環境が'
            '母集団平均、種内偏差が0、ランダム効果が0の基準点から、1変数のみを'
            '+1標準偏差動かした記述用予測。光源校正済み反射色や色素量ではない。'
            '強く相関する環境を独立に動かす条件付き比較であり、自然に起きる気候変化の予測ではない。','',
            '|環境|尺度|ΔL*［95%区間］|ΔC*［95%区間］|色相差（度、記述用）|',
            '|---|---|---|---|---:|']
    colour=result['models']['colour/drop0']
    for v in ENV:
        for suffix,label in [('__between','種間'),('__within','種内')]:
            term=v+suffix;cells=[]
            for outcome in ['corolla_lab_lightness','corolla_lab_chroma']:
                e=colour['components'][outcome]['effects'][term];lo,hi=e['raw_ci95']
                cells.append(f"{e['raw_units_per_environment_sd']:+.2f} [{lo:+.2f}, {hi:+.2f}]")
            hue=colour.get('conditional_hue_contrasts',{}).get(term,{}).get('change_degrees')
            h=f'{hue:+.2f}' if hue is not None else '評価不能'
            lines.append(f'|{names[v]}|{label}|{cells[0]}|{cells[1]}|{h}|')
    lines+=['','### 向き','',
            '|環境|尺度|画像基準角の変化（度）［95%区間］|','|---|---|---:|']
    orientation=result['models']['orientation/drop0']['components']['orientation_image_vertical_angle']
    for v in ENV:
        for suffix,label in [('__between','種間'),('__within','種内')]:
            e=orientation['effects'][v+suffix];lo,hi=e['raw_ci95']
            lines.append(f"|{names[v]}|{label}|{e['raw_units_per_environment_sd']:+.2f} [{lo:+.2f}, {hi:+.2f}]|")
    lines+=['','## 特定種への依存','',
            '除外対象は形質の結果を見ず、全自生域観測数の上位2種・5種で固定した。'
            '下表の除外後P値は未補正の感度分析値。主結果のHolm補正Pと混同しない。','',
            '|形質群|環境・尺度|全種（未補正P）|上位2種除外|上位5種除外|',
            '|---|---|---:|---:|---:|']
    for module,label in labels.items():
        for v in ENV:
            for suffix,level in [('__between','種間'),('__within','種内')]:
                cells=[]
                for drop in (0,2,5):
                    p=result['models'][f'{module}/drop{drop}'].get('tests',{}).get(v+suffix,{}).get('p')
                    cells.append(f'{p:.4g}' if p is not None else '評価不能')
                lines.append(f'|{label}|{names[v]}・{level}|'+ '|'.join(cells)+'|')
    lines+=['','上位5種：'+ '、'.join(result['models']['colour/drop5']['excluded_taxa'])+'。','',
            '## 結論の限界とGEBへの位置づけ','',
            '- 総苞片の指標は輪郭突出・広がりであり、解剖学的な開出角や防御効果ではない。'
            '明確な関連がない場合も、総苞片の防御仮説そのものの否定にはならない。',
            '- ランダム切片分散が0の境界に達した成分はresults.jsonに記録した。'
            'その成分について種間異質性が明確に推定できたとは扱わない。',
            '- 種内と種間は異なる尺度の条件付き傾向。種間の差を同一個体・同一種の'
            '温暖化応答や可塑性に読み替えない。すべての種が同じ傾きを持つとも断定しない。',
            '- 主解析は広域空間構造・系統を同時調整していない。土壌を含む追加検証も別工程。'
            '本結果のみで投稿準備完了とはしない。',
            '- GEBの主題は、連続的な写真形質によって種カテゴリー内の変異まで比較でき、'
            '異なる生物学的尺度の環境関連を区別できること。解析手法の利用自体を新規性とはしない。',
            '',
            '[GEB Aims and Scope](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/aims.htm)',
            '',
            '再現情報：入力・測定可否・環境相関・コードSHAはresults.json、'
            '全効果量はcomponent_effects.csv、全共同検定はjoint_tests.csv。']
    lines+=['','## 実際のモデル行列における環境相関','',
            '下表は明度の解析対象行。種間は観測数による重み付きの診断である。'
            '相関・VIFは説明変数の重なりを表し、未測定交絡の除去や原因の識別を保証しない。'
            '全座標の診断はresults.jsonに保存した。','',
            '|尺度|温度–VPD相関|温度VIF|日射VIF|VPD VIF|降水VIF|','|---|---:|---:|---:|---:|---:|']
    diag=result['actual_design_diagnostics']['endpoints']['corolla_lab_lightness']
    for key,label in [('between_observation_weighted','種間'),('within','種内')]:
        a=diag[key]
        if a['status']=='estimated':
            lines.append(f"|{label}|{a['correlation'][0][2]:.3f}|"+
                         '|'.join(f'{x:.2f}' for x in a['vif'])+'|')
    full=diag['full_eight_terms']
    if full['status']=='estimated':
        lines+=['',f"8項全体の標準化条件指数：{full['standardized_condition_index']:.2f}。"
                '8項全体でのVIFはresults.jsonを参照。']
    lines+=['','## 温度とVPDの寄与を分けられるか','',
            '4環境モデルの温度・VPDを共同で検定し、日射・降水を調整した'
            'この2変数ブロックの関連を評価した。これは独立した温度効果やVPD効果の証明ではない。'
            '以下は追加6検定内でのHolm補正。主解析24検定の補正とは別である。','',
            '|形質群|種間 温度＋VPD|種内 温度＋VPD|','|---|---:|---:|']
    for module,label in labels.items():
        block=result['models'][module+'/drop0'].get('joint_temperature_vpd',{})
        cells=[]
        for suffix in ('__between','__within'):
            p=block.get(suffix,{}).get('holm_p_six_block_tests')
            cells.append(f'{p:.4g}' if p is not None else '評価不能')
        lines.append(f'|{label}|'+ '|'.join(cells)+'|')
    lines+=['','温度除外・VPD除外モデルは、各座標の解析対象行を固定して再推定した。'
            '省略すると他方の効果に共通の関連が混ざり得るため、有意性を理由に主モデルを'
            '置き換えない。係数・区間・共同検定はresults.jsonのattribution_sensitivityに保存した。']
    (Path(output)/'conclusions_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('table','receipt','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();execute(a.table,a.receipt,a.output)

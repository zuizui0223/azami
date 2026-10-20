"""One fixed environmental response analysis, without specification searching."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from analysis.v4.environment_candidate_audit import matrix_summary
from analysis.v4.observation_models import MODULES, fit_module, reuse_inventory

LABELS={'colour':'色全体','orientation':'向き','bract_projection':'総苞片の輪郭突出'}
ENV_LABELS={'chelsa_bio01':'温度','chelsa_rsds_mean':'日射','chelsa_vpd_mean':'VPD',
            'chelsa_bio12':'降水','chelsa_bio04':'温度季節性','chelsa_bio15':'降水季節性',
            'chelsa_sfcwind_mean':'風速','chelsa_gsp':'生育期降水量','chelsa_npp':'NPP'}


def execute(table, receipt_path, output, selection_path):
    receipt=json.loads(receipt_path.read_text())
    digest=hashlib.sha256(table.read_bytes()).hexdigest()
    if digest!=receipt['prepared_table_sha256']:
        raise ValueError('Input checksum mismatch')
    selection=json.loads(selection_path.read_text())
    if selection['input']['input_sha256']!=digest or selection['status']!='selected_before_trait_models':
        raise ValueError('Selection does not match verified observations')
    predictors=selection['selected']
    if not isinstance(predictors,list) or not predictors or len(set(predictors))!=len(predictors):
        raise ValueError('Invalid selected predictor list')
    if not set(predictors)<=set(selection['input']['candidates']):
        raise ValueError('Selected variable outside candidate universe')
    frame=pd.read_csv(table,low_memory=False)
    if len(frame)!=receipt['native_observations'] or not frame.obs_id.is_unique:
        raise ValueError('Native population changed')
    scales={}
    for v in predictors:
        frame[v]=frame[v].where(np.isfinite(frame[v]))
        scales[v]={'mean':float(frame[v].mean()),'sd':float(frame[v].std(ddof=0))}
        if not np.isfinite(scales[v]['sd']) or scales[v]['sd']<=0:
            raise ValueError('No finite variation: '+v)
        frame[v+'__between']=frame.groupby('accepted_key')[v].transform('mean')
        frame[v+'__within']=frame[v]-frame[v+'__between']
    diagnostics={}
    complete=np.isfinite(frame[predictors]).all(axis=1)
    for outcomes in MODULES.values():
        for endpoint in outcomes:
            rows=frame.loc[complete & frame[endpoint+'__eligible'].eq(True) & np.isfinite(frame[endpoint])]
            diagnostics[endpoint]={'n':len(rows),'designs':{
                level:matrix_summary(rows,columns) for level,columns in [
                    ('between',[v+'__between' for v in predictors]),
                    ('within',[v+'__within' for v in predictors]),
                    ('full',[v+s for v in predictors for s in ('__between','__within')])]}}
    models={}
    for module,outcomes in MODULES.items():
        print('Fixed environmental response: '+module,flush=True)
        models[module]=fit_module(frame,outcomes,scales,predictors,required_environment=predictors)
        for endpoint,c in models[module]['components'].items():
            if c['n']!=diagnostics[endpoint]['n']:
                raise ValueError('Endpoint-specific sample count changed')
    tests=[models[m].get('tests',{}).get(v+s,{})
           for m in MODULES for v in predictors for s in ('__between','__within')]
    for test,q in zip(tests,multipletests([t.get('p',1.) for t in tests],method='holm')[1]):
        if test.get('status')=='estimated':
            test['holm_p']=float(q)
    crosswalk=Path(__file__).with_name('trait27_functional_crosswalk_v1.json')
    inventory=reuse_inventory(frame,receipt,json.loads(crosswalk.read_text()))
    result={'question':'Which environmental gradients are associated with head traits, '
                        'conditional on the other included environmental dimensions?',
            'input_sha256':digest,'input_receipt':receipt,
            'git_sha':os.environ.get('GITHUB_SHA'),
            'run_id':os.environ.get('GITHUB_RUN_ID'),
            'model_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'fit_code_sha256':hashlib.sha256(Path(__file__).with_name('observation_models.py').read_bytes()).hexdigest(),
            'selection':selection,'selection_sha256':hashlib.sha256(selection_path.read_bytes()).hexdigest(),
            'family_size':len(tests),'predictors':predictors,'scales':scales,'models':models,
            'environment_parameterization':'within_between',
            'environment_diagnostics':diagnostics,'reuse_inventory':inventory,
            'inference_status':'estimated' if all(t.get('status')=='estimated' for t in tests)
                               else 'incomplete_inference',
            'scope':'Retrospective conditional associations after environment-only vifstep selection. No '
                    'causal identification, spatial correction or phylogenetic correction.'}
    output.mkdir(parents=True,exist_ok=True)
    (output/'environment_response.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    rows=[]
    for module,model in models.items():
        for endpoint,c in model.get('components',{}).items():
            for term,e in c.get('effects',{}).items():
                rows.append({'module':module,'endpoint':endpoint,'term':term,
                             'n':c['n'],'species_n':c['n_species'],**e})
    pd.DataFrame(rows).to_csv(output/'environment_response_effects.csv',index=False)
    pd.DataFrame(inventory).to_csv(output/'measurement_reuse_inventory.csv',index=False)
    report(result,output)
    return result


def report(result,output):
    if result.get('environment_parameterization')=='common_slope':
        return report_common(result,output)
    lines=['# 頭花形質の環境勾配応答','',
           '問い：他の環境条件を同時に考慮したとき、どの環境勾配に沿って頭花形質が変わるか。','',
           '6群の全9候補から環境データのみでvifstep選定し、採択セットを固定して同時投入：'
           +'、'.join(ENV_LABELS[v] for v in result['predictors'])+'。'
           '観測値を平均せず、環境を種平均と種内偏差に分け、成分ごとに種ランダム切片を推定した。'
           '色全体は明度・彩度・円周色相の共同検定。各成分の利用可能な観測を使う。','',
           '過去の結果を閲覧済みの改訂観測解析であり、事前登録した独立検証ではない。'
           f"{result['family_size']}共同検定でHolm補正。成分の95%区間は点別である。",'',
           '主解析は色全体と向き。総苞片は測定数・被覆範囲と機能的指標の制約からSIに置く。'
           '移動を理由に48検定の補正をやり直していない。','',
           '|形質|環境|種間P（Holm）|種内P（Holm）|','|---|---|---:|---:|']
    for module in ('colour','orientation'):
        label=LABELS[module]
        for v in result['predictors']:
            cells=[]
            for s in ('__between','__within'):
                p=result['models'][module].get('tests',{}).get(v+s,{}).get('holm_p')
                cells.append(f'{p:.4g}' if p is not None else '評価不能')
            lines.append(f'|{label}|{ENV_LABELS[v]}|'+ '|'.join(cells)+'|')
    lines+=['','## 各成分の変化','',
            '全自生域母集団の環境1標準偏差当たりの変化。単独の成分検定で'
            '色全体の関連を判定しない。色相sin/cosは共同で読む。','',
            '|成分|観測数|種数|環境・尺度|変化［95%区間］|','|---|---:|---:|---|---|']
    for module in ('colour','orientation'):
        model=result['models'][module]
        for endpoint,c in model.get('components',{}).items():
            for term,e in c.get('effects',{}).items():
                lo,hi=e['raw_ci95']
                lines.append(f"|{endpoint}|{c['n']}|{c['n_species']}|{term}|"
                             f"{e['raw_units_per_environment_sd']:+.3g} [{lo:+.3g}, {hi:+.3g}]|")
    lines+=['','## 環境の重なりと解釈の限界','',
            '相関・VIF・欠測数はJSONに全成分分を保存。高いVIFを理由に結果を見て'
            '追加で変数を削除しない。選定は母集団の環境値で行ったため、種間・種内や'
            '測定対象群では相関が残ることがある。相関が強い場合、個別係数は不安定であり、'
            '独立した環境作用を識別できたとは限らない。','',
            '種内の関連は可塑性の証明ではなく、種間の関連も適応の証明ではない。'
            '種内と種間の係数差を直接検定していないため、片方だけの有意性を'
            '尺度間の差の証拠にはしない。'
            '空間・系統・未測定環境を調整した結論ではない。'
            '写真の色は色素量ではなく、向きは画像基準角、総苞片は輪郭突出の代理指標である。','',
            f"再現情報：Actions run {result['run_id']}、commit {result['git_sha']}。"]
    (output/'environment_response_ja.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    si=['# SI：総苞片の輪郭突出・広がり','',
        '既存測定1,303観測の補足解析。主解析より観測規模・対象被覆が限定されるためSIに置く。'
        'この判断は有意性による選別ではなく、測定範囲と指標の制約に基づく。'
        '既存の48検定Holm補正を維持する。','',
        '輪郭突出・広がりは個々の総苞片の開出角や棘による防御効果そのものではない。'
        '気候との関連が明確でなくても、食害者・訪花者に対する防御機能の否定にはならない。','',
        '|環境|種間P（Holm）|種内P（Holm）|','|---|---:|---:|']
    model=result['models']['bract_projection']
    for v in result['predictors']:
        values=[model.get('tests',{}).get(v+s,{}).get('holm_p') for s in ('__between','__within')]
        si.append('|'+ENV_LABELS[v]+'|'+'|'.join(f'{p:.4g}' if p is not None else '評価不能' for p in values)+'|')
    si+=['','|成分|観測数|種数|環境・尺度|変化［点別95%区間］|','|---|---:|---:|---|---|']
    for endpoint,c in model.get('components',{}).items():
        for term,e in c.get('effects',{}).items():
            lo,hi=e['raw_ci95']
            si.append(f"|{endpoint}|{c['n']}|{c['n_species']}|{term}|"
                      f"{e['raw_units_per_environment_sd']:+.3g} [{lo:+.3g}, {hi:+.3g}]|")
    si+=['','全係数・診断・来歴はenvironment_response.jsonに保存。再測定はこの結果の成立条件にしない。',
         f"推定元：Actions run {result['run_id']}、commit {result['git_sha']}。"]
    (output/'environment_response_si_ja.md').write_text('\n'.join(si)+'\n',encoding='utf-8')


def report_common(result,output):
    for filename,modules,title in [
        ('environment_response_ja.md',('colour','orientation'),'主解析：色と向き'),
        ('environment_response_si_ja.md',('bract_projection',),'SI：総苞片の輪郭突出')]:
        lines=['# '+title,'',
               'モデル：観測ごとの形質 ～ 採択環境の同時投入 ＋ (1 | 種)。',
               '形質も環境も種平均には集約せず、環境の種内・種間分解も行わない。'
               '係数は種内・種間双方の情報に基づく共通傾きであり、どちらか一方の効果ではない。','',
               '採択環境：'+'、'.join(ENV_LABELS[v] for v in result['predictors'])+'。',
               f"主解析とSIを合わせた{result['family_size']}共同検定でHolm補正。"
               '種ランダム切片と種クラスタに基づく区間・共同検定を使用。','',
               '|形質|環境|共同検定P（Holm）|','|---|---|---:|']
        for module in modules:
            model=result['models'][module]
            for v in result['predictors']:
                p=model.get('tests',{}).get(v,{}).get('holm_p')
                lines.append(f"|{LABELS[module]}|{ENV_LABELS[v]}|"
                             +(f'{p:.4g}' if p is not None else '評価不能')+'|')
        lines+=['','## 効果量','',
                '環境が自生域母集団の1標準偏差増えたときの変化。区間は点別95%。'
                '成分の有意性を色全体の共同検定の代わりにしない。','',
                '|成分|観測数|種数|環境|変化［95%区間］|','|---|---:|---:|---|---|']
        for module in modules:
            for endpoint,c in result['models'][module].get('components',{}).items():
                for term,e in c.get('effects',{}).items():
                    lo,hi=e['raw_ci95']
                    lines.append(f"|{endpoint}|{c['n']}|{c['n_species']}|{ENV_LABELS.get(term,term)}|"
                                 f"{e['raw_units_per_environment_sd']:+.3g} [{lo:+.3g}, {hi:+.3g}]|")
        lines+=['','種の違いをランダム切片で扱うことは、系統・空間や未測定交絡の調整を保証しない。'
                'ランダム切片と説明変数の独立性などのモデル仮定の下での関連であり、'
                '因果効果、適応、可塑性や個々の種に共通する応答の証明ではない。',
                '総苞片は測定規模と代理指標の制約からSIとする。アントシアニン濃度・防御効果は未測定。',
                f"推定元：Actions run {result['run_id']}、commit {result['git_sha']}。"]
        (output/filename).write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('table','receipt','output','selection'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    execute(a.table,a.receipt,a.output,a.selection)


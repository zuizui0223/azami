"""Audit already-published anonymous QC summaries; never access original data.

This does not fit an RF, create a significance test, pool seeds, or authenticate
biological trait values. It distinguishes two existing uncertainty diagnostics.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import zipfile

TRAITS=('angle','lightness','chroma','hue','aspect','circularity','solidity','width_profile')
SEEDS=(20261101,20261103,20261105)
SOURCE_NAME='rf_photo_observer_disjoint_original_S1_1deg_5deg_aggregate.json'
SOURCE_SHA256='ac58ce0e58b7cfaa40fcab994e0852b5c06899673f06c96b4434df2c5d5696f5'
LABELS={'angle':'画像基準の頭花角度','lightness':'明度 L*','chroma':'彩度 C*',
        'hue':'色相（sin/cos共同）','aspect':'輪郭縦横比','circularity':'円形度',
        'solidity':'Solidity','width_profile':'幅プロファイルCV'}


def _finite(value):
    value=float(value)
    if not math.isfinite(value):raise ValueError('NONFINITE_AGGREGATE_VALUE')
    return value


def _range(values):
    values=list(values)
    return [min(values),max(values)]


def audit(data: dict) -> dict:
    if (data.get('schema')!='azami_original_S1_global_observer_disjoint_split_half_v1'
        or data.get('source_original_S1_observations')!=313350
        or data.get('distinct_authenticated_source_observers')!=92839):
        raise ValueError('UNEXPECTED_ORIGINAL_QC_SOURCE')
    rows=data.get('results',[])
    expected={(t,g) for t in TRAITS for g in (1,5)}
    keys=[(r.get('construct'),r.get('degrees')) for r in rows]
    if len(keys)!=16 or set(keys)!=expected:raise ValueError('INCOMPLETE_OR_DUPLICATE_ALL8_GRAINS')
    output=[]
    for row in rows:
        splits=row.get('splits',[])
        if len(splits)!=3 or {s.get('seed') for s in splits}!=set(SEEDS):
            raise ValueError('INCOMPLETE_OR_DUPLICATE_SEEDS')
        values=[]
        for s in splits:
            if s.get('status')!='estimated_exploratory_reliability':
                raise ValueError('SOURCE_CONDITION_NOT_ESTIMATED')
            n=int(s['n_taxa']);r=_finite(s['taxon_equal_spatial_alignment']['q50'])
            lo,hi=map(_finite,s['median_alignment_taxon_bootstrap_95'])
            null=_finite(s['within_taxon_cell_labels_permuted_null_median_alignment']['q95'])
            fraction=_finite(s['fraction_taxa_alignment_at_least_0_5'])
            positive=_finite(s['fraction_taxa_alignment_positive'])
            count=round(n*fraction)
            if (n<1 or not -1<=r<=1 or not -1<=lo<=hi<=1 or not -1<=null<=1
                or not 0<=fraction<=1 or not 0<=positive<=1
                or not math.isclose(n*fraction,count,abs_tol=1e-8)):
                raise ValueError('INVALID_TAXON_DISTRIBUTION_OR_COUNT')
            values.append(dict(seed=s['seed'],n=n,cells=int(s['n_qualified_cells']),
                photos=int(s['n_contributing_photos']),median=r,lo=lo,hi=hi,
                null_q95=null,large_count=count,large_fraction=fraction,positive=positive))
        output.append({'construct':row['construct'],'degrees':int(row['degrees']),
            'n_splits':3,'n_taxa_range':_range(v['n'] for v in values),
            'n_cells_range':_range(v['cells'] for v in values),
            'n_contributing_photo_records_range':_range(v['photos'] for v in values),
            'median_alignment_range':_range(v['median'] for v in values),
            'median_above_zero_count':sum(v['median']>0 for v in values),
            'median_above_matched_permutation_q95_count':sum(v['median']>v['null_q95'] for v in values),
            'taxon_bootstrap_lower_above_zero_count':sum(v['lo']>0 for v in values),
            'taxa_alignment_ge_0_5_count_range':_range(v['large_count'] for v in values),
            'taxa_alignment_ge_0_5_fraction_range':_range(v['large_fraction'] for v in values),
            'positive_taxon_fraction_range':_range(v['positive'] for v in values),
            'fixed_split_summaries':values})
    output.sort(key=lambda r:(TRAITS.index(r['construct']),r['degrees']))
    return {'schema':'azami_all8_existing_observer_QC_evidence_audit_v1',
            'new_RF_fit_executed':False,'seeds_are_independent_experiments':False,
            'original_result_values_changed':False,
            'scope':'Recomputation of summary counts from previously executed anonymous QC only.',
            'limitations':[
                'Permutation-reference exceedance is not an adjusted P value.',
                'Cell-label permutations do not preserve geographic autocorrelation.',
                'Taxon bootstrap includes taxon-composition uncertainty but not phylogenetic independence.',
                'Different seeds/grains can contain different taxa; high-correlation identities cannot be tracked from aggregates.',
                'No biological measurement-error, adaptation, or environmental-prediction claim follows from these QC summaries.'
            ],'conditions':output}


def markdown(report: dict) -> str:
    lines=['# 全8形質：既存の撮影者非重複QCを再集計（2026-11-06）','',
        '**新しいRF適合ではない。** 2026-11-01の匿名集計を原本SHA256で検証し、既存の点推定・帰無分布・bootstrapを別々に整理した。',
        '', '|形質|粒度|分類群|一致度中央値（3分割範囲）|対応する置換帰無95%点より上|bootstrap下限>0|一致度≥0.5の分類群数|',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in report['conditions']:
        lo,hi=r['median_alignment_range'];n=r['n_taxa_range'];k=r['taxa_alignment_ge_0_5_count_range']
        lines.append(f"|{LABELS[r['construct']]}|{r['degrees']}°|{n[0]}–{n[1]}|{lo:+.3f}–{hi:+.3f}|{r['median_above_matched_permutation_q95_count']}/3|{r['taxon_bootstrap_lower_above_zero_count']}/3|{k[0]}–{k[1]}|")
    lines.extend(['', '## 解釈', '',
        '頭花角度は1°・5°の各3分割すべてで、点推定が対応するセル置換帰無の95%点を上回る。一方、分類群bootstrapの95%区間はすべてゼロを含む。固定された観測集合のセル対応に関する探索的な手がかりと、分類群を越えた安定性は別の問いであり、前者を新しい有意な一般則にしない。',
        '', '花色の3構成は両診断で各粒度3/3。輪郭4構成は全分割で支持されず、円形度5°の一部分割などの例外を一般化しない。向きや形を以後の解析から除外する理由にも使わない。',
        '', '一致度≥0.5の分類群数は記述的な集計であって有意な分類群の数ではない。分類群IDを含まないため、別分割でも同じ分類群が高いかは未判定。高い分類群だけを選んで再解析する設計にはしない。',
        '', '3分割は重なるため独立した3実験ではない。1°と5°の分類群集合も異なる。セルラベル置換は地理的自己相関を保存しない。ここでの点検は未校正画像の空間パターンのQCであり、環境による予測や色素・姿勢の生物学的妥当性を検証してはいない。',
        '', '## 出典', '',
        '既存実行: https://github.com/zuizui0223/azami/actions/runs/41261048854',
        '既存artifact: 13065605960',
        f'元匿名JSON SHA256: `{SOURCE_SHA256}`',
        '', '入力に個別観測、場所、撮影者・分類群ID、秘密情報はない。凍結原稿・既存RFの数値は変更していない。'])
    return '\n'.join(lines)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-zip',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    with zipfile.ZipFile(args.artifact_zip) as archive:
        if archive.namelist().count(SOURCE_NAME)!=1:raise ValueError('MISSING_OR_DUPLICATE_SOURCE_MEMBER')
        raw=archive.read(SOURCE_NAME)
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256:raise ValueError('SOURCE_AGGREGATE_SHA_MISMATCH')
    result=audit(json.loads(raw));result['source_aggregate_sha256']=SOURCE_SHA256
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'all8_observer_qc_evidence_audit_20261106.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (args.output_dir/'20261106_ALL8_OBSERVER_QC_EVIDENCE_AUDIT_JA.md').write_text(markdown(result),encoding='utf-8')
    print('Verified aggregate-only evidence audit: 16 conditions, 48 split summaries; no RF fit.')


if __name__=='__main__':main()

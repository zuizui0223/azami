"""Summarize source-authenticated four-V4-model uncorrected, BH and Holm tests.

No raw observation identifiers or pixel/image data. Hue sin/cos are ONE
biological endpoint via the joint two-degree-of-freedom cluster Wald test.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path
from analysis.v4.s1_s2_vif5_si10_alltraits import SPEC, MAIN, SI

EXPECTED={
    "S1": {"four_bio12":80,"four_bio18":80,"main_vif5":100,"si_vif10":140},
    "S2": {"four_bio12":128,"four_bio18":128,"main_vif5":160,"si_vif10":224}}


def domain(c):
    if c=="orientation_image_vertical_angle":
        return "orientation"
    if c in ("corolla_lab_lightness","corolla_lab_chroma","hue_joint_sin_cos"):
        return "colour"
    if c.startswith("capitulum_"):
        return "outline"
    if c.startswith("bract_"):
        return "bract_photo_proxy"
    raise ValueError("UNREGISTERED_OBSERVABLE_TRAIT "+c)


def analyze(d):
    if (d.get("status")!="SOURCE_VERIFIED_S1_S2_FOUR_MODELS_ALL_AVAILABLE_PHENOTYPES"
        or not d.get("V4_NOT_V3")
        or d["S1"]["original_native_n"]!=313350
        or d["S2"]["original_wild_source_n"]!=297622):
        raise ValueError("WRONG_SOURCE_OR_ANALYSIS")
    if (set(d["selection"]["selection"]["main_vif5"]["retained"])!=set(MAIN)
        or set(d["selection"]["selection"]["si_vif10"]["retained"])!=set(SI)):
        raise ValueError("VIF_SELECTION_CHANGED")
    summary, grouped={}, {"S1":defaultdict(dict),"S2":defaultdict(dict)}
    for study, models in EXPECTED.items():
        for model, expected in models.items():
            f=d["multiplicity"][study+"/"+model]
            if f["n_tests"]!=expected or len(f["rows"])!=expected:
                raise ValueError("EXPECTED_BIOLOGICAL_TEST_COUNT_MISSING "+study+model)
            types=defaultdict(lambda:{"tested":0,"raw":0,"BH":0,"Holm":0})
            for r in f["rows"]:
                if r["coordinate"] in ("corolla_hue_sin","corolla_hue_cos"):
                    raise ValueError("HUE_COMPONENTS_WRONGLY_COUNTED")
                key=(r["group"],r["coordinate"],r["environment"],r["scale"])
                if model in grouped[study][key]:
                    raise ValueError("DUPLICATE_ENDPOINT_MODEL")
                t=domain(r["coordinate"])
                types[t]["tested"]+=1
                types[t]["raw"]+=int(r["raw_lt_0_05"])
                types[t]["BH"]+=int(r["BH_lt_0_05"])
                types[t]["Holm"]+=int(r["Holm_lt_0_05"])
                grouped[study][key][model]=r
            source_vif=d["selection"]["fixed_model_set_diagnostics"][model]
            summary[study+"/"+model]={
                "tested":expected,"raw":f["n_unadjusted"],
                "BH":f["n_BH"],"Holm":f["n_Holm"],
                "BH_vif_eligible":f["n_BH_with_model_VIF_pass"],
                "raw_vif_eligible":f["n_raw_with_model_VIF_pass"],
                "max_source_vif":source_vif["max_vif"],
                "source_vif_pass_at_5":source_vif["max_vif"]<=5,
                "by_trait_group":dict(types)}
    shared={}
    for study,records in grouped.items():
        r=[]
        for (group,coordinate,exposure,scale),model_rows in records.items():
            if len(model_rows)<2:
                continue
            est=[float(v["beta_std"]) for v in model_rows.values()
                 if v["beta_std"] is not None]
            sign_consistent=(not est or all(x>0 for x in est)
                             or all(x<0 for x in est))
            r.append({
                "group":group,"coordinate":coordinate,
                "environment":exposure,"scale":scale,
                "trait_group":domain(coordinate),
                "n_models":len(model_rows),
                "sign_consistent":sign_consistent,
                "n_raw_sig":sum(v["raw_lt_0_05"] for v in model_rows.values()),
                "n_BH_sig":sum(v["BH_lt_0_05"] for v in model_rows.values()),
                "values":{model:{
                   "beta_std":v["beta_std"],
                   "ci95_wald":v["ci95_wald"],
                   "p_raw":v["p_original"],"q_BH":v["p_BH"],
                   "p_Holm":v["p_Holm"],
                   "n":v["n"],
                   "trait_cohort_vif_pass":v["trait_cohort_vif_pass"]}
                   for model,v in model_rows.items()}})
        shared[study]=r
    return {"status":"VERIFIED_FOUR_MODEL_CORRECTED_VS_UNCORRECTED",
        "model_order":list(SPEC),
        "model_predictors":{k:list(v) for k,v in SPEC.items()},
        "counts":summary,
        "shared_effects":shared,
        "S1_joint8_n":d["S1"]["joint_8_constructs_28_pairs"]["n"],
        "S1_joint8_taxa":d["S1"]["joint_8_constructs_28_pairs"]["taxa"],
        "S1_28_pairs":d["S1"]["joint_8_constructs_28_pairs"]["models"],
        "S1_bract_proxy_pairs":d["S1"]["bract_to_all_8_constructs"],
        "S2_bract":d["S2"]["optional_bract_source"]["status"],
        "limitations":{"pointwise_95CI_not_simultaneous":True,
            "same_observations_across_models_not_independent_replication":True,
            "S1_spatial_phylogeny_not_controlled":True,
            "S2_native_introduced_not_evolutionary_before_after":True,
            "BIO12_four_var_S2_arvense_introduced_VIF_above5":True,
            "hue_as_joint_2df_only":True,
            "corrected_families_per_model_posthoc":True}}


def table(result):
    o=["# V4 S1–S2：4変数×2 / VIFstep≤5 / VIFstep≤10 の比較",
       "","全取得形質について、補正前P<0.05、BH-FDR q<0.05、Holm P<0.05を並列比較。",
       "色相はsin/cos成分を二重計上せず、2自由度共同検定で1形質とする。",
       "","## モデルとVIF","",
       "| モデル | 説明変数 | 最大VIF |","|---|---|---:|"]
    for m,cols in result["model_predictors"].items():
        v=result["counts"]["S1/"+m]["max_source_vif"]
        o.append(f"| {m} | {', '.join(cols)} | {v:.3f} |")
    o+=["","## 未補正と補正の比較","",
        "| 層 | モデル | 検定数 | 未補正P<0.05 | BH q<0.05 | Holm P<0.05 |",
        "|---|---|---:|---:|---:|---:|"]
    for layer in ("S1","S2"):
        for model in result["model_order"]:
            x=result["counts"][layer+"/"+model]
            o.append(f"| {layer} | {model} | {x['tested']} | {x['raw']} | {x['BH']} | {x['Holm']} |")
    o+=["","## 形質群別",""]
    for layer in ("S1","S2"):
        o += [f"### {layer}","",
            "| 形質群 | モデル | 検定数 | 未補正 | BH | Holm |",
            "|---|---|---:|---:|---:|---:|"]
        for t in ("orientation","colour","outline","bract_photo_proxy"):
            for m in result["model_order"]:
                x=result["counts"][layer+"/"+m]["by_trait_group"].get(t)
                if x:
                    o.append(f"| {t} | {m} | {x['tested']} | {x['raw']} | {x['BH']} | {x['Holm']} |")
        o.append("")
    o+=["## 共通する環境変数における効果量のモデル依存性","",
       "以下は少なくとも3モデルで推定し、1モデル以上BHを通過した関連。",
       "4変数BIO12型しか持たないBIO12の係数は5/7変数モデルと直接比較できない。",""]
    def disp(v,field):
        if v is None:
            return "—"
        x=v[field]
        return "joint" if x is None else (f"{x:+.3f}" if field=="beta_std"
                                            else f"{x:.3g}")
    for study in ("S1","S2"):
        out=[x for x in result["shared_effects"][study]
             if x["n_models"]>=3 and x["n_BH_sig"]>=1 and x["sign_consistent"]]
        out.sort(key=lambda x:(-x["n_BH_sig"],-x["n_raw_sig"],
                               x["group"],x["environment"]))
        o += [f"### {study}","",
          "| 対象/形質 | 環境/スケール | モデル数 | β（BIO12/BIO18/VIF5/VIF10） | 元P | BH q |",
          "|---|---|---:|---|---|---|"]
        for x in out[:55]:
            r=x["values"]
            b=" / ".join(disp(r.get(m),"beta_std") for m in result["model_order"])
            p=" / ".join(disp(r.get(m),"p_raw") for m in result["model_order"])
            q=" / ".join(disp(r.get(m),"q_BH") for m in result["model_order"])
            o.append(f"| {x['group']}/{x['coordinate']} | "
                     f"{x['environment']}/{x['scale']} | "
                     f"{x['n_models']} | {b} | {p} | {q} |")
        if not out:o.append("該当する比較はありません。")
        o.append("")
    o += ["## 解釈上の制約","",
      "VIFは共線性の確認、P補正は多重比較の確認であり目的が異なる。",
      "Holm/BHは各層×各モデルで計算した。モデル間で最も小さいPを",
      "事後選択すること自体の誤差を補正したものではない。",
      "モデルごとに共変量と検定件数が異なるので、通過件数を性能順位としない。",
      "S1は地理・系統・画像姿勢の残存交絡があり、S2は2種の異なる分布域の写真データ。",
      "総苞片はS1の画像突出・広がりproxyのみであり刺・粘液ではない。",
      "標準化βと95%CIを含む全結果はこの成果物のJSON参照。",
      ""]
    return "\n".join(o)+"\n"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--aggregate",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    d=analyze(json.loads(args.aggregate.read_text(encoding="utf-8")))
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"V4_S1_S2_four_model_comparison.json").write_text(
        json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    (args.out/"V4_S1_S2_four_model_comparison_JA.md").write_text(
        table(d),encoding="utf-8")
    print(json.dumps({"status":d["status"],"counts":{
        k:[v["raw"],v["BH"],v["Holm"]] for k,v in d["counts"].items()}},
                     ensure_ascii=False),flush=True)

if __name__=="__main__":
    main()

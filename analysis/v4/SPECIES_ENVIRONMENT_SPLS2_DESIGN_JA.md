# V4 元S1の全8花形質×気候22変数：Sparse PLS2 共変動軸探索（2026-11-07）

**目的**：気候変数を単独に順位づけするRFやEWASとは異なり、頭花の向き、花色（明度・彩度・色相）、輪郭4構成を**同時に**説明する気候潜在軸を探索する。元のS1/S2/H3c結果を上書きしない。探索的・事後的な研究であり、独立確認ではない。

## 母集団と尺度

- 元V4暗号化S1、WCVP自生域313,350観察を既存保護GitHub ActionsでAES-GCM・SHA・元環境画素検証して復号。元データ表・撮影者・分類群名・撮影位置・予測の個別行は公開しない。
- 元の**同一写真で全8構成がQC適格かつ有限値**である行だけを種平均の元とし、分類群ごと最低5写真。**8地理セル/40写真の条件は使わない。** そのため、向きのみの178分類群や各形質別177分類群と同じNを保証しない。実Nは実行時に報告する。
- 従属形質は画像鉛直角度、L\*、C\*、単位色相sin/cos、2D縦横比・円形度・solidity・幅プロファイルCV。色相は2成分を合わせて1形質と数え、算術角度平均は使わない。色相の種平均合成長0.2未満は不定とみなす。
- 説明変数は既存の元画素**BIO01–BIO19、CHELSA RSDS、VPD、風速の22変数**。**標高・DEM地形5変数・土壌6変数は、このモデルに入れない。** 欠測気候は全22変数の完全ケースにまとめて削らず、同じ適格写真内で各気候変数ごと最低5有効値を使って種平均化。学習時点での中央値補完のみ。
- 分類群均等重み（各分類群1行）。写真数の多い種の地理的観測偏り、撮影者バイアスと種平均の測定誤差は本解析だけでは除去されない。

## 実装したSparse PLS2（R mixOmicsではない）

- Pythonの実装では、Xを各学習fold内で中央値補完→中心化・標準化。Yも学習fold内で中心化・標準化する。色相の2成分は共同の合成標準偏差で割り、円周色相**全体で他の各形質と同じ1単位**の分散重みとする（色相を2形質として過大に扱わない）。
- X–Y共分散の上位特異方向を交互更新し、X方向の重みをL1 soft-threshold（絶対値が上位のkeepX個まで非零）、Y方向は全成分使用。PLS2のX/Y回帰deflationで次成分を求める。**mixOmicsパッケージの結果の再現と呼ばない**。keepX=22では非疎なPLS2となる。
- 候補keepX：**2, 4, 7, 11, 16, 22**、成分数：**1または2**。選択は外側評価foldの**訓練分類群だけ**で4-fold内側交差検証。最小誤差、同点は少ない変数・成分を優先。外側は分類群5-fold×2反復、種平均ベースラインより予測誤差を何割減らしたかをSkillで報告。
- 各形質の予測Skill、8形質合計Skill、外側各foldの選択数、環境変数の外側選択頻度、全データで定めた構成に条件づけた100回分類群bootstrap選択頻度、X重みとY loadings を出す。符号は潜在軸の任意向き、**物理的な気候効果の正負ではない**。
- 経度帯3分割による未知地域評価では、二つの経度帯だけで内側CVを行い、残り一帯で検証する。総Skill・形質別Skill・各地域のSkillを両方報告。これは独立した3系統の反復ではない。
- 生データは保護作業ディレクトリに限る。匿名集計JSONのみartifactで公開し、フィールド名を走査して`obs_id`、`accepted_key`、`taxon`、`latitude`、`longitude`、`user_id`等がないことを確認する。

## 位置づけと結果の解釈

本sPLS2は、全8構成の**共通観察集合の種平均**に限る。以前の別々の177/178分類群RFのSkillとは、応答形質の組・平均写真の集合・多変量標準化・損失の尺度が違うので、絶対値を直接比較して一方が上位と宣言しない。必要なら同じ共通集合でRFを再適合する。

生態学的主張は「どの気候軸にどの頭花画像構成の種間差が共変動するか」という**観測的パターン**。系統・地域・観測努力を条件づけない潜在軸は、独立した気候選択や遺伝的モジュールを識別しない。全地域外Skillが負なら、群内の正の関係を全球的に一般化しない。花色は画像proxyで、色素定量ではなく、頭花画像鉛直は重力基準の下垂ではない。

## 実行方法

`analysis/v4/species_environment_spls2.py` は既存 `species_environment_existing_inputs.load_existing_inputs` を再利用する。元S1と既存7+9環境表・ハッシュ一致レシートがある保護環境で以下を実行する。コード自体には復号、外部通信、公開原データ取得機能がない。

```bash
python -m analysis.v4.species_environment_spls2 \
  --native work/rf21/private/decoded/observations.csv.gz \
  --decoded-receipt work/rf21/private/decoded/measurement_receipt.json \
  --seven work/rf21/private/environment/original_seven.csv.gz \
  --nine work/rf21/private/environment/chelsa_extra_nine.csv.gz \
  --nine-manifest work/rf21/private/environment/chelsa_nine_manifest.json \
  --source-receipt work/rf21/public/bio21_source_receipt.json \
  --source-receipt-sha256 26b80b279fa2b372099ef945972d3e4fdf58127fe9f8fd92530ac08a931d9545 \
  --output work/rf21/public/v4_original_s1_species_joint_spls2_climate22_aggregate.json \
  --bootstrap 100
```

生成物は全8構成・22気候の匿名集計のみであり、未実行のモデルを実測値として埋めない。

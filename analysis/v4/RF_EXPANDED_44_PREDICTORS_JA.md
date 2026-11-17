# V4独立探索PR: Random Forestの候補を12→44環境変数に拡張

**2026-10-31 状態:** ソース取得レシピ、登録表、段階的RF推定器、QAと合成データテストを実装済み。**元の非公開S1観測における追加32層のラスタ抽出／モデル適合／予測値は未実行**。本書は新しい肯定的な実測結果ではない。

## 研究目的と凍結境界

新しい問いは「**頭花の向き、花色、2D輪郭の局所的な画像形態が、どの環境領域を追加することで未知の分類群・地域にも予測可能になるか**」。既存論文の3.1降水と画像上の頭花角度、3.2二種の在来/導入域における花色、3.3気候類似地点の多形質類似性、原S1/S2のHolm結果・H3c探索的家族は変更しない。生物学的なメカニズムを同定したとは主張しない。

RFの入力は環境説明変数だけ。**分類群名、観察ID、経緯度、国、撮影者、元の花色・頭花形態から派生した変数は説明変数にしない**。これらはサンプル認証・グループ分割・観測単位の指定だけに用いる。

## 段階的に増やす変数群

| 比較の右側のモデル名 | 追加するグループ | 累計変数数 | 出典と定義 |
| --- | --- | ---: | --- |
| `archived12` | 原S1のCHELSA BIO01/04/05/06/12/13/14/15/17/18、RSDSとWorldClim標高 | **12** | 元からソース確認済み。ここまでの4→12は既存RFパイロットのまま |
| `complete_bioclim21` | **BIO02/03/07/08/09/10/11/16/19** | **21** | CHELSA-bioclim v2.1、1981–2010。追加9層は今回まだ取得していない |
| `water_balance26` | VPD、PET、AET、気候水分不足、土壌水分 | **26** | TerraClimateの月値から1981–2010平均/平均年合計を再現して作成 |
| `wind28` | 平均風速・季節内風速変動 | **28** | 同じTerraClimateの1981–2010月値、再現可能な集約が必要 |
| `terrain33` | 傾斜、起伏、TPI、northness、eastness | **33** | 元と同じWorldClim 2.1 DEMから固定近傍サイズで派生レイヤを作成 |
| `soil39` | 土壌pH、SOC、粘土、砂、bulk density、CEC | **39** | SoilGrids v2.0、**全て深さ0–5cm、Q0.50**で統一。原格納値→物理単位変換必須 |
| `vegetation44` | NDVI平均、NDVI月季節変動、森林/農地/草地の被覆率 | **44** | MOD13A3 v061の2001–2010品質管理済み月値＋Copernicus LC100 Collection 3の2015年被覆率 |

論文の現行4変数（BIO18、BIO05、RSDS、標高）と既存12変数は完全に別名のモデルとして残す。**BIO03=100×BIO02/BIO07など、一部のBIOCLIM指標は他列から決定的に導出される**。RFに投入しても独立した因果情報が44個あるわけではない。値が足りないときの追加変数補完、無根拠な値の推定、観察日の天候との混同は禁止。

元の写真に近い年代を一貫して復元できない環境層は、**時間尺度が異なる感度解析**として扱う。NDVI 2001–2010と土地被覆2015を「1981–2010の気候」や「撮影当日の植生」と呼ばない。SoilGridsは現地土壌実測ではなく予測表層であり、作成時の気候・地形等が共変量になっている。雲・影・色補正・画像角度は未校正の写真測定限界のまま。

### 外部データの一次資料

- [CHELSA-bioclim（全19 BIOCLIM、v2.1）](https://www.chelsa-climate.org/datasets/chelsa_bioclim)
- [TerraClimate variable definitions（VPD、PET、AET、DEF、soil moisture、wind）](https://www.climatologylab.org/terraclimate-variables.html)
- [ISRIC SoilGrids 250m](https://isric.org/explore/soilgrids)・[depth/scaling](https://docs.isric.org/globaldata/soilgrids/SoilGrids_faqs_01.html)
- [WorldClim v2.1 original DEM used by the project](https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip)
- [NASA MOD13A3 v061 monthly NDVI](https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD13A3)
- [Copernicus LC100 Collection 3（2015–2019）](https://land.copernicus.eu/en/products/global-dynamic-land-cover/copernicus-global-land-service-land-cover-100m-collection-3-epoch-2015-2019-globe)

これらは**入手可能な公式データ製品の資料**であって、今回の新規32変数が既に元S1へ抽出された証拠ではない。

## 同じデータ同士を比較する仕組み

一度の解析では隣接する**2モデルのみ**を同一の形質別・環境完全ケース集合に当てる。例: 12対21、21対26、26対28、28対33、33対39、39対44。比較ペアの右側で必要な説明変数だけ完全ケースを要求し、他の未導入の層に欠測があっても右側モデルの候補集団から不必要に落とさない。欠測により分類群が8セル/40観察のゲートを失った場合は分類群ごと除外し、どちらのモデルでも復活させない。**異なるペアの絶対スキルを並べて「変数数の効果」と解釈しない**。

- 観測: 生写真の8形質（9数値座標）と環境を統合し、各形質で同じ分類群×5°地理セルの写真中央値にする。色相はsin/cos共同2次元の円形誤差。
- 検証: 未知の分類群に対するGroupKFoldと、地理組成に感度が高い**未知の3経度セクター**Leave-One-Group-Out。後者は同じ分類群が訓練・評価に存在し得るため、未知分類群テストと同義ではない。広域テストは1セクター5分類群以上を要求する。
- 指標: 学習fold平均に対するheld-outの重み付き二乗スキル、右モデル−左モデルのペア差、**環境領域をまとめた**held-out permutation損失。全評価で分類群均等重みと凍結ハイパーパラメータを保つ。新しいP値や「最重要の独立因子」は発表しない。
- 注意: 相関した説明変数群の置換は無現実な組合せを作り得る。変数重要度は因果・媒介効果・最適適応の証明ではない。重要度が低いことは生態学的に無関係という意味でもない。

## 実行できる状態と外部依存

新しいコード:

- `analysis/v4/rf_extended_environment.py`：44軸登録、元12軸の保全、層/観測ID/マニフェストのチェック、隣接ティアのRF解析、公開可能な集計JSON
- `analysis/v4/rf_extended_raster_sampler.py`：保護された元S1座標に対して**事前に作成された32候補ラスタから指定ティア分のみ**ピクセル抽出。ファイルSHAを検証し、無効値をNaNのまま残す。外部データの自動ダウンロードなし
- `analysis/v4/rf_multitrait_predictability.py`：既存4対12を変更せず、ペアごとに要求変数集合/変数群を渡せる後方互換の拡張
- `tests/test_v4_rf_extended_environment.py`、`tests/test_v4_rf_extended_raster_sampler.py`：合成テストのみ

作業例（**ソース未取得のため実行済みではない**）:

```bash
# 1. まず追加のBIO9だけを用意するレシピを生成。空欄は未取得で実行拒否される。
python -m analysis.v4.rf_extended_raster_sampler \
  --upper-tier complete_bioclim21 \
  --create-recipe-template /secure/bio21_recipe.json

# 2. 公式CHELSAから9層を取得し、実ファイルSHA、物理単位変換、
#    元の測定の解像度・期間をレシピに記録したあと、保護環境内でのみ抽出。
python -m analysis.v4.rf_extended_raster_sampler \
  --upper-tier complete_bioclim21 \
  --native-table /secure/original_s1.csv \
  --native-receipt /secure/original_s1_receipt.json \
  --recipe /secure/bio21_recipe.json \
  --output-table /secure/bio21_extra9.csv \
  --output-manifest /secure/bio21_manifest.json

# 3. 元S1の既存7追加環境CSVと合わせ、12対21を同一観測で評価。
python -m analysis.v4.rf_extended_environment \
  --upper-tier complete_bioclim21 \
  --native-table /secure/original_s1.csv \
  --native-receipt /secure/original_s1_receipt.json \
  --archived-seven-table /secure/original_extra7.csv \
  --archived-seven-sha256 ACTUAL_SEVEN_SHA256 \
  --extra32-table /secure/bio21_extra9.csv \
  --extra32-sha256 ACTUAL_NINE_SHA256 \
  --layer-manifest /secure/bio21_manifest.json \
  --output /secure/aggregates/rf_12_vs_21.json
```

`extra32-table` は歴史的な引数名で、**そのティアまでに追加した列だけ**を含む（21変数ティアでは9列）。ティア拡張時には既存レシピの対象層も再収録する。マニフェストには取得物のURL、時期、単位、導出規則、ラスタ単位SHA、表全体SHAを記載。ラスタ実体の一意性は専用サンプラーの実行中に検証され、RF側は受け取ったマニフェストを確認する。**原写真、座標、写真ID、種ID、予測行、原ラスタをPRや公開artifactにアップしない。**

**次の科学的判定**: 12→21の未知分類群・未知地域の予測差が成立するかを先に評価し、追加領域（水収支、風、地形、土壌、植生）の効果は同じゲートで順次確認する。支持が出ても探索的SI候補として別判断するまで、今のGEB本文・図・係数を触らない。
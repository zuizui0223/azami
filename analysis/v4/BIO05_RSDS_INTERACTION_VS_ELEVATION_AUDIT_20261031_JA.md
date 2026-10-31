# V4 3.1・3.2：標高調整とBIO05×RSDS交互作用を比較する（2026-10-31）

**優先：科学的な結果の取捨選択ではなく、どの予測を検証するかで説明変数を決める。** ３変数でHolm支持が増えたという事実を、解釈したくないから削除してはいけない。これらの結果が**標高の調整、写真コホートの違い、調整する交互作用によってどう変わるか**を、元の同じ写真で検証する。

## モデル（3.1と3.2の同じ元S1/S2写真で比較）

| モデル | 説明変数 | 生態学的な役割 |
| --- | --- | --- |
| **M3** | BIO18・BIO05・RSDS | 標高を外すと現れた関連の診断。現在の論文主解析には選ばない |
| **M4 元の凍結解析** | **BIO18・BIO05・RSDS・標高** | 今の主解析の**第一候補**。標高による地理・地形環境の調整を残す |
| **M3I** | BIO18・BIO05・RSDS・**BIO05×RSDS** | **標高の代わりに交互作用を入れる案**。M4と４モデル項の同数で比較可能 |
| **M4I** | BIO18・BIO05・RSDS・標高・**BIO05×RSDS** | 標高を調整したまま、気温が日射の関係を修飾するかを検証 |

**3.3は現行のBIO12・BIO01・RSDSのみ**（1°＋0.04096366、5°＋0.09433109、元探索的Holm-2 P 0.010）を継続し、モデル変更しない。

交互作用を加える理由は**「気温により日射と頭花形態の関係が違うか」**であり、標高が気候の代理をする問題への解決ではない。BIO05とRSDSはそれぞれ**主効果を残す**。S1は**分類群間の温度×日射（taxon mean×taxon mean）**の１項のみを追加。S1の分類群内交互作用を検証したものと混同しない。S2は、２種×在来・導入４群それぞれにBIO05×RSDSの１項を加える。

## なぜ３変数では別の相関が増えたのか

まず交絡と多重補正を区別する。前の旧M4（標高あり）は形質10×環境４×分類群内／間の**80件**がHolm対象、新３変数は形質８×環境３×分類群内／間の**48件**。追加の写真（旧S1向き168,907、新S1向き169,952）も違う。**P値の変化だけを標高の効果とは読めない**。

同じ旧全12環境完備コホートで、M3 vs M4（標高だけ）、M3 vs M3I（交互作用だけ）、M4 vs M4I（標高調整後の交互作用）、M3I vs M4（別モデルの適合）を比較する。また、同じ３変数M3を旧全12環境完備と新６気候完備の２コホートで再推定して、**写真採択の差**を別項目にする。写真が異なる２コホート間でAICを直接比較しない。

**原S1 BIO18×画像角度 β +0.045446843、旧168,907／178分類群、S2元６係数と対象数**を元データから厳密再現することを実行ゲートとする。形質は画像向き、明度、彩度、色相（sin/cos２自由度共同検定）、縦横比、円形度、solidity、幅変動の８構成。２次元色相を２つの独立な検定と数えない。

## 交互作用の評価基準

- **固定したモデル・固定した対象写真・固定したWLS重み・クラスタ頑健SE**で適合比較。S1分類群間の１項／８形質の**新８検定**、S2各群の１項／４群×８形質の**新32検定**に個別Holm補正をする。M3IとM4Iはそれぞれ別の探索的ファミリー。
- S1/S2の「元M4」と「M4I」、および「M3」と「M3I」で**対応するモデル間AIC差・残差減少率**を記録。非入れ子のM3I対M4についても同じコホートのAICは参考値になるが、調整した因果関係を識別しない。
- BIO05×RSDSが支持されても、**温暖化による花の光防御、花粉損傷、種の適応的応答を実証したわけではない**。元写真の色と向き、地域・系統構成は依然として観測上の条件付き相関。
- **原2026-10-24のS1 80／S2 128検定とその公開用集計は来歴として不変**。新たに得る８・32個の相互作用検定は後付けの探索的感度分析で、主張を置換しない。

## 実装と次の判断

[コード](s1_s2_bio05_rsds_interaction_source.py)／[合成データテスト](../../tests/test_v4_s1_s2_bio05_rsds_interaction.py)／[元データ専用Actions](../../.github/workflows/v4-s1-s2-temperature-radiation-interaction.yml)。

**現在地：** 今回は元の４変数を主解析の第一候補に戻し、交互作用は科学的に妥当な追加仮説として定義した。**交互作用の元データ検定値は、このワークフローが完走するまで未確定。** 後付けの有意性から主モデルを変更せず、地域・分類群を考慮した追加検証で意思決定する。Results 3.3の年間気候類似性結果は変更しない。


## Source-verification outcome (31 October 2026)

The implementation passed original-source tests in [GitHub Actions #41182224216](https://github.com/zuizui0223/azami/actions/runs/41182224216), with public aggregated-only artifact `13035453002`. See [the full source-verified scientific report](results/20261031_S1_S2_BIO05_RSDS_INTERACTION_SOURCE_VERIFIED_JA.md) and [machine summary](results/20261031_S1_S2_BIO05_RSDS_INTERACTION_SUMMARY.json). **S1: 0/8 interactions after Holm in both M3I and M4I. S2: 1/32 after Holm for introduced *C. vulgare* chroma (M4I β−0.088454, Holm P=0.011458).** Source-matched S1 radiation×chroma coefficients show that elevation removal, not BIO05×RSDS itself, is responsible for the major coefficient difference. Therefore **retain M4 original as main**, and present M4I as exploratory sensitivity. Original H3c is unchanged; original four-variable manuscript/figures remain correctly matched to main results.

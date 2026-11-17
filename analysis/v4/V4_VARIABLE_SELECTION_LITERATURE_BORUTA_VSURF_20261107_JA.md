# V4変数選択の先行研究監査：Boruta・VSURF・Sparse PLS2（2026-11-07）

**位置づけ：査読済み方法論・生態学的適用例の文献確認に基づく選択プロトコル。Boruta/VSURFによるV4元データの変数選択結果は**未実行**。本稿のS1 80検定・S2 128検定・H3cは変更しない。研究対象は分布予測（SDM）ではなくアザミ属の種間／種内の画像花形質であり、先行SDM研究の性能値は転用しない。**

## 原著で確認できること

1. **Boruta: all-relevant selection.** Kursa, M.B. & Rudnicki, W.R. (2010), *Feature Selection with the Boruta Package*, Journal of Statistical Software 36(11), 1–13, DOI: https://doi.org/10.18637/jss.v036.i11 。元変数の重要度を無作為にシャッフルしたshadow変数と繰り返し比較してConfirmed／Rejected／Tentativeを判定。大きく相関した説明変数がまとめて関連候補になる可能性がある。「1変数だけを選ぶ最小集合」とは別の目的。2026年のRパッケージBoruta 10.0.0は数値応答（回帰）も扱える：https://cran.r-project.org/web/packages/Boruta/index.html 、https://cran.r-project.org/web/packages/Boruta/refman/Boruta.html 。
2. **VSURF: interpretationとpredictionを区別。** Genuer, R., Poggi, J.-M. & Tuleau-Malot, C. (2015), *VSURF: An R Package for Variable Selection Using Random Forests*, The R Journal 7(2), DOI: https://doi.org/10.32614/RJ-2015-018 。(a) RF重要度の安定した順位と低重要度の除外、(b)冗長性をある程度許す解釈用集合、(c)予測力をほぼ保ち冗長性を減らした予測用集合、を出す。著者は、選択工程中のOOB誤差を最終モデルの独立した予測誤差と呼ばず、別のCVが必要と明記。
3. **Sparse PLS: 多応答の共変動による選択。** Lê Cao, K.-A. et al. (2008), *A sparse PLS for variable selection when integrating omics data*, Statistical Applications in Genetics and Molecular Biology, https://pubmed.ncbi.nlm.nih.gov/19049491/ 。PLSに疎化を加えて2ブロックの共変動と変数選択を併用。Rohart, F. et al. (2017), *mixOmics*, PLoS Computational Biology 13:e1005752, https://doi.org/10.1371/journal.pcbi.1005752 はkeepXの交差検証と選択安定性を強調。V4のPython実装はmixOmicsソフトそのものではない。
4. **相関変数の重要度は単独の因果効果でない。** Strobl, C. et al. (2008), *Conditional variable importance for random forests*, BMC Bioinformatics 9:307, https://doi.org/10.1186/1471-2105-9-307 。通常のRF重要度は説明変数の相関で歪み、条件付きPermutationが一部を改善。Dormann, C.F. et al. (2013), *Collinearity: a review of methods to deal with it*, Ecography 36:27–46, https://doi.org/10.1111/j.1600-0587.2012.07348.x は生態学的な近接性と機構に基づく変数の選択を重視。VIFなどの統計のみで「正しい原因」は確定しない。
5. **空間／階層CVが必要。** Roberts, D.R. et al. (2017), *Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure*, Ecography 40:913–929, https://doi.org/10.1111/ecog.02881 。種・地域・系統の依存性を無視した検証は楽観的になりうる。Meyer, H. et al. (2019), *Importance of spatial predictor variable selection in machine learning applications*, Ecological Modelling 411:108815, https://doi.org/10.1016/j.ecolmodel.2019.108815 も、地理的位置情報による過剰適合と空間保留の重要性を実例で示す。

## 実際の生態学的適用研究と、相反する教訓

- Hanberry, B.B. (2024), *Practical guide for retaining correlated climate variables and unthinned samples in species distribution modeling, using random forests*, Ecological Informatics 79:102406, https://doi.org/10.1016/j.ecoinf.2023.102406 。北米80樹種のSDMでは温度7・降水7の計14気候変数を保持したRFが、2変数への事前縮約より高精度。**「相関を避けるために温度1・雨1へ必ず絞る」が正解ではない**ことを示す一事例。ただし本V4は種間の連続花形質を予測しており、性能結果を直接移植できない。
- Cushman, S.A. et al. (2024), *Simulating multi-scale optimization and variable selection in species distribution modeling*, Ecological Informatics 83:102832, https://doi.org/10.1016/j.ecoinf.2024.102832 。Boruta・MRMR・VIF等を比較するシミュレーションで、MRMR／univariate scaling＋GLMのモデル群が総合的に良好で、Boruta単独が万能なわけではなかった。**試す手法を一つに固定せず、同一外部評価で比較する**根拠となる。これもSDMでありアザミ花形質の実測ではない。
- 厳密な生態学変数選択は「結果が有意な変数を後から採用する」ことと別。本文確定済みH3c（同一S1の花形質を使う）は、同じS1全体で花形質を使って選んだ軸・変数で再検定すると、選択バイアスを含む。独立の評価分類群または選択自体を帰無反復に含めることが必要。

## V4の元S1で得られている状況（確認済み実測、Boruta未実行）

- 元在来S1 313,350観測だが、モデルの**独立な解析行は約173–178分類群**であり313,350独立反復ではない。RF（標高なし気候22）:色彩度の上位にRSDS・BIO08、色相にBIO08、向きでは個別1位BIO07だが**降水変数群をまとめた重要度**が大きい。元のBIO18×向きはHolm陽性だが強い地理調整で減衰する。
- Sparse PLS2:全8形質173分類群、成分1にBIO08／RSDS／BIO10／BIO01、成分2にBIO13／BIO16／BIO18／BIO12。外側10fold中6foldでkeepX=22が選ばれ、**全変数選択の頻度が高いことを疎で安定した選択の証拠とはしない**。未知経度帯の全8形質予測Skillは負。
- 2026-11-07の3候補に限るRF・PLS置換P（3候補内Holm）は**事後探索**であり、22変数×8形質の全選択過程を補正した確認的なP値ではない。
- 元3.1回帰の変数はBIO18・BIO05・RSDS・標高。元3.2も同じ4説明変数で両種×在来／移入の観察群を比較し、導入域C. arvense標高×明度のHolm陽性がある。元3.3気候アナログはBIO12・BIO01・RSDS。**この元結果は変更しない**。S2のBIO08・他追加気候軸の全観測レイヤ一致は未確認。

## V4の採用判断：二段階で変数数を選ぶ

**目的A（関連候補の発見）**: Boruta all-relevant、sPLS2の負荷量・条件付き再標本化、RF気候ファミリー重要度を併記し、温度・降水・日射・湿度／風ごとに関連候補集合を報告。BorutaではConfirmed/Rejected/Tentativeを残し、shadow比較で残った数を事前に3・4と決めない。Boruta 8形質は基本7つのスカラー応答で独立に実行し、**色相sin/cosの2成分は円周1形質として同時評価する多応答shadowの設計が必要**。単純なsin,cos各別の重要度を1形質の有意性と混同しない。

**目的B（本論文3.1–3.3の簡潔な共通変数）**: VSURF interpretation setとprediction setを分ける。変数の数を固定せず、同一元分類群・学習foldで3・4・5・6・8・22変数またはVSURF選択集合を評価し、学習側のみで特徴選択。外側の未知分類群および地理ブロック保留でSkill・各形質Skill・VIF・環境サポートを比較。小モデルは最良モデルと実質的に同等（foldの不確実性を考慮）なら選ぶ。科学的に本来必要なBIO18などの降水軸を「向きRFの全体Skillが小さいから」無条件に削除しない。

**実行ゲート（固定結果の保持）**:
1. Boruta＋VSURFの結果と各変数の選択頻度を元S1の保護runnerで実測。B​​orutaの選択ステータス・VSURFの2集合・気候クラスター別の評価を匿名集計で出す。**このファイルだけではBorutaのConfirmed数や新たなP値は未確定**。
2. 「気候3変数BIO08/BIO18/RSDS」は**比較候補**の1つに据える。元のBIO05/BIO18/RSDS＋標高の4変数は基準、標高なし気候3→+標高を同じ標本の感度とする。標高を純粋な気候軸の一つとしてBoruta選択に混ぜない。
3. 3.2の2種×2域の全新気候レイヤの一致を事前確認し、不完全なら新候補をそのまま移入域へ転移したと主張しない。
4. 3.3は元のBIO12/BIO01/RSDS気候アナログ肯定結果を保持。花形質で新選択した気候変数で同じS1地点を使う追加検定は、事後探索／評価分離として報告し、既存Holm Pと同一帰無と扱わない。

**現時点の判断：変数数は未決定。Boruta単独ではなく「all-relevant発見→冗長性の評価→VSURF等による簡素化→外部CV」を採用し、3.1–3.3の元解析は維持。**

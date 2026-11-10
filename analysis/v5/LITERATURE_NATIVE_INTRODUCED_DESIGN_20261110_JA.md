# V5：在来域・移入域を比較する理由と先行研究準拠の検証設計（2026-11-10）

> **状況**：既存V5の計算を「正になったものだけ」を主結果へ移すのではなく、侵入生態学の公開済み研究で使われている *native–introduced comparison* の問い・利用可能環境の調整・対立仮説に合わせて整理する。下記の新規背景補正地域ペア検定は**提案仕様であり未実行**。既存ソースやV4の凍結結果は変更しない。ユーザー指定により気候**希少度の連続指標を現行V5の主解析に使わない**。

## 1. なぜ在来域と移入域を比較するのか

「在来域→移入域で適応進化したこと」を確認するためではない。同じ種が人為的に地理的障壁を越えたとき、(a) *占める気候*、(b) *頭花表現型の多変量分布*、(c) *気候と頭花の対応* が維持されるか、どの成分が変わるかを、それぞれ独立に検証するためである。

在来域は比較の参照分布だが、移入時の直接の祖先集団ではない。在来・移入の違いは、利用可能な気候、導入元・創始者、敵・送粉者、可塑性、地域内の進化、写真・観察の偏りから生じ得る。遺伝的適応や新しい基本ニッチを証明する「自然実験」ではなく、**歴史的に異なる分布域を利用した観察的・同種内比較**とする。2種は独立した大量の導入実験の代替ではない。

**Part Iとの論理接続**：Part Iは「異なる分類群・地理の気候類似性と頭花類似性が対応するか」。Part IIは「分布史を変えても、*同じ種*のその対応が保たれるか」。前者だけでは種間の歴史と生態が混じる。後者は同じ種の比較である一方、地域・導入史の交絡が残る。両者が補完するのが採用理由である。

## 2. 先行研究の方法をどう継承するか

| 原著（確認した要点） | 論文での実際の比較 | V5で採用する部分／できない部分 |
|---|---|---|
| Hulme & Barrett (2013), *J. Ecol.* [DOI:10.1111/1365-2745.12009](https://doi.org/10.1111/1365-2745.12009) | 在来／移入の**環境範囲・形質平均・形質～環境の傾き**を分離し、range×environment交互作用と複数のシナリオを提案。開花日の温度関係を例示 | 主たる理論テンプレート。同じ気候でも形質平均・傾きが違うか、新規気候でも既存形質分布に収まるかを分ける。V5の写真から遺伝・可塑性は特定できない |
| Broennimann et al. (2012), *Global Ecol. Biogeogr.* [DOI:10.1111/j.1466-8238.2011.00698.x](https://doi.org/10.1111/j.1466-8238.2011.00698.x) | 発生地点と背景環境から共通PCA環境空間を作り、出現密度と利用可能環境の偏りを考慮しDとニッチ類似性を評価 | 全域の背景補正済み気候Dを継承。観測写真だけの気候Dと明示的に別にする。V5の実装はR ecospatと**bitwise同一ではない** |
| Guisan et al. (2014), *Trends Ecol. Evol.* [DOI:10.1016/j.tree.2014.02.009](https://doi.org/10.1016/j.tree.2014.02.009) | 類似した利用可能気候（analog background）で **stability, expansion, unfilling** を分解。非類似気候を別扱い | 共有背景の占有拡大候補と、そもそも気候が一方にしか利用可能でない場合を分ける。「95%域外＝新しい適応ニッチ」と呼ばない |
| Petitpierre et al. (2012), *Science* [DOI:10.1126/science.1215933](https://doi.org/10.1126/science.1215933) | 50植物種の在来・移入比較。analog climateに限ると大きな気候ニッチ転換は少ないと報告 | 共有気候への制限が結論を変える可能性を事前に認識。域外割合を全移入密度分母と共有背景分母で混同しない |
| Atwater et al. (2018), *Nat. Ecol. Evol.* [DOI:10.1038/s41559-017-0396-z](https://doi.org/10.1038/s41559-017-0396-z) | 815植物種・5大陸、ニッチ変化は指標と地域の気候利用可能性に依存。大陸間の変化と在来域内の比較も考察 | 在来の複数地域を外した「自然な地域差」を基準にする。50種研究と対立する割合をそのままV5の事前期待にしない |
| Rönnfeldt et al. (2026), *Nature Commun.* [DOI:10.1038/s41467-025-68023-6](https://doi.org/10.1038/s41467-025-68023-6) | **316植物種×1,566移入事例、8地域**を分けて同種の地域依存性・滞在期間・背景気候を比較。類似背景でニッチ変化を評価 | 「全球の移入域ひとまとめ」を唯一の比較にしない。可能な限り移入地域ごとに検証し、20°地理ブロックは独立の導入イベントではないことを明記 |
| *Rumex* 3種の複数移入地域比較 (2023, *Biological Invasions*) [DOI:10.1007/s10530-022-02893-5](https://doi.org/10.1007/s10530-022-02893-5) | 同種でも導入地域ごとにニッチシフト方向が異なる | 地域をプールしたDだけから一律の分化・保守性を言わない |

**同じ植物群を扱う実証的な先例もある。** Guggisberg et al. (2012) [J. Biogeogr., DOI:10.1111/j.1365-2699.2012.02746.x](https://doi.org/10.1111/j.1365-2699.2012.02746.x) は *C. arvense* の欧州と北米の58集団・1,522個体に遺伝データと気候比較を用い、複数導入・祖先集団の非一様性と、比較的保存された気候ニッチを示した。したがってV5は**在来全域＝実際の導入元**として扱えない。Mráz et al. (2011) [*International Journal of Plant Sciences*, DOI:10.1086/658151](https://doi.org/10.1086/658151) はキク科 *Centaurea stoebe* の在来・移入由来の植物を同一温室条件で育てて40以上の形態形質を多変量比較し、頭花形状も含めて倍数型の影響と移入域の差を区別した。V5の写真だけの観察的多形質空間比較が、**遺伝的分化を検証する共通圃場と異なる**ことを具体的に示す先例である。Hodgins et al. (2020) [*Plant Communications*, DOI:10.1016/j.xplc.2020.100116](https://doi.org/10.1016/j.xplc.2020.100116) も *C. arvense* の在来・移入集団で共通条件とストレス処理を比較したが、V5の可視頭花C*の機能を検証した研究ではない。

上表は**先行研究の方法の継承点**であり、これらの原著が「頭花3形質Dと気候Dの地域ペア相関」をすでに立証したという意味ではない。**気候空間と写真ベース頭花空間を同じセルに結びつける部分はV5側で新たに検証する仮説**である。

## 3. 検証の主系列（気候希少度は使わない）

### 3.1 気候ニッチの変化を成立させる（Broennimann・Guisan型）

対象は *C. arvense* と *C. vulgare*。各種のWCVPによる在来・移入の野外記録を種ごとに比較し、BIO01・BIO04・BIO12・BIO15・RSDSの**5気候軸**を用いる。標高は主モデルに入れない。

- 在来・移入のWCVP土地背景（面積重み）で共通PCA2軸を作り、出現密度を各分布域の気候利用可能性で補正。
- Dを出し、*共通背景*でのstability／expansion／unfillingと、共通背景外の非類似な利用機会を別区分で記述する。占有95%・背景99%を主閾値、90/99%・背景感度を記録。
- 「在来主要95%占有外」の占有拡大候補と「在来背景にもない気候」を混同しない。WCVPの背景は現実のアクセス可能性・移入経路を完全には表さない。
- **既存実測**：背景補正気候D=0.474／0.442；在来主要95%域外（全移入密度分母）18.8%／44.0%；共通背景分母12.6%／20.4%。別分母なので無理に一致させない。
- 「当該種で移入域への到達可能性が変わった」ことと「生理学的ニッチが進化した」ことを分ける。

### 3.2 頭花形質空間の変化を同じ分布史で比較（Hulme & Barrettの形質視点を拡張）

同一の頭花適格記録の**観察セル**を在来・移入で揃え、(1) 在来 vs 移入全域、(2) 在来 vs 移入の在来類似気候、(3) 在来 vs 移入の在来主要95%外の新規気候、の頭花Dを比較。**向き・CIELAB彩度C*・画像輪郭指数**は同じ写真・セルの観測値だが、輪郭指数には回転依存性があるため向き＋C*の2次元比較を重要感度とする。Dは平均差・色素量・遺伝率ではない。

- 地域・セル数を揃えたwithin-20° region 1:1マッチを主たる頑健性評価とし、無調整の全セルDを記述値として分離。
- **既存実測**：3形質D（類似→新規）の地理マッチ後の差は +0.005（Holm 0.953）／+0.057（Holm 0.096）。帰無支持がないことは「同じである」証明ではない。
- 平均形質差と傾き差を必要な場合のみ追加。移入によって形質が進化したという断定はしない。

### 3.3 気候ニッチDと頭花Dが連動するか（V5の**新しい問い**）

**同一種の在来地域×移入地域ペア**を解析単位とし、同一写真適格セルで気候Dと頭花Dを並べる。元のV5では20°地域ごと8セルを32回抽出し、17×8、11×7の地域ペアを得ている。しかし一つの地域は多数のペアに現れるので、136ペアや77ペアを独立サンプルとは扱わず、行・列の地域ラベル置換を行う。

- **元の実測結果**：rho=+0.071／−0.108、2種Holm P=1.0／1.0。これは撮影地点の**出現気候D**と頭花Dであり、地域の**利用可能気候で補正したニッチD**と頭花Dの結合の検証ではない。
- **次の必須拡張（未実行）**：同じ地域ペアの利用可能背景気候（地域内のWCVP背景と共通PCA）を取得し、背景補正した地域別気候Dを推定。同じ観測セル・写真QC・地域の反復数・頭花Dを保ったまま、行列置換・leave-one-region-out・平滑化幅0.35/0.5/0.75で rho の頑健性を検定する。
- 地域の背景が不十分な場合は**NOT_ESTIMABLE**として止める。背景補正後の数値を、現在のrhoで代用しない。
- 同所性の補助検定は、気候Dと頭花Dの**地域ペア相関**とは別の推定対象。既存の気候非重複と頭花非重複の重み差は−0.8［−1.9,+1.1］／−1.2［−2.4,+1.6］ポイントで、明瞭な同所集中を示していない。

### 3.4 在来域内の自然な地域差を超えるか（Atwater・Rönnfeldt型の対照）

在来の1地域を保留して残りの在来地域でDを計算し、同数・同じ空間分解能の移入地域のDと比較。対象により観測の季節・写真QC・気候構成が異なることを明記。

- **既存記述値**：頭花全体の地域ずれ（1−D）中央値は *arvense* で在来0.390・移入0.363、*vulgare*で在来0.364・移入0.367。*vulgare* の彩度ずれは在来0.201・移入0.233。ただし地域ベース記述値で、適応の統計検定ではない。
- **解釈**：移入側のずれが在来地域差の範囲内なら「既存の表現型分布で新しい気候へ到達した」可能性（ecological fitting）を検討できる。これは遺伝的変化がないことの証明ではない。

## 4. 事前に整理すべき4つの対立パターン

| 気候占有の類似性 | 頭花空間の類似性 | 解釈候補（因果ではない） |
|---|---|---|
| 高い | 高い | 気候と形質の両方が保たれた |
| 低い | 高い | 違う気候でも既存の頭花分布が維持された・ecological fittingの候補 |
| 高い | 低い | 類似気候でも表現型の構成が変わった・気候以外の要因／由来の候補 |
| 低い | 低い | 両方変化したが、変化を同じ原因へ帰せない |

この表は**概念的シナリオ**であり、気候D=0.4と頭花D=0.8という値の閾値を直接比較して種を分類しない。異なる次元・背景補正・写真の選択という基準の違いを押さえ、地域ペア相関と帰無・感度分析で検証する。

## 5. Figure・Methods・Discussion の採用仕様

1. **Fig II-A**：在来／移入の**利用可能な気候背景**と出現占有密度を同じPCAに描き、類似背景・非類似背景、占有95%を区別。D・全移入密度分母／共通背景分母の値を混同しない。
2. **Fig II-B**：同じ観測セルの在来・移入**頭花多形質空間**（角度、彩度、輪郭）と向き＋彩度の感度。Dは実測3DかPCA2軸かを必ず表示。
3. **Fig II-C**：地域ごとの**気候D vs 頭花D**（現在は観測地点ベース。背景補正版は未実行）。2種別、地域再利用を考慮したrho・置換P、n地域を付記。単にDを差し引かない。
4. **Fig II-D / SI**：移入の在来類似気候 vs 在来主要95%外の新規気候で、在来頭花とのDの差。20°同地域マッチと2種Holmを優先表示。彩度などは内訳。
5. **Methodsに宣言**：上記は既存V5結果を見た後の探索的再構成。撮影者・セル・地域・分類群・2種を混ぜて有効サンプル数を水増ししない。帰無の独立単位は地理地域であって導入イベントとは限らない。

## 6. 投稿で許される結論と残る実験

**支持できる範囲**：2種の移入域が在来と異なる実現気候分布を占める。写真由来の頭花多形質分布は大きく重なり、既存の地域ペア比較から気候Dと頭花Dの一貫した正相関は検出されない。新規気候の頭花D低下は地理マッチ後、補正有意には至らない。これは先行研究が問題にした**ニッチ変化と形質の変化の非自明な対応**を示す、制約の明確なケーススタディになる。

**支持できない**：ニッチ変化と形質変化が完全に独立、適応的放散、花色の乾燥適応、新規遺伝変異、移入後の進化速度、二種から植物全体への一般化。「non-significant」は「same／no relationship」ではない。

**将来の実証**：移入元を遺伝的に同定し、同じ由来の在来・移入集団で共通圃場・相互移植・反射スペクトル・花粉生存／結実を測定し、遺伝分化／可塑性／選択を区別する（Hulme & Barrett 2013に沿う）。

### 既存V5の参照

- [現行V5 Part I／Part II](results/RESULTS_JA.md)
- [元の地域ペアの気候・頭花D](results/paired_overlap/METHODS_JA.md)
- [全域の頭花Dと気候D](results/whole_range_overlap/RESULTS_JA.md)
- [類似／新規気候での頭花D・地域置換](results/novel_head_overlap/RESULTS_JA.md)
- [在来地域を保留した基準](results/invasion_design/RESULTS_JA.md)
- [現行の主線変更と数値の帰属](results/CLIMATE_HEAD_OVERLAP_DECISION_20261109_JA.md)


## 7. 論文本文用の英文草案（原著から内容を学んで独自に執筆。投稿版では未確定）

**Introduction / rationale (working draft):**

> Comparisons between native and introduced ranges offer a within-species opportunity to examine whether climatic occupancy and phenotype–environment associations persist across major biogeographic changes. Studies of plant invasions have shown that climatic niche shifts depend on environmental availability, geographic context and introduction history, and that differences between ranges do not necessarily indicate evolved changes in climatic tolerance (Petitpierre et al., 2012; Guisan et al., 2014; Atwater et al., 2018; Rönnfeldt et al., 2026). Likewise, changes in phenotypic traits need not parallel changes in occupied climates: pre-existing variation, environmental filtering, plasticity and nonrandom introduction of source populations can generate different combinations of niche and trait dynamics (Hulme & Barrett, 2013). We therefore treated the native ranges of two widely introduced *Cirsium* species as contemporary reference distributions and asked whether changes in climatic occupancy are accompanied by changes in multivariate, photograph-derived capitulum phenotypes.

**Methods / test logic (working draft):**

> For each species, we quantified climatic occupancy overlap between its native and introduced ranges using a shared five-variable climatic ordination and accounted for differences in geographically available climate. We then compared the distributions of image-referenced capitulum angle, photographic corolla chroma and projected outline shape using taxon- and photograph-linked records. Climatic and phenotypic overlap were examined separately for whole ranges and, using exactly matched eligible geographic cells, across pairs of native and introduced regions. Region-label permutations accounted for the reuse of regions across pairwise comparisons. We additionally tested whether capitulum phenotypes in introduced occurrences outside the native range's main 95% climatic occupancy region differed from those in climatically analogous introduced occurrences, using within-region matched-cell contrasts. Comparisons restricted to angle and chroma assess sensitivity to the rotation dependence of the projected outline metric. We interpret these tests as observational associations, not as evidence of selection or genetic evolution after introduction.

**Discussion / evidence-limited formulation (working draft):**

> The climatic environments occupied in the introduced ranges differed from those occupied in the native ranges, whereas photographed capitulum phenotypes retained substantial overlap. However, climatic and phenotypic overlap across paired regions did not show a detectable monotonic association, and lower phenotypic overlap in climatically novel introduced occurrences was not supported after spatial matching and multiplicity correction. These findings are compatible with colonization of different realized climatic environments without large, consistently detectable changes in the measured capitulum phenotype space, but they do not establish phenotypic invariance or ecological fitting as a causal mechanism. Inference remains limited by the number of focal species, uncertain introduction sources, uneven geographic replication and photographic measurement error.

上記のDiscussion草案は今後の背景補正済み地域ペア解析や独立形質QCで結果が変われば修正する。特に現在の**地域ペアrhoは観測地点の気候D**であって背景補正済みニッチDではないため、Methods/Resultsに2つの推定対象を明記して混同しない。

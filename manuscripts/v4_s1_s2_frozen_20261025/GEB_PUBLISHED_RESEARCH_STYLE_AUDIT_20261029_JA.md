# GEB 掲載論文の比較読解と本論文の位置づけ

更新日：2026-10-29。対象：PR159 の `5bc37b22adf4d8411fd7814df1aaa58de3d94fcb` を基点とする原稿改訂。

## 1. 本論文の位置づけ

**公開写真から連続的な頭花形質、同じ写真記録における形質の組合せ、その観測位置を保持する計測・比較の仕組みを構築し、それによって繁殖器官の外観が気候に沿って地理的に配置されることを示すマクロ生態学研究。**

GEB の Research Article として、方法と生態学の貢献を一つの論理で結ぶ。方法の価値は、8 形質を計算した数そのものだけでなく、連続値・形質の共起・観測場所を同時に保持した点にある。生態学的な中心は、その情報を用いた H3c が、各分類群の選択されたセルにおける多形質分布を保持した帰無比較に対して、気候の似たセル間でより大きな頭花形態の類似を検出した点にある。

この構造では、3.1 の降水と頭花角度、3.2 の日射・標高と花色、3.3 の局所的な多形質類似が、それぞれ違う比較単位から気候と繁殖器官の対応を記述する。3.2 が示すのは各分布域内の関係であり、在来／導入間の傾きの差を検定した結果ではない。3.3 の主要な問いは「実際に似た気候を占める地点同士では、頭花形態も似るか」で統一し、幅・重心・共分散・KDE の補足結果は SI に保持する。

以下の執筆方針は掲載論文の比較から導いた**編集上の判断**であり、他論文の結論を本研究の実証結果に移したものではない。

## 2. 比較読解した文章とアクセス範囲

比較対象は **14 論文**。うち **12 本は主本文の Introduction、Methods、Results、Discussion を確認**した（GEB 11 本、Ecology and Evolution 1 本）。残る 2 本は公開範囲の限定された資料として扱う。GEB の掲載論文は全体で 12 本であり、「14 本すべてが GEB」「14 本すべてを全文読解」とは数えない。全論文の SI を読解・再解析したという意味でもない。

### A. 主本文を確認した 12 本

| 論文・一次資料 | 問いと観測の組み立てから学んだ点 | 今回の反映 |
| --- | --- | --- |
| **Kuppler et al. (2020)**. *Global gradients in intraspecific variation in vegetative and floral traits are partially associated with climate and species richness*. GEB 29:992–1007. [DOI](https://doi.org/10.1111/geb.13077) | 種平均で残る問いから、種内変異という応答量を定義し、形質ごとの環境関係へ進む。考察は主要なパターンを述べてから種・形質ごとの差を説明する。対象は種内変異の「大きさ」であり、形質値の環境勾配と同一ではない。 | 序論で、変異の幅に加えて「その表現型がどこにあるか」という観測情報を明示。 |
| **Govaert et al. (2024)**. *Trait–micro-environment relationships of forest herb communities across Europe*. GEB 33:286–302. [DOI](https://doi.org/10.1111/geb.13789)／[著者最終稿](https://backoffice.biblio.ugent.be/download/01HGCZ5NSJF0Q5K8X397ZCC56M/01HKW771C397GKW8N6TSNMV3GD) | 局所測定を保持することの意義を、環境と形質の対応として具体化する。観測の細かさを、研究の問いと同じ段落で説明する。 | 地理的に反復された形質測定の必要性を序論第2段落へ。現地の微環境測定と本研究の気候平年値は区別する。 |
| **Hanz et al. (2022)**. *Climatic and biogeographical drivers of functional diversity in the flora of the Canary Islands*. GEB 31:1313–1331. [DOI](https://doi.org/10.1111/geb.13507) | 同じ環境要因でも、植物群の由来と地域史を含む比較によって解釈が定まる。序論で広い問いを置き、対象系がその問いをどう具体化するか説明する。 | 地域史と環境を併記して開始し、Cirsium の地域的放散を研究系の理由へ接続。 |
| **Wang et al. (2023)**. *Disentangling thermal from alternative drivers of reflectance in jewel beetles: A macroecological study*. GEB 32:408–420. [DOI](https://doi.org/10.1111/geb.13632)／[掲載PDF](https://devistuartfox.com/wp-content/uploads/2023/07/2022-wang-disintangling.pdf) | 追加の波長域を測る理由を生態学的な識別問題から説明し、方法の冒頭でも測定の利点を機器仕様より先に示す。考察は得られた結果と測定が可能にした解釈を結ぶ。 | 写真計測が保持する情報と H3c の比較可能性を直結。較正された反射率と公開 RGB 写真の外観の違いも整理。 |
| **Kamoske et al. (2022)**. *Towards mapping biodiversity from above: Can fusing lidar and hyperspectral remote sensing predict taxonomic, functional, and phylogenetic tree diversity in temperate forests?* GEB 31:1440–1460. [主本文](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13516) | 各応答次元の生態学的意味を先に説明し、観測情報を組み合わせる必要性へ進む。考察・結論で方法と生態学への貢献を両方回収する。 | 3 領域・8 構成概念の意義を、単なる多次元化ではなく局所的な組合せを比較するための情報として記述。 |
| **Hakkenberg & Goetz (2021)**. *Climate mediates the relationship between plant biodiversity and forest structure across the United States*. GEB 30:2245–2258. [主本文](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13380) | 序論末の問いを分析と考察に対応させ、条件によって変わる関係自体を成果として扱う。観測の価値を述べた後で機構解明に必要な実験へ進む。 | 3 つの問いと 3 つの結果節を対応させる。ただし本研究の S2 には交互作用検定がないため、同論文の「条件依存性を検定した」という主張は移さない。 |
| **Kellenberger, Winner & Jetz (2026)**. *The Performance and Potential of Deep Learning for Predicting Species Distributions*. GEB 35:e70184. [主本文](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.70184) | 技術の価値をモデル名や流行ではなく、公平な比較で示した能力に結び付ける。数値性能と生態学的理解を区別する。 | YOLO の使用自体を新規アルゴリズムとせず、実際に保持・比較できた情報を方法論の貢献とする。 |
| **Ringelberg et al. (2020)**. *Biomes as evolutionary arenas: Convergence and conservatism in the trans-continental succulent biome*. GEB 29:1100–1113. [主本文](https://onlinelibrary.wiley.com/doi/full/10.1111%2Fgeb.13089) | 離れた地域の類似について競合する解釈を立て、各解釈に対応する予測と比較を区別する。方法の役割を説明した後に地域史の議論へ戻る。 | 「離れた地域での類似」という広い位置づけを明示。本研究はセル単位の表現型対応を示し、系統を用いた進化的収束の検証は区別する。 |
| **Riera et al. (2026)**. *From Purpose to Persistence: How Intentional Introductions by Humans Influence Ecological Niches of Non-Native Plants*. GEB 35:e70245. [主本文](https://onlinelibrary.wiley.com/doi/10.1111/geb.70245) | 考察冒頭で問いへの答えと効果の大きさを述べ、次にどの次元で成立したかを特定する。由来を示す群分けと因果的な過程を分ける。 | 支持された関係を先に述べ、その対象・比較単位をすぐ添える。分布域区分から導入後の変化を自動的に推論しない。 |
| **Anadon-Rosell et al. (2025)**. *Intraindividual Variability as a Large Source of Trait Variation in Clonal Tundra Dwarf Shrubs Along Elevation and Latitude Gradients*. GEB 34:e70134. [主本文](https://onlinelibrary.wiley.com/doi/10.1111/geb.70134) | 種間・種内・個体内という観測階層を順に導入し、変異を保持したこと自体の生態学的意味を示す。各形質で同じ環境関係が成立する必要はない。 | 観測単位を明確化し、各形質の同一勾配を要求せずに、多形質の局所対応を論じる。 |
| **Chen et al. (2025)**. *‘Slow-Fast’ Plant Trait Spectra Are Associated With Ecological Niches Across Global Climatic Gradients*. GEB 34:e70115. [主本文](https://onlinelibrary.wiley.com/doi/10.1111/geb.70115) | 序論末の複数の問いを結果・考察節に対応させる構造が参考になる。予測と異なる結果を明示してから説明を検討する。 | 問いと節の対応という構造のみ参照。葉形質とニッチ幅の方向は要旨と Results に不一致があり、本研究の事実的根拠に使わない。進化的因果表現も採用しない。 |
| **Laitly et al. (2021)**. *Is color data from citizen science photographs reliable for biodiversity research?* Ecology and Evolution 11:4071–4083. [主本文](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.7307) | 測定法の有用性を検証対象と観測階層に対応づけ、結果の利点から適用範囲へ進む。 | 写真による色の定量化の先行研究として使用。鳥類での測定対応・必要写真数を Cirsium の検証精度へ転用しない。 |

### B. 公開範囲を限定して確認した 2 本

| 論文 | 実際に確認できた範囲 | 使用方法 |
| --- | --- | --- |
| **Schiller et al. (2021)**. *Deep learning and citizen science enable automated plant trait predictions from photographs*. Scientific Reports 11:16395. [一次資料](https://doi.org/10.1038/s41598-021-95616-0) | 要旨、一次資料の検索索引に含まれる Methods／Discussion 抜粋、出版情報。主本文ページのアクセス制限と約40MBのPDF制限により、今回の全文読解数に含めない。 | 種名を通じて結び付けた訓練情報による形質予測の先行研究として位置づける。繁殖器官の検出や頭花の直接幾何計測を実証した文献とは書かない。 |
| **Kim, Bae & Kang (2026)**. *Plumage Polymorphism Varies Along Climatic and Urban Gradients Across Global Cities*. GEB 35:e70313. [出版社ページ](https://onlinelibrary.wiley.com/doi/10.1111/geb.70313) | 公開要旨・出版情報・SI の説明。主本文の Introduction／Discussion は未確認。2026年9月29日公開で、今回の作業日より前の掲載論文。 | 公開写真による色彩モルフと環境のマクロ生態学的先行研究として認める。本文の文章構成を読んだとは扱わず、要旨で確認できる内容だけを引用。 |

## 3. 文章構成として採用したこと

### 序論：生態学的な問い → 欠けている観測情報 → その情報を保持する方法 → 対象系 → 比較

第一段落は「似た環境で似た形が現れるか」という地域をまたぐ問題から始める。続いて、種平均や変異幅では保持されない地理的な配置を特定する。第三・第四段落では写真研究の先行例を認めたうえで、連続値、同じ記録内の組合せ、位置情報を同時に持つ必要性を説明する。最後に Cirsium の共有構造と地域的広がりを示し、3 つの問いを結果節に対応させる。

この順序は、Wang の「測定によって識別できる問い」、Kamoske の「観測次元ごとの意味」、Ringelberg の「離れた地域の比較を生物学的問いから導く」という構造を、本研究の実証範囲に合わせて応用した。原文の特徴的な表現は転写していない。

### 方法：装置・アルゴリズム名より、得る値と比較単位を先に明示

2.2 は検出 → 連続計測 → 写真記録への集約 → 地理的比較という順に改訂した。YOLO は候補領域を検出し、色・幾何量はその後の計算で得る。色相は HSV、lightness/chroma は CIELAB である。領域ごとに適格な頭花が異なり得るため、共通単位は「同じ写真記録」とする。8 構成概念を9数値座標で表す理由と、各解析が異なる適格コホートを使う理由を明記した。

詳細な測定定義は新設の [Table S4](SI_PHOTOGRAPHIC_TRAIT_MEASUREMENT_DEFINITIONS_20261029.md) に整理した。H3c の気候3軸、セル集約、標準化、rarefaction、whole-vector permutation、original Holm-2 family は変更していない。

### 結果：計測によって組み立てられた観測情報と、生態学的な発見を両方述べる

3.1 の冒頭を、156,683 写真・174 分類群で3領域の連続値を共通記録に保持できたという実証的な到達点から始める。続く BIO18 の係数・信頼区間・検定と spatial sensitivity は変更しない。3.2 は全文を保持する。3.3 は地理的な比較単位の表現を正確化し、肯定的な主結果・数値と、その統計的安定性を限定する実際の感度は保持する。

### 考察：得た答え → 何が新しく分かったか → 比較単位に即した解釈 → 次の検証

冒頭では、地理的配置の気候的構造と、それを調べられる計測の仕組みを一つの貢献として述べる。4.3 は、分類群内で完全な形質ベクトルを並べ替えても保たれるものと、並べ替えによって失われる場所との対応を説明する。4.4 は **“What public photographs add to floral macroecology”** とし、共通コホート上の28形質ペアと気候類似地点比較という、実際に可能になった解析を示す。較正・視点・照明の問題は、その情報を野外形質に接続する具体的な課題として置く。

### 結論：方法の能力と生態学的発見を回収

最後は、地理的に反復された多形質測定によって、個別形質の勾配とセル単位の多形質類似を接続できた、と結ぶ。今後の課題だけを最終的なメッセージにせず、今回残った実証的な貢献を明示する。

## 4. 本論文で主張する方法論の成果と根拠

| 方法論の成果 | 確認済みの根拠 | 生態学的な利用 |
| --- | --- | --- |
| カテゴリー内の中間状態を連続値で保持 | 角度、lightness、chroma、円周 hue、4 種の輪郭量という8構成概念。 | 方向・色・輪郭の個別環境勾配を比較。 |
| 同じ写真記録で複数領域を対応づける | 156,683 写真・174 分類群で全28構成概念ペアを同じ記録上に保持。 | 地理的に分散した別々の種平均から合成せず、観測された組合せを比較。 |
| 適格性と観測階層を明示 | 頭花→写真記録→分類群×セル。各応答の元のフラグと環境完全性を保持。 | 観測数と分類群・セル数を区別した推論。 |
| 分類群の標本分布を保った局所比較 | H3c は選択された8セルの3形質ベクトルを分類群内で並べ替え、多変量分布と共変動を保持。 | 気候との地理的対応が分類群ごとの標本分布に加える情報を評価。 |
| 数値結果と実装を追跡できる | 元の係数・source hashes・208 回帰検定・H3c aggregate と source-linked scripts。 | 文章を更新しても測定値・比較family・分析単位の対応を維持。 |

これらは実装・データ構造・実際の比較から示せる成果である。新しい検出ネットワークの発明や、生物学的な測定誤差の独立較正は別の主張になる。旧版の反転画像・摂動試験や小規模な植物学的カテゴリー比較は、異なるコホートの技術的検討であり、現在の V4 の精度検証値へ置き換えない。

## 5. 生態学的成果を最大限に述べるための比較単位

| 主結果 | 明確に述べる生態学的知見 | 解釈の単位 |
| --- | --- | --- |
| 3.1 BIO18–angle | より湿潤な暖候期の気候を占める在来分類群ほど、写真上でより下向きの頭花を示した。 | 分類群間。重力基準の頭花角・降雨時の機能は別測定。 |
| 3.2 radiation/elevation–colour | 2 種の分布域内で、日射・標高に沿う花色の連続的な変化を捉えた。 | 各 species×range stratum 内。有意性の違いだけでは在来／導入の傾き差を示さない。 |
| 3.3 climate analogue–joint resemblance | 実際に気候が似るセル同士では、分類群ごとの標本形態分布を保った帰無期待より頭花形態が似た。 | 異なる地理セクターに割り当てられた分類群間のセル比較。地理的配置が加える情報。 |

H3c の局所対応は、3形質それぞれが同じ機構・同じ方向で変化するという主張を必要としない。5°で統計的支持がより安定していても、コホートと標準化が異なるため「生物学的効果が約2倍」「5°が最適な生態学的スケール」とは解釈しない。分類群の所属セクターは優占する占有セル数で決めるため、個々の気候一致セルが必ず異なるセクターにあるわけではない。これらの範囲を保つことで、主結果を適応・収束の証明へ広げずに、生態学的な意味を明確にできる。

## 6. 改訂ファイルと確認方針

- [主原稿](GEB_RESEARCH_ARTICLE.md)と blinded 版：要旨、序論、Methods 2.2、Results 3.1 の計測成果段落、考察の導入・4.3・4.4、結論、関連文献、Figure 1 の構成を改訂。H3c の地理的比較単位と帰無比較が保持する分布の表現を精密化。
- [Table S4](SI_PHOTOGRAPHIC_TRAIT_MEASUREMENT_DEFINITIONS_20261029.md)：測定の操作的定義、同じ写真記録という単位、元コードの対応を新設。
- [Cover-letter contribution draft](GEB_COVER_LETTER_CONTRIBUTION_PITCH.md)：方法と生態学の貢献を接続する説明へ改訂。送信・投稿はしていない。
- 3.2 全文、3.1 の気候関連係数段落、3.3 の統計結果段落を基点と比較し、元の数値・解析コード・source JSON・Holm family を保持する。
- 既存の manuscript contract tests を実行し、見出し変更に対応する期待文字列のみ更新する。測定精度・GEB 採択可能性を文章整合性テストの合格から推論しない。

出版用の図、実データからの採集地図、写真の権利確認、匿名の永続的データアクセスは既存の投稿準備項目として残る。本改訂は原稿の学術的位置づけと方法・結果の説明を完成させるものであり、それらの完成を装わない。

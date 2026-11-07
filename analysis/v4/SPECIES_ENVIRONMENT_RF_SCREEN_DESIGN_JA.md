# Cirsium 種平均ランダムフォレスト：探索用環境変数ランキング（2026-11-07）

**実行目的**：新しい因果・確認的検定を加えるのではなく、**どの気候・地形変数、あるいは相関の強い変数群が、種間の頭花向き・色・輪郭の違いの予測に情報を持つか**を探索する。従来の26種RF（種内セルを8個以上必要とした）と推定対象が異なる。

**今回の状態**：コードと合成テストは完成。313,350観察に基づく種平均RFの**実測は未実行**。公開GitHub成果物は既存モデルの匿名集計のみで、分類群×形質×環境の結合行列を含まない。根拠のない変数重要度や順位を作らない。

## 解析単位と候補

- **1分類群=1行**。元S1のQC済み写真から各形質平均を計算。頭花の画像角度、L*、C*、円周色相sin/cos、縦横比、円形度、Solidity、幅のCVの**全8構成形質**。
- 各分類群で最少5適格写真。**8地理セル・40写真の閾値は用いない。** 244は元S1の**向き適格分類群数**で、8形質共通Nではない。結果のNは形質別に記録する。元撮影者非重複・種平均信頼性QCは未計算であり、RFで代用しない。
- 既存、復元・検証済み23候補：CHELSA BIO01–BIO19（19）、RSDS（1）、VPD（1）、地表風（1）、WorldClim標高（1）。**23変数全部に登録行**を設ける。未取得列・有効分類群70%未満・定数は「未使用と理由」を残す。土壌34/44提案の未取得フル種平均層をあるように扱わない。
- 環境の種平均は、**同じ形質のQC適格写真集合**に対して各環境変数ごとに計算する。環境の部分欠損で形質平均や他環境の支持観察を削らない。種平均の観察重みをもつ写真観測が地理的に偏り得ることは留保。
- 色相は元写真のsin/cosを共同単位化→種平均→方向単位化（長さ0.2未満の不定な種は除外）。角度0–360を直線として回帰しない。

## 多変量RFと出力

- 30分類群以上の解析対象に、**固定設定** RandomForestRegressor（120 trees、min_samples_leaf=4、max_features=0.75、1 thread）。形質ごとにハイパーパラメータ探索を行わない。
- 分類群単位5-fold×3反復を固定seedで繰り返す。全分類群がテスト側に回るたび、**訓練分類群の平均形質**を参照とするSkill=1−RF損失/参照損失。負のSkillは「RFの外部予測が悪い」であって、環境との関連がゼロの証明ではない。
- 欠測環境の補完は**各学習foldだけの中央値**。テスト側の情報で補完しない。学習結果の環境依存性を見るため、未知分類群のテスト側の変数をシャッフルし、損失増分を参照損失で規格化した**held-out permutation importance**を計算。
- **三層の重要度を併記**：①各変数単独入替（相関変数間で寄与を奪い合う）、②生態学的な固定群（温度・降水・日射・VPD・風・標高）を同じ行順で一括置換、③訓練側の|Spearman r|≥0.80／complete linkageから得た相関群を同時置換。後者は各foldで変数の束が変わり得るので、同じ組合せに対して観測したfold数も表示する。
- 変数順位の頻度（top 5に何割入るか）、正のPermutation importanceの割合、平均と分位点を表示。**高い重要度はそのRFへの依存であり、独立な効果・生態機構・P値ではない**。全体Skillが負なら候補としての解釈を弱める。
- 元位置の経度が許可済み入力に存在し、経度3帯それぞれ10分類群以上なら、**優占する観察経度帯ごとに3分割の保留評価**を補助表示。経度帯は系統や大陸の独立反復ではない。座標や分類群名は出力しない。
- 表は匿名集計のみ。元S1・撮影者・場所・種名・個別予測・RFモデル本体を公開しない。従来のV4論文・S1/S2 Holm・H3cや既存RF、再現性QCを改変しない。

## 既存原本からの実行

既存の `species_environment_existing_inputs.load_existing_inputs` をそのまま使用し、元S1のデコード表＋外部7＋追加9と正規receiptをチェックする。データ取得、復号、Actions操作、秘密鍵取得・権限変更は**このコードには含まれない**。

許可済み環境で以下を実行（パスは説明用、実在未確認）。

```bash
python -m analysis.v4.species_environment_rf_screen \
 --native /authorized/decoded/observations.csv.gz \
 --decoded-receipt /authorized/decoded/measurement_receipt.json \
 --seven /authorized/environment/original_seven.csv.gz \
 --nine /authorized/environment/chelsa_extra_nine.csv.gz \
 --nine-manifest /authorized/environment/chelsa_nine_manifest.json \
 --source-receipt /authorized/public/bio21_source_receipt.json \
 --source-receipt-sha256 26b80b279fa2b372099ef945972d3e4fdf58127fe9f8fd92530ac08a931d9545 \
 --output /authorized/public/species_environment_rf_screen_aggregate.json
```

**公開するのは** `species_environment_rf_screen_aggregate.json` のみ。既存receiptの表SHAに合わない原本を無理に通さない。出所が異なる版は、その版の正規レシートを上流で発行して区別する。

## 結果の読み方

報告するのは形質ごとに：N分類群、使えた環境変数、未知分類群Skill、変数別・固定環境群別・相関群別の予測損失増分、15foldの順位頻度、必要なら3地域除外Skill。実測でない仮想データの順位を生態学的発見としない。

**既知の比較対象**：V2/V3では BIO12×向き、RSDS×彩度に回帰で残る関連がある。ただしV4のBIO18×向きは地理構成で大きく弱まるので、今回のRF結果を既存の降水–向きの因果的確証へすり替えない。新しいRFはあくまで探索用であり、本文の主解析ではなく、後の独立したモデル・新データで候補を評価する。

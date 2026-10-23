# T3：実資料から増補したSenecioneae頭花向き候補（2026-10-23）

## 確認した実資料と現在の生物学的結果

従来のT3には、CardueaeのCirsium核系統研究に対して他連SenecioneaeのCremanthodium campanulatum人工頭花操作の単一機能例があった。今回は、**公開された種別の頭花（茎ではなくcapitulum）形態原記載**を独立に照合し、[9種のエビデンス表](20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv)に登録した。

**形態記載上の6下垂候補**：Cremanthodium campanulatum、C. arnicoides、C. spathulifolium、Ligularia lingiana、L. fangiana、Senecio eliseae。

**形態記載上の1直立候補**：Xenophyllum esquilachense（capitula erect）。

**判定しない2種**：Xenophyllum staffordiaeの「rather nodding」は程度を限定する語なので二値Dとして未採用。一方、Senecio humbertiiの原文「at first erect, later nodding」は発生過程で状態が変わるため種の固定U/Dにできない。

この追加により、**Senecioneae内にも直立と下垂の形態記載がある**こと、および種内の発生的姿勢変更を無視してU/D固定化すると進化起源数の推定が誤り得ることを具体的な公刊資料で確認した。ただし、**9種の記載 ≠ 9回の独立進化**であり、属内の複数下垂種をそれぞれ独立な起源として数えない。

## 原文へのリンクと根拠の区別

- **直接実験**：[Chen et al. 2013：Cremanthodium campanulatum](https://doi.org/10.1080/17550874.2012.702793)。頭花を人為的に上向きにした処理で痩果形成減少、雨・UV-Bが花粉に有害な可能性。これは当該種の機能的利益であり、科全体での反復起源の因果実証ではない。
- **Flora of China**：[Cremanthodium arnicoides](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200023734)、[C. spathulifolium](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200023794)、[Ligularia lingiana](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200024246)、[L. fangiana](https://efloras.org/florataxon.aspx?flora_id=2&taxon_id=200024218)。
- **種原記載**：[Senecio eliseae](https://doi.org/10.11646/phytotaxa.230.1.6)。
- **分類改訂**：[Xenophyllumの直立例と限定的下垂例](https://pmc.ncbi.nlm.nih.gov/articles/PMC7483337/)。
- **開花段階の違い**：[Senecio humbertii Flora of China](https://efloras.org/florataxon.aspx?flora_id=2&taxon_id=242348317)。

機械的な監査は [t3_source_evidence_gate.py](../t3_source_evidence_gate.py) で実施。現状では9件いずれも植物誌原文の種名と分子標本の同一性、開花時点、独立2名採点、根付き核種系統のtip照合が未完了のため、**T3祖先状態・Mk・降水選択検定は実行していない**。既存T2の核樹結果を含めても現時点で連間収斂を主張できない。

## 次の核系統入力の候補（実データ取込はまだ）

- [Mandel et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6628808/)：Asteraceaeの256分類群・935核座位の基幹系統。
- [2024年Asteraceae核系統研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11211554/)：706種・41連の核系統。

**いずれも今回の9種とCirsiumの現行対象名が完成した機械可読Newickに同時に含まれることをまだ検証していない**。属の配置だけを使って種内の起源回数や降水応答を推論しない。次に、種レベルtip被覆を実物treeで監査し、対象が少ない場合は東アジア/南米の同属独立核樹への拡張を行う。

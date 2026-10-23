# T3：実資料から増補したSenecioneae頭花向き候補（2026-10-23）

## 確認した実資料と現在の生物学的結果

従来のT3には、CardueaeのCirsium核系統研究に対して他連SenecioneaeのCremanthodium campanulatum人工頭花操作の単一機能例があった。今回は、**公開された種別の頭花（茎ではなくcapitulum）形態原記載**を独立に照合し、[10種のエビデンス表](20261023_T3_SENECIONEAE_POSTURE_SOURCE_CANDIDATES.csv)に登録した。

**形態記載上の6下垂候補**：Cremanthodium campanulatum、C. arnicoides、C. spathulifolium、Ligularia lingiana、L. fangiana、Senecio eliseae。

**形態記載上の1直立候補**：Xenophyllum esquilachense（capitula erect）。

**種の固定U/Dへ変換しない4件の論点**：
- *Xenophyllum staffordiae* は「rather nodding」に加えて、**植物体が成長すると下垂する傾向**が記載される。最初の状態・同一頭花の時系列追跡は不明なので、依然として `?`、ただし年齢依存フラグを追加した。
- *Senecio humbertii* は「at first erect, later nodding」、*S. nigrocinctus* は「at first erect, then nodding」と記載。**少なくとも2種で生育段階に沿う方向変化の記載**がある。ともに `V` として固定U/Dから除外した。
- *Senecio eliseae* の2015年原記載は下垂だが、2018年のSalomónらのsynonymyと、2022年のCalvo & Grandaによる異議がある。元の花姿勢を削除せず、分類概念の対応未確定を明示した。

この追加により、**Senecioneae内にも直立と下垂の形態記載がある**こと、および種内の発生的姿勢変更を無視してU/D固定化すると進化起源数の推定が誤り得ることを具体的な公刊資料で確認した。ただし、**10種の記載 ≠ 10回の独立進化**であり、属内の複数下垂種をそれぞれ独立な起源として数えない。

## 原文へのリンクと根拠の区別

- **直接実験**：[Chen et al. 2013：Cremanthodium campanulatum](https://doi.org/10.1080/17550874.2012.702793)。頭花を人為的に上向きにした処理で痩果形成減少、雨・UV-Bが花粉に有害な可能性。これは当該種の機能的利益であり、科全体での反復起源の因果実証ではない。
- **Flora of China**：[Cremanthodium arnicoides](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200023734)、[C. spathulifolium](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200023794)、[Ligularia lingiana](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=200024246)、[L. fangiana](https://efloras.org/florataxon.aspx?flora_id=2&taxon_id=200024218)。
- **種原記載**：[Senecio eliseae](https://doi.org/10.11646/phytotaxa.230.1.6)。
- **分類改訂**：[Xenophyllumの直立例と限定的下垂例](https://pmc.ncbi.nlm.nih.gov/articles/PMC7483337/)。
- **開花段階の違い**：[Senecio humbertii Flora of China](https://efloras.org/florataxon.aspx?flora_id=2&taxon_id=242348317)、[S. nigrocinctus Flora of China](https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=242348425)、[X. staffordiae 分類改訂](https://pmc.ncbi.nlm.nih.gov/articles/PMC7483337/)。
- **原学名・種境界の不一致**：[Salomón et al. 2018](https://phytotaxa.mapress.com/pt/article/view/phytotaxa.340.1.1)は *Senecio eliseae* を *S. rhizomatus* に統合。一方 [Calvo & Granda 2022](https://doi.org/10.15553/c2022v771a1)はその異名処理を支持しない。T3では今のところ分類上の二つの見解を併記し、どちらかを独断でaccepted概念へ強制しない。

機械的な監査は [t3_source_evidence_gate.py](../t3_source_evidence_gate.py) で実施。現状では10件いずれも植物誌原文の種名と分子標本の同一性、開花時点、独立2名採点、根付き核種系統のtip照合が未完了のため、**T3祖先状態・Mk・降水選択検定は実行していない**。既存T2の核樹結果を含めても現時点で連間収斂を主張できない。

## 次の核系統入力の候補（実データ取込はまだ）

- [Mandel et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6628808/)：Asteraceaeの256分類群・935核座位の基幹系統。
- [2024年Asteraceae核系統研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC11211554/)：706種・41連の核系統。

**いずれも今回の10種とCirsiumの現行対象名が完成した機械可読Newickに同時に含まれることをまだ検証していない**。属の配置だけを使って種内の起源回数や降水応答を推論しない。次に、種レベルtip被覆を実物treeで監査し、対象が少ない場合は東アジア/南米の同属独立核樹への拡張を行う。

## 修正済み集計と科学的な停止規則（2026-10-23）

真正の独立著者形態資料10件を保持したうえで、固定下垂記載6件、固定直立記載1件、弱い年齢依存傾向で二値化不能1件、発生段階でU→Dと記載された2件に分けた。少なくとも1件の著者間分類概念異説がある。複数種の同じ属・異名があっても独立系統イベントへ数えない。

[機械可読の監査コード](../t3_source_evidence_gate.py) は、(i) 年齢依存・開花段階変化のレコードを固定U/Dに転記、(ii) 著者間の異名論争を無視して同一tipのaccepted conceptへ昇格、(iii) 盲検二重採点・種concept・根付き独立核樹の照合前にT3の厳密比較へ昇格、をエラーで停止する。この追加は**実在する形態知識の適切な入力整理**であり、同じ頭花の経時的姿勢を実測した結果、古気候を復元した結果、独立反復を同定した結果ではない。

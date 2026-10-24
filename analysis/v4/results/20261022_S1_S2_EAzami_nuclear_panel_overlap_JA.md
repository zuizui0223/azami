# S1/S2とEAzami公開核系統・気候ニッチ標本の実際の重複照合（2026-10-22）

**結論：EAzamiの東アジア20-tip樹はS1の連続画像角度244分類キーと厳密に5分類群だけ重なり、S2の2種とは0種。だがMoreyra et al. 2025の全球核データの公開標本台帳にはS2の両種がある。** 「公開分子データがない」という意味ではない。現在のS1/T2用52配置樹とEAzami実核系統樹、原著のMoreyra全球最終樹は別資源。

## 根拠を固定

- **S1実際の元画像と環境：** [source-verfied GitHub Actions #39584744194](https://github.com/zuizui0223/azami/actions/runs/39584744194) の公開aggregate artifact **12339413151** (`S1_T2_joint_244_194.json`)。`T2_independent_state_inventory.taxon_source_audit` の `s1_validated_taxon_match=true`、`s1_observations`、`s1_phylogeny_194_key_match` を使用。原分類名をWCVP accepted keyへ厳密対応。写真角度や観察数は実データ集計であり、文献の植物学的U/D状態とは別。
- **東アジア20-tip核系統と文献姿勢：** [EAzamiの20-tip source crosswalk](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/orientation_comp1061_20tip_source_crosswalk_v1.csv)、[20-tip核系統Newick](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/full20_comp1061_primary_tree_v1.nwk)。この樹はEAzami側の独自Comp1061再構築で、原著Moreyra 350核座位の最終Newickではない。
- **気候ニッチに使った9/12分類群：** [EAzami source-frozen taxon coverage result](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/chapter2_orientation_occurrence_coverage_audit_result_v1.json)。9 = 1種あたり10以上、12 = 5以上のsource-thinned元気候環境記録というEAzami側の条件。これとS1元写真枚数の閾値を混同しない。
- **S2主比較：** [元S2全導入域比較](20261021_unified_abcde/S2_ALL_INTRODUCED_DIRECT_COMPARISON_20261021_JA.md) の *C. arvense* と *C. vulgare*。
- **全球Moreyra原公開標本：** EAzamiの[258件public SRA元標本のsource-linked reconciliation](https://github.com/zuizui0223/EAzami/tree/main/data/evidence/moreyra2025_cirsium_reconciliation_v1) に *C. arvense*（変種名を除く1件、変種名2件）・ *C. vulgare*（原tree_code `Cirsium vulgare` 1件と番号付き`Cirsium vulgare 1–5` の5件）の登録がある。両種の公開SRA runが元標本と紐付く。ただし原著最終木のtip/枝長/支持値がこのrepositoryに保存されていることは意味しない。[未回収の最終全球樹・最終350座位監査](https://github.com/zuizui0223/EAzami/blob/main/docs/MOREYRA_2025_FINAL_TREE_RECOVERY_LOG_2026-08-10.md)。

## 科学的な重複行列

| EAzami 20-tipにある分類群で S1 と exact accepted-key matching | S1元画像角度の観察数 | S1 194-tip配置樹にもある？ | EAzami日本7種パネル | 東アジアn>=10・9種 | 東アジアn>=5・12種 |
|---|---:|---|---|---|---|
| *C. kamtschaticum* (accepted key 206742945) | **385** | はい | はい | はい | はい |
| *C. kawakamii* (206744933) | **23** | はい | いいえ | はい | はい |
| *C. tatakaense* (206738470) | **95** | はい | いいえ | いいえ | はい |
| *C. morii* (206738938) | 1 | はい | いいえ | いいえ | はい |
| *C. suffultum* (206740495) | 1 | **いいえ** | はい | はい | はい |
| **厳密な重複合計** | **5分類群／505元観察** | **4分類群** | **2** | **3** | **5** |

EAzami核樹全20分類群のうち残り15分類群は、今のS1への**独立受容分類群としての元画像角度照合には成功していない**。とくに *C. brevicaule* と *C. irumtiense* の元ソース名はS1内で**同じ受容キーに統合される衝突**があり、両種を別のS1 tipと認定しない。*C. yezoense* も別の原学名との衝突で厳密照合に不適格。*C. alpicola* や *C. gyojanum* などには元S1核系統に独立な写真角度結果がない。このことは当該種に写真や野外分布データが一切存在しないことの証明ではない。

| 分析対象 | 東アジア20-tip | S1の244写真適格 | S1の194-tip配置 | S2の2種 |
|---|---:|---:|---:|---:|
| 東アジア核系統樹との厳密な一致 | 20 | **5** | **4** | **0** |
| 東アジア日本7種ニッチ群からのS1一致 | 7 | **2** | **1** | **0** |
| 東アジアn>=10・9種ニッチ群からのS1一致 | 9 | **3** | **2** | **0** |
| 東アジアn>=5・12種ニッチ群からのS1一致 | 12 | **5** | **4** | **0** |

## BM/OUと「ニッチ放散」の次の判定

- EAzamiの **9/12分類群に対するBM/OU比較**は、その日本・台湾群での**現生CHELSAニッチ値の系統モデル**を検査できても、種集合が主に異なるため、**azami S1全球244分類群やS2導入2種に対する直接検証にはならない**。
- 実際のS1画像角度とEAzami核樹を両方持つ5分類群のうち、**2分類群の画像観察は各1件**。写真10件以上で使えるのは *kamtschaticum/kawakamii/tatakaense* の**3分類群**だけ。この重複だけのOU vs BMモデル比較から、頑健な適応最適値や放散様式を推定することはできない。
- S2との共通の核系統情報の経路はEAzami東アジア20-tip樹ではなく、**Moreyra全球の公開配列・標本と正式な枝長付き樹**。原著最終樹の公開回収監査は2026-08-10時点で未完了。公開原配列から新たな再構成を行う場合も、**原著の最終350核座位系統を再現した**とは自称しない。
- ブラウン運動（BM）によるPGLSや祖先値補間を実行したことは、**BMがOUより適合すること**や、過去のニッチが放散したこと・適応的レジーム転換があったことを証明しない。OU比較には同じtip集合・形質/ニッチ定義・枝長/根の根拠・適切な小標本ペナルティとモデル診断が必要。

**次の実装の判断：** T2核系統樹のスコープを元S1/S2の対象種との重複表から先に定義し、S1/2と接続しない局所少数種のみのOU結果を「全球S1–S2を説明した」と主張しない。

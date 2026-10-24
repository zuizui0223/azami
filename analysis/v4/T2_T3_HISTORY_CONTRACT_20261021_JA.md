# T2/T3 実装契約｜降雨と頭花向きの進化史（2026-10-21）

> **2026-10-24現行科学方針（旧U/D主線の停止）：** [S1の連続画像角度を用いるT2・T3の主解析](S1_CONTINUOUS_T2_T3_RESEARCH_MAINLINE_20261024_JA.md)。この文書の植物誌直立/下垂U/Dと二値Mk解析・D型追加は**現在の主要推定対象ではありません**。過去の計算と出所を保存するために残しているだけです。S1の元170,066枚とS2の2種の連続写真角度が主形質で、T3も同じ実測方法が揃うまで延期します。


**PR #144 S1 → S2 → T2 → T3 の主線。** この文書は実際に動く trait-history パイロットの入力契約を定める。**実測の祖先状態・独立起源数・年代は未推定。** S1の +0.4382°/100 mm は画像縦基準の種間相関であり、T2/T3の形質真値や進化的事前確率に使用しない。V3の8段階・16参照値とS1/S2の元データは不変。

## 作業と停止条件

1. **T2-A 植物誌独立採点:** 対象は開花中の頭花（capitulum）。直立 erect、点頭/下垂 nodding、記述上の多型 variable、判断不能 unknown の4状態。必ず文献根拠引用、節/ページ、原学名、accepted concept IDを保持。植物学的評価者2名を**気候・系統・現在の効果結果から盲検化**。不一致は第三者裁定まで停止。矛盾する異なる文献は多数決で真値化せず停止。蕾と開花後は別タグ。
2. **T2-B 独立系統:** Cirsium核マルチ座位樹の実ファイル、出所 DOI、SHA256、tip↔accepted concept 一対一対応、根の根拠、枝長単位を要求。**52本の既存grafted-placement感度樹は禁止。** 枝支持・交雑/倍数性・核樹代替を監査し、未解決tipは仮に置換しない。年代には独立の校正データと不確実性が必要。
3. **T2-C 条件付き祖先推定（実装済み）:** 二値erect/noddingのみを使い、各監査済み樹で対称二状態Mk尤度・最尤速度を推定し、枝内の複数変化を許す条件付き確率マッピングで推定遷移分布を出す。variable/unknownおよび不明な写真角度を0/1に置換しない。出力区間は**樹と速度点推定に条件付き**であり、系統/年代/測定の全不確実性を含まない。
4. **T2-D 雨×温度の進化モデル（未実行）:** BIO12とBIO1、標高、地域、空間、開花段階、系統の不確実性を同時考慮して初めて「湿潤下で独立に下向き化したか」を検証する。単純なtipの雨分布、T2-Cの遷移数、52感度樹の広がりを雨よけ選択と解釈しない。BIO15は二次。
5. **T3 別連への拡張（入力ゲートのみ実装）:** Cirsium → Cardueae他属 → Cremanthodiumを含む別の連の順。**開花中の頭花自体の同等な姿勢、分類と系統**をレビューし、2連以上あっても、それだけで収斂と呼ばず、根と各クレードの祖先状態・独立遷移・雨×気温の比較モデルを求める。不可ならT3を止め、T2を弱めない。

## 利用候補の公開データ（存在確認 ≠ tree-ingested）

| スコープ | 実在する参照 | 未充足 |
|---|---|---|
| Cirsium核種系統 | Moreyra et al. 2025 DOI 10.1016/j.ympev.2025.108285（350核座位、248 Cirsium種、SRA PRJNA957074） | Newickの厳密なファイル、tip照合、枝長単位/根/年代校正がPRに未登録 |
| Cardueae核基幹系統 | Herrando-Moraira et al. 2019 DOI 10.1016/j.ympev.2019.05.001; データ DOI 10.17632/bhvv6rmyt6.1（76種、樹・年代付録あり） | Cirsium属内の密な反復の代用にはしない。別連を含むAsteraceae全体の系統でもない |
| 他連の機構の先行例 | Chen et al. 2013 DOI 10.1080/17550874.2012.702793 （Cremanthodium campanulatum 操作実験） | Cirsiumの雨保護の因果推定や独立収斂の確証にはしない |

EAzamiには moreyra_2025_retrieved_samples.csv と build_azami_eazami_trait_tip_crosswalk.py の候補台帳が存在する。EAzami/data/phylogeny の topology-fragment CSV を本研究の実枝長樹として扱わない。

## PGLS既存52配置樹の再利用：T2-P（独立系統検証ではない）

**技術的には再利用できる。** 現行PGLS樹は GBOTB.extended.LCVP／V.PhyloMaker2 由来で、S1/S3と50回のランダムなS2配置を含む。枝長は存在するが、**アザミ属の実際の分子配列で各分岐を推定したものではない**。特に属内へ追加された種の相対的位置・枝長は推定シナリオに依存する。対象216分類群すべてがGBOTBに直接載っているわけではないため、52本を独立な52回の進化とみなさない。

[新しい探索的スクリプト](t2_pgls_placement_screen.py)は、既存のZenodo DOI 10.5281/zenodo.22703215内で SHA256 固定された**同じ52本**を復元し、Siniscalchi et al.（2023; DOI 10.1086/724310）の著者公開 HeadPosition（0 erect、1 pending）63分類群（Git固定コミット c925d4af86160064448f794a9715f70b917a2786、blob fa127f0587b676e6e004ea121424d26069586bcb）と**名前完全一致**した分類群だけを解析する。

- この63分類群は既報の**著者による独立した形態コード**である。ただし現行T2-Aが要求する新規の盲検二重採点ではなく、生育段階の個別検証も含まない。PGLSの写真から推定した角度や環境値を正解ラベルに転用しない。
- 各配置樹について、同じ共有tipの実際の状態数、根を自由・直立・下垂と置く最小遷移回数、観測支持が十分な場合だけ対称Mkの条件付き遷移シミュレーションを計算する。実行入口は12以上の共有tip、下垂2以上・直立2以上。**不足の場合は結果を偽装せず parsimony-only/blocked** と出力する。
- 52シナリオの差は「バックボーンへの配置の仮定を変えた場合の推定値の感度」であり、**独立起源数の信頼区間でも、絶対年代でもない**。被覆が少ない場合、検出力がないこと自体を記録する。元の植物誌全属採点・独立核種樹への接続は継続課題。
- 数値から「雨の環境下で下向き化が多かった」と解釈するには、BIO12とBIO1、標高、地域、共通祖先、開花段階およびカメラ姿勢を別途評価する必要がある。現スクリプトの数値は**雨の適応進化を検定しない**。
- T3の複数連への一般化に、属内のみの52系統樹を使わない。

**手元に既存PGLS treeがあるから捨てる、ではなく、探索的な祖先復元として使い、厳密なT2とは役割を分ける。**

## 入力ファイルと実行

分析は analysis/v4/rain_history.py。空ヘッダー雛形は analysis/v4/rain_history_templates に置く。

- **scores.csv** ヘッダー: record_id,accepted_taxon,concept_id,source_taxon,source_id,source_locator,evidence_quote,head_stage,annotator_id,annotation_role,blind_climate,blind_phylogeny,state。原則1根拠passageにprimary2名、必要に応じadjudicator1名。両blindフィールドは yes、head_stageは anthesis。
- **crosswalk.csv** ヘッダー: tip,accepted_taxon,concept_id,decision,tribe,phenotype_reference,head_stage,homology_reviewed。decision=verified_exact のみが状態付tips。T3では phenotype_reference=botanical_head_posture、homology_reviewed=yes、2連以上が必要。
- **tree.nwk と tree_manifest.json**: manifestに tree_sha256, source_doi, tree_type=independent_nuclear_species_tree, root_method, branch_unit=substitutions_per_site または million_years, review_status=independently_reviewed を記録。million_yearsには calibration_doi, calibration_method, calibration_uncertainty も必要。
- **climate.csv（任意）**: concept_id,chelsa_bio12,chelsa_bio01,elevation_m,region,sampling_frame。結合品質のみ計算し、この版では気候の回帰は実施しない。

コマンド:

    python -m analysis.v4.rain_history --stage T2 --scores scores.csv --crosswalk crosswalk.csv --tree tree.nwk --manifest tree_manifest.json --climate climate.csv --output result.json --draws 1000

T3は --stage T3 に変更する。既定20以上の採点済み二値tips、両状態各2以上を要する。synthetic testsは小さいサンプルで動作を検証するが実在の証拠として登録しない。出力ステータスは conditional_trait_only_Mk_pilot または blocked_no_historical_inference だけであり、どちらも**rain-adaptation effect verified**ではない。模型を通すために架空の樹・枝長・状態を作らない。

**実データ依存の未完了:** 全対象種の植物誌二重採点、ファイル実体とSHA256で確認した分子樹、分類照合、雨・気温同時の進化モデル、連間の独立起源の判定。完了前はT2/T3の進化回数・年代・雨の選択について数値主張しない。

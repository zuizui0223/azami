# GEB v2 figure and table consistency audit

Scientific source: current `main` at `692a5eb98efce6dc8b8a2ba49a3e0eacf92d46d7`, the frozen full-27 result ledger, `manuscript/final_claims.json`, and the machine-readable outputs under `analysis_outputs/v2_full27_environment_atlas_2026-08-27/`. The audit reshapes frozen outputs only; it adds no model, P value, selection rule or exploratory hypothesis.

## 1. 現在の図表構成の問題点

監査開始時のMainは4図で、旧Figure 1は模式的な箱だけだったため、実写真から数値までの測定過程を検証できなかった。地理的なサンプリング密度、フィルタリング、Europe + North Americaへの集中、実現気候領域もMainに無かった。旧Figure 4は最終2候補の棒だけを示し、234行から2行へ至る逐次ゲートを追跡できなかった。Main Tableは無く、Supporting Informationにも全結果を一括追跡できる投稿用表が無かった。

これらのうち、凍結current v2から再構成できる範囲は解消した。Mainは測定、実現地理、taxon-mean information loss、scale atlas、sequential robustnessの5図となり、Table 1が研究段階と図表を結び、Table 2が最終2候補の正確な証拠鎖を示す。Supplementは診断図4点をMainと重複させず、完全表S1–S10をワークブックに収録する。

ただし、Figure 2の地図は最終46,276観察の実現分布である。current v2の公開freezeには406,582件のdetector-positive source streamを座標付きで再構成できる表が含まれないため、source-stage geographyは流れ図の件数としてのみ表示した。旧6,626-head artifactの座標で代用するとv1 cohort混入になるので採用しなかった。detector-positive source geographyを別レイヤーで示すには、current v2 source-stage coordinate artifactをfreezeへ追加する必要がある。

## 2. 不足しているMain figures/tables

監査時に不足していたFigure 1、Figure 2、Tables 1–2を追加した。

- Figure 1は、open-licensed実写真、実bbox、tight/context crop、corolla mask、outline/convex hull、image-vertical orientation、代表的な連続値、27 registered / 22 measured / 5 unexecutedを一枚に示す。写真部分はpresentation/source provenance onlyであり、v1結果値を持ち込まない。
- Figure 2は、2-degree disclosure-safe cellへ集約した最終46,276観察の世界密度、406,582から46,276へ至る固定コホート流れ、92.26%の地域集中、BIO1–BIO12実現領域を示す。representativenessは主張しない。detector-positive source geographyは上記artifact欠落のため地図化せず、source-stage countのみである。
- Table 1は8解析段階を4列へ圧縮したportrait版で、範囲、統計操作・判定ゲート、Main/Supplement表示を対応付ける。完全な8列版はTable S10へ保持する。
- Table 2は最終2候補について、global among-taxon、sampling-composition、broad-space、52-tree historical-placementの凍結値を一行で追跡可能にする。

Mainは5図 + 2表 = 7 display piecesとなり、GEB Research Articleの目安である6–8点に収まる。凍結出力だけで追加可能なMain図表は閉じた。唯一の未充足要素は、Figure 2におけるdetector-positive source geographyの地理レイヤーであり、これはcurrent v2 source-stage coordinate artifact待ちである。Mainに完全係数表や全シナリオを移す必要はない。

## 3. 不足しているSupplement figures/tables

Supporting Figures S1–S3は役割が独立している。さらに、既存の104 tree-level rowsを新しい推論なしで表示するFigure S4を追加し、52-tree placement sensitivityの係数、P値、lambdaの分布を可視化した。欠けていた完全表も生成した。

- Table S1: 27-endpoint contract、27行。
- Table S2: cohort/geographic filtering、5行。
- Table S3: complete variance decomposition、27行。
- Table S4: complete within-taxon atlas、234行。
- Table S5: complete among-taxon atlas、468行。
- Table S6: complete cross-scale classification、234行。
- Table S7: sampling robustness、134のselected-row × scenario-family要約行、シナリオ数合計674。
- Table S8: spatial + historical sequential gate、36行。
- Table S9: secondary whole-capitulum frozen headline、2行。
- Table S10: full eight-column analytical roadmap、8行。

Table S9はResultsに127 taxa、18 endpoints、42.3%、Spearman 0.3663という数値主張があるため追加した。全表はCSVと書式付きXLSXの双方を保持する。

## 4. 解析段階 ↔ 図表対応表

| 解析段階 | source output | Main | Supplement | status |
|---|---|---|---|---|
| photo localization and continuous measurement | endpoint inventory + Figure 1 source provenance | Figure 1; Table 1 | Table S1 | sufficiently documented |
| detector-positive source geography | current v2にはsource-stage countのみ、座標付きfreezeなし | Figure 2のflow countのみ | Table S2 | figure missing |
| frozen cohort and realized domain | final claims + environment report + disclosure-safe cell aggregation | Figure 2; Table 1 | Table S2 | sufficiently documented |
| taxon-mean information loss | variance decomposition | Figure 3; Table 1 | Table S3 | sufficiently documented |
| within-taxon atlas | complete within table | Figure 4; Table 1 | Table S4 | sufficiently documented |
| among-taxon atlas | complete among table | Figure 4; Table 1 | Table S5 | sufficiently documented |
| 234-row scale classification | cross-scale table | Figure 4; Table 1 | Table S6 | sufficiently documented |
| 674 sampling scenarios | scenario + summary tables | Figure 5; Table 1 | Table S7; Figure S2 | sufficiently documented |
| broad-space gate | spatial within + among tables | Figure 5; Table 1 | Table S8; Figure S3 | sufficiently documented |
| 52-tree placement gate and final candidates | historical models + summary + final claims | Figure 5; Tables 1–2 | Table S8; Figure S4 | sufficiently documented |
| secondary whole-capitulum synthesis | final claims | none | Table S9 | sufficiently documented |

段落単位の完全対応は `Manuscript_claim_figure_table_mapping.csv` に保存した。

## 5. v1 contamination audit

アクティブな本文8ファイル、Supporting Information、図表生成script、DOCX生成scriptを検索した。`9 primary endpoints`、`6,626`、`3,725`、`36 primary CHELSA tests`、`SPDE-INLA`、旧PCAの`32.9%` / `23.2%` / `69.3%`は検出されなかった。旧BIO1/BIO12係数は再利用していない。

`50 randomized`はcurrent v2の52-tree sensitivityを構成する内訳としてのみ残る。必ずtwo deterministic + 50 randomized = 52 audited treesの文脈で記載し、旧50-tree結果とは扱わない。Figure 1の実写真と検出例は上流のproduction-photo provenanceで、caption、source metadata、figure build reportの三箇所でpresentation/source provenance onlyと明示した。Figure 2のNatural Earth地図は表示用basemapであり、推論入力ではない。

リポジトリ内の歴史文書にはv1語が残り得るが、投稿bundleは明示列挙されたactive source、current v2 outputs、current figures/tablesだけから構築される。validatorは上記v1固有tokenがactive submission sourceへ戻るとfailする。

## 6. 追加・削除・移動すべき図表の確定案

追加はMain Figures 1–2、Main Tables 1–2、Supporting Figure S4、Tables S1–S10。旧模式Figure 1は実写真pipelineへ置換した。旧Main Figures 2–4は内容を保ったままFigures 3–5へ再設計・改番し、旧ファイル名は削除した。Figure 5には234 → 10 → 2 → 2の逐次縮約と、各候補のbase、sampling minimum ratio、broad-space、residual Moran、52/52 treesを統合した。current v2 source-stage coordinate artifactが凍結された時点でのみFigure 2へdetector-positive geographyレイヤーを追加する。旧v1座標の転用はしない。

Supporting Figures S1–S4はSupportingに置く。Main図のSupplement再掲、系統図のMain復帰、v1 SPDE/PCA図の復活は行わない。完全係数、missing/unexecuted、sampling、spatial、historical、secondary synthesisはSupplement tablesに限定する。

GEB/Wiley要項に合わせ、図番号とタイトルは図画像内から削除し、図と同じページのstand-alone legendにのみ置いた。パネル表記は(a)、(b)、(c)へ統一し、PDF vector版と300 dpi PNG review版を保持する。全MainページはA4 portraitに統一し、旧landscape Table 1ページは廃止した。

## 7. 投稿時の最終Main/Supp figure/table numbering

- Main Figure 1: Measurement pipeline.
- Main Figure 2: Geographic sampling and analytical domain.
- Main Figure 3: Information loss by taxon means.
- Main Figure 4: Within/among scale atlas.
- Main Figure 5: Sequential robustness and final candidates.
- Main Table 1: Current v2 analytical roadmap.
- Main Table 2: Frozen evidence chain for the two candidate patterns.
- Supporting Figures S1–S4: endpoint support, full selected-row sampling audit, spatial diagnostic surface, historical-placement stability.
- Supporting Tables S1–S10: endpoint contract, cohort filtering, variance decomposition, within atlas, among atlas, cross-scale classification, sampling robustness, sequential spatial/historical gate, secondary synthesis, full analytical roadmap.

科学的statusは`submission_hold_pending_independent_validation`のまま。図表の追跡可能性は閉じたが、外部validation gateが閉じたことは意味しない。

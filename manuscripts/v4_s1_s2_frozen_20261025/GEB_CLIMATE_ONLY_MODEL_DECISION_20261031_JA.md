> **SUPERSEDED MODEL CHOICE — historical record (31 October 2026).** The choice below to use BIO18+BIO01+RSDS across all Results was subsequently reconsidered. The current author-directed recommendation after authenticated BIO05×RSDS evaluation is **BIO18+BIO05+RSDS+elevation for Results 3.1/3.2 and original BIO12+BIO01+RSDS for Results 3.3**. See [the current primary-model decision](GEB_PRIMARY_MODEL_AFTER_INTERACTION_AUDIT_20261031_JA.md). This old document remains for numerical/source provenance; **do not cite its earlier model selection as current**.

# GEB主解析の説明変数確定とResults 3.2の科学的結論（2026-10-31）

**確定した説明変数：CHELSA BIO18（最暖四半期の月平均降水量）＋BIO01（年平均気温）＋RSDS（短波日射量）。標高を主モデルに使わない。**

- 3meanの３変数をResults 3.1～3.3をつなぐ共通環境軸とする。
- BIO18＋BIO05＋RSDS（3hot）は温度時間尺度の感度分析、BIO18＋BIO01＋RSDS＋BIO04＋BIO15（5seasonal）は気温・降水季節性の感度分析としてSIに残す。
- 2026-10-24凍結の標高あり４変数結果を削除・改名・新モデルへの有意性の流用はしない。**これは解析結果を確認した後のモデル選定なので探索的であり、事前登録された独立の確証的証拠ではない。**
- 論文本体 `GEB_RESEARCH_ARTICLE.md` と既存図3・4・表2は2026-10-31現在、**旧４変数の数値を掲載している**。新しい３変数に切替える本文/図表/SI/要旨の改稿は、全体整合性を確認した別コミットで行う。研究データのZenodo既存候補の数値も黙って書き換えない。

## 実測で確定したResults 3.2：S2全８形質×２種×在来/導入×３環境変数

出典：元のS1/S2係数を認証後に実行した[Actions #41179235813](https://github.com/zuizui0223/azami/actions/runs/41179235813)、artifact `13034028217`。新S2 familyは**96検定、Holm補正後２件のみ支持**。今回は全て標高を説明変数・写真の欠測除外条件に使っていない。

| 種・分布域 | 関係 | 対象写真・地理セル | 推定と有意性 | 結論 |
| --- | --- | --- | --- | --- |
| *C. vulgare* 導入域 | RSDS×可視画像花色明度 | **27,554写真・88セル** | 標準化β **−0.129204**；点別95%CI **[−0.187539, −0.070869]**；未補正P `1.41814e-5`；**Holm P=0.00136142** | 日射の高い地点ほど、花が画像上で暗く写る |
| *C. arvense* 在来域 | RSDS×円周色相（sin/cosを１組として２自由度共同検定） | **23,130写真・186セル** | χ²(2)=**16.83377**；未補正P `0.000221102`；**Holm P=0.0210047** | 日射と色相の間に関連がある。ただし色相変化の方向はこの共同検定だけでは特定できない |

**解釈上の重要な非支持点：**

- *C. vulgare* 導入域のRSDS×彩度：β **+0.081303**、Holm **P=1.0**。従来の「日射↑で鮮やかかつ暗い」を主結果として残せない。「暗い」だけが支持された。
- *C. arvense* 導入域のBIO01×明度：β **+0.090168**、未補正P=0.000669、BH P=0.021406、**Holm P=0.062881**で未支持。元の標高×明度は主モデルに標高がないので主結果ではない。
- *C. vulgare* 在来域のRSDS×明度：β **−0.044582**、Holm **P=1.0**。**導入域と在来域の傾きの差そのものは検定していない**ので、導入後の進化・導入域特異的な応答を実証したとは言えない。
- 同familyのBIO18（降水）またはBIO01（気温）に対する花色・輪郭等の項目は、Holm 5%基準を通過していない。２件の支持関連はいずれも日射（RSDS）との関係である。

### 3.2の本文に採用する結論案

**「同種内の花色に見られる気候との地理的な関係は、２種の異なる分布群で日射量と対応した。ただし、導入域の *C. vulgare* では明度、在来域の *C. arvense* では円周色相が関係し、普遍的な暗色化・高彩度化や導入後の適応進化を示すものではない。」**

短波日射は花色の可視写真指標と統計的に関連したのであり、UV吸収、熱吸収、色素量、適応度を測ったものではない。これらは今後の機能仮説と位置づける。

## 他のResultsとのつながり・投稿時の処理

- **3.1（S1）** 新３変数のBIO18×画像頭花角度：**169,952写真／178分類群、β=+0.046221、Holm P=0.001132**。地理・系統交絡の問題は残る。
- **3.2（S2）** 上記２件。４分布群の比較は独立の野外実験ではなく、分布域×環境交互作用は未検定。
- **3.3（H3c）** BIO18・BIO01・RSDSでは、1°および5°双方で分類群を超えた局所的な形態類似を確認（別の探索的Holm-4 P=0.020）。季節性２変数を追加しても形態類似の増分は未支持。
- 新３変数モデルのテストファミリーと異なる探索的Holm-4を混合して、全論文の単一補正Pと主張しない。本文・Figure 3/4・Table 2・SI・アブストラクトの改稿を同期させる。

[全８形質の元データ記録](../../analysis/v4/results/20261031_S1_S2_CLIMATE_ONLY_3_VS_5_SOURCE_VERIFIED_JA.md) · [機械可読な元データ集計](../../analysis/v4/results/20261031_S1_S2_CLIMATE_ONLY_3_VS_5_SUMMARY.json) · [H3c４モデル感度](../../analysis/v4/results/20261031_H3C_FOUR_ELEVATION_FREE_CLIMATE_AXES_SOURCE_VERIFIED_JA.md)

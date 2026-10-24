# T3｜S1連続画像角度を別のキク科系統へ展開するための実写真ソース監査（2026-10-24）

> **2026-10-24 UPDATED S1→S2→T2→T3 SCIENCE:** この15分類群の小規模初回結果は新しい[全131核ITS原名／48画像分類群／720四気候完全写真／115画像の盲検QC台帳](../S1_S2_T2_T3_FOUR_AXIS_ECOLOGICAL_MAINLINE_20261024_JA.md)の下位資料。下記の15分類群をT3本体全データと誤解しない。現在のT3にS1と等価な実測画像角度は0件で、BIO18などの形質検定は未実行。


**元131分類群の真正核rDNA配列名に一致する15分類群を、形態文献U/Dや環境値を見ずに選び、実際のGBIFを照会しました。**

- [実データActions #40057656854](https://github.com/zuizui0223/azami/actions/runs/40057656854) **SUCCESS、ソース照合・1テスト・15分類群GBIF照会完走**。
- [実行コード](../t3_s1_continuous_external_photo_source_pilot.py)・[CI](../../../.github/workflows/v4-t3-s1-continuous-external-photo-source-pilot.yml)・[検査](../../../tests/test_v4_t3_s1_continuous_external_photo_source_pilot.py)。
- 真正ITS核rDNA原配列の131-taxon照合元：公開2025配列取得Actions #39960458391、artifact 12504855926、ZIP SHA256 `64419785872135e4232eae8f96d6976f22da958300764209619951223fd9a55e`。新GBIF成果はartifact **12545992902**、ZIP sha256 **`b3034f9c2617b823a4fd61b958ad054fe302c11421c5ca31f45d548ea68a5a16`**。内包 `T3_PUBLIC_IMAGE_SOURCE_ELIGIBILITY.json` sha256 `6526cb09472cfc3871c645f92351d1bed069dae40894293be8c3f82b1613da56`、126写真URL・権利者・ライセンス台帳 `T3_PUBLIC_IMAGES_LICENSED_SOURCE_CANDIDATES.csv` sha256 `db1f95c59e7f83f937b3b21dbc0ed72114c2090f849337f10fcc3ecd3bb45539`。

## 取得状況

厳密GBIF学名照合成功**15/15**。これらの先頭30件までの画像付き・座標付き・HUMAN_OBSERVATIONの公的検索から、**7分類群の126個の実写真付き観察**を認証。0件は**先頭30件・指定GBIF絞り込みで検出なし**という意味であり、ネット上に写真が存在しない証明ではない。

| 元131核rDNAの選択種 | 元GBIFで認証した画像記録 |
|---|---:|
| *Cremanthodium brunneopilosum* | 5 |
| *Cremanthodium campanulatum* | 4 |
| *Cremanthodium decaisnei* | 5 |
| *Parasenecio adenostyloides* | 22 |
| *Senecio scandens* | 30（検索上限） |
| *Tussilago farfara* | 30（検索上限） |
| *Crassocephalum crepidioides* | 30（検索上限） |
| **合計** | **126** |

残り8原名は、初回の条件付き照会では0件。対象は*Ligularia atroviolacea*、*L. botryodes*、*L. caloxantha*、*Parasenecio ainsliaeiflorus*、*P. albus*、*Sinosenecio baojingensis*、*S. bodinieri*、*S. changii*。これらを生物学的に「写真なし」と断定しない。

写真のGBIFライセンス欄はCC BY-NC 4.0 **111件**、CC BY 4.0 **15件**。画像本体は**一切取得・再公開していない**。元URL・権利者・ライセンス・元GBIF occurrence keyを保持し、ライセンスの条件に従った取得／画像処理判断を後工程に分離。

## 2026-10-24 続報：S1オリジナル検出器の長期アーカイブ回収と真正性検査も完了

旧GitHub Actions artifact 8076736948は404だったが、既存の[再現性アーカイブ台帳](../../../reproducibility/durable_archive_manifest.json)で**同一artifact ID 8076736948のGoogle Drive長期保存ZIP**を特定し、元ZIPを実取得した。**ZIP SHA256 `2bc7cc49c2f213d4a5c5dab96e42fcdda7ceda23ffc3130520a5dd44b8c76c05`** は台帳と完全一致。内部 `recovery/model/weights/best.pt` の**SHA256 `4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed`** も凍結台帳と完全一致（サイズ5,457,498 byte）。

これにより、**T3の外部分類群写真で検出器転移を評価するためのモデルファイルの所在地・真正性の問題は解消**した。ただし**モデルがCremanthodium等にも有効であるとは未検証**で、元GBIF写真126件の画像本体取得・同一撮影条件の監査・頭花ROI/PCAと花柄の検出・S1と等価な連続角度計測はいずれも未実施。モデルは既存の非共有アーカイブで保持し、公開Gitリポジトリへ無断転載しない。

## 解析を「実行した」と言わない境界

**T3の連続角度測定件数はまだゼロ。** 元S1で計測した表現型は `orientation_image_vertical_angle`（EXIF補正画像縦軸への0–180°頭花軸角度）。今回の画像は**元S1のYOLO検出器／頭花と花柄の画像幾何測定をまだ通していない**。文献姿勢D9/U0を新写真の角度やラベルに変換しない。

同一表現型でT3へ進むための実作業：
1. 現在の元S1の凍結したcapitulum detector＋幾何測定スクリプトにアクセスし、別属15分類群のうち写真が取得できる7分類群で、**頭花が写っているか、花柄が見えるか、角度の成否、撮影姿勢・開花期**を監査する。元モデルは*アザミ属*特化なので別属へ即転移できると仮定しない。
2. 元検出器の[凍結配布台帳](../../../reproducibility/actions_artifact_catalog.json)は旧GitHub artifact **8076736948**・モデルSHA256 `4078e0510532852681b65ee529cd82237b649ec99b17c4ca5f1da460a62d2bed`を記録。2026-10-24のGitHub artifact直接取得は404（保存期限切れ）。台帳の**Google Drive長期アーカイブ**に回収記録があるので、そちらから原実体を検証してはじめて使用する。新しいモデルの性能として誤記しない。
3. S1との測定不変性／共通花期／種間撮影構成の感度が成立したとき、2連以上（CirsiumのCardueae、今回Senecioneae）の**年降水量と連続画像角度の種間・種内/群内傾向**を同じ推定量で比較する。その後、真正核種樹＋地理補正を導入する。**今回の131核rDNAは単一の連鎖マーカー遺伝子樹**で、真正の多座位核種樹ではない。

**今回の進捗は：T3が写真入手の段階へ到達。実角度の推定・環境との対応や収斂進化を新たに確認した結果ではない。**

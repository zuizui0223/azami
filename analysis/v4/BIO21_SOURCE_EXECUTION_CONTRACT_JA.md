# V4-RF 12→21 実ソース検証ゲート

**作業境界:** PR #167 のみ。現行 GEB V4 main manuscript、S1 80検定・S2 128検定、H3c、元係数、表・図は不変更。**このファイルは実測結果ではない。**

## 何を検証するか

元 S1 の checksum-verified 観測表 313,350 件を GitHub Actions の保護ジョブ内のみで復元する。元の5変数デコード表と既存の追加7変数を別々にチェックし、CHELSA v2.1 1981–2010 の未使用 BIO02、BIO03、BIO07、BIO08、BIO09、BIO10、BIO11、BIO16、BIO19 を元観測点へ直接ピクセル抽出する。**BIO03 は isothermality（%）**、BIO16/19 は四半期の平均月降水量であり、他の温度形質と単位を取り違えない。公式原ラスタの source URL・実ダウンロード bytes・SHA-256・座標参照系・解像度・raw scale/offset を収集時に監査する（事前SHA確定とは区別）。

9変数を成功裏に抽出できた場合、**同一の環境完全ケース・同一の分類群×5°セル・同一のfold**で12対21変数のRFを比較する。対象は第一段階として `angle`, `chroma`, `aspect` の3構成形質。モデルは固定設定80 trees、min leaf 5、max features .75、held-out taxa GroupKFold（5分割）と longitude-sector LeaveOneGroupOut（3分割）。4→12と12→21は別の完全ケース集合になるため、差を横断比較しない。残り5形質は第一段階のデータ支持を確認してから実行する。

**判定:** `skill = 1 − RF weighted SSE / train-fold-mean weighted SSE`。正値でも因果・適応・気候追随・遺伝的形態進化の証拠ではない。負値、種別のデータ不足、ホールドアウトでの失敗も必ず報告する。温度・降水・放射・標高の相関グループに対応した置換重要度を併記するが、独立因果効果とは呼ばない。追加BIOCLIMには他の列と数学的・強相関の座標が含まれる。

## 実行と監査

- [保護ソース実測 workflow](../../.github/workflows/v4-rf-bio21-source.yml) は **PR #167 の同一所有者・同一ブランチ**に限る。鍵はステップのみに渡され、S1をAES-GCMで復号し、元のSHA receiptと照合する。
- [原CHELSAバイトと実観測ピクセル抽出](rf_bio21_source_real.py) は16層（既存7＋新9）を1層ずつ取得・照合・抽出し、ラスタを削除する。ファイル改変や空間グリッド不一致・BIO03等の物理単位逸脱があれば停止する。外部層の初回SHAは**取得した実バイトから計算**する。初回成功後にSHAを固定して再実行する。
- すべての非公開 `obs_id`、元座標、種名、復号表、新環境表、予測行、ラスタ、モデルweightは `work/rf21/private` と `work/rf21/temp` のみに置く。成功・失敗どちらも終了時に消去し、公開 Actions artifact は **`bio21_source_receipt.json` と `rf_12_vs_21_anchor3_aggregate.json` だけ**。
- GitHub上の「CI成功」は**実データのpositive findingの同義ではない**。2つの出力JSONに source-normalise、3形質ごとの持越し可能性、変数群寄与、欠測による支持コホート低下が全て含まれて初めて科学的判定する。

CHELSA公式：<https://www.chelsa-climate.org/datasets/chelsa_bioclim>。URLの候補は既存ソースのV2.1経路を用いる。ユーザに別個の非公開ファイルを要求しないで、既存保護アーカイブの利用可否を最初に確認する。

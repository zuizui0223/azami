# V4｜BIO1・BIO15とVPDは同じ「湿潤機構」なのか？（2026-10-24）

## 現在の結論

**VPDは重要な生理・微気候の代替仮説だが、既存解析で「BIO1＋BIO15が強い理由はVPDだ」と確認できていない。** 特に、CHELSAのVPDは飽和水蒸気圧と実際の水蒸気圧の差（大気の乾燥要求、Pa）であり、雨量そのもの、BIO15の降水月変動係数、頭花の雨濡れ時間とは別の物理量である。

- **EAzami H1（固定済み東アジア9/12）：** 植物誌での真正U→D遷移重み×現生BIO15上昇＋BIO1低下の複合対応（厳格9種4/126、12種16/792）。その解析でVPDを置き換えた同じtip・同じ内部枝・同じ帰無分布の**直接比較は実施済みではない**。後からVPDを加えてH1の有意性を強化したことにしない。
- **EAzamiの別の真正common9環境アトラス：** 原著公表Actions #34208671960の元`taxa.csv`・`univariate_rows.csv`を含むartifact **10049165416**（元ZIP sha256 `b92694b1527fe46856bda972420802dbe7468588d347c10d0dc02e74b5eae502`）を2026-10-24に取得。これは**日本出所17分類群（D5/U12、1種3点以上の間引き済み出現）であり、EAzami H1の東アジア厳格9分類群とは異なる標本集合**。この違いを忘れてVPDをH1に帰属させない。
- 17分類群の公表元数値表から再集計すると、**年平均BIO1とVPDのPearson r=+0.9728（17分類群）**、**BIO15とVPDのr=+0.7914**。この標本集団では大気乾燥度が温度と強く共変動しているという**共線性の発見**であり、形質との原因関係・媒介効果の検証ではない。
- 同じ17分類群の平均VPDは植物誌D型 **498.030 Pa（n5）**、U型 **528.825 Pa（n12）**。**D−U = −30.795 Pa = −0.344785 SD**。著者による全6188個のU/D数固定ラベルマップの**VPD二側rank 3287/6188=0.53119**、27形質×環境のBH q=**0.89638**。**下垂型が低VPD側という方向は見えるが、統計的な支持は弱い**。
- 対照としてこの17群では**BIO15のD−Uも−4.151（元値）**。これはEAzami東アジア厳格9群の「BIO15高い方にD」という向きと異なる。母集団・地点の出所・形質採点が異なるので、同じ仮説の直接矛盾・外部再現とは呼ばない。
- **azami V3の元9変数環境アトラス**（S1オリジナル5環境モデルとは別）：142分類群の種間の連続画像角度×VPDは `beta_std=-0.119421`、`p=0.1519`、`q=0.31793`（非支持）。一方種間BIO12は `beta_std=+0.304359`、`q=0.0036`。同種内の画像角度×VPDは37,196観察／199分類群で`beta_std=+0.014165`、単独P=0.04466、`q=0.15178`（多重調整後非支持）。**種間と種内の向きまで一貫したVPD効果は未証明**。
- **S2:** *C. arvense* と *C. vulgare* の元在来／全導入域のモデルはBIO1/BIO4/BIO12/BIO15/RSDSの5軸で、VPDはそのままでは含まれない。**導入後の角度分布差をVPDで説明できたとはまだ言えない**。

## どういう機構を識別すべきか

| 仮説 | 予測すべき直接曝露（現地で測る量） | BIO1/BIO15/VPDの限界 |
|---|---|---|
| **雨濡れからの生殖器官の保護** | **開花期間の降雨頻度・強度、頭花の花粉・柱頭の濡れ時間、結実** | 年平均VPDが低くても降雨回数・花期の濡れ時間とは同一ではない |
| **大気の乾燥・蒸散・微気候の調節** | **開花中の昼間VPD・花面温度・花冠／総苞の水収支、頭花周辺の湿度** | BIO1が気候背景としてVPDと高相関で、統計に双方を足すだけでは独立原因と解釈できない |
| **乾季／雨季と花期の一致** | **花期月の月降水と月VPDの同時系列、乾湿季の位相差** | BIO15は12か月降水の変動係数なので、**どの季節に開花するか／乾季が開花期か**を含まない |
| **送粉者または花食者を介した差** | 訪花・柱頭接触・花粉輸送・花食被害とその際の気象 | CHELSAの年平均VPDから相互作用の実時間の環境曝露は同定できない |

**現段階で最も識別的な次の比較（確認的H1とは別枠の事後探索）：**

1. **同一tip・同一出所の種別環境値だけ**で、EAzami厳格9/12のBIO1・BIO15・BIO12・VPD・開花期降水を揃える。ただしVPDが欠測なら欠測と示し、既存common9の17群をH1の9群の代用としない。
2. 別々の機構を比較する：**（A）花期の実降雨・濡れ曝露**、**（B）花期のVPD（乾燥要求）**、**（C）両者の組合せ**。CHELSA月別pr/VPDは長期平均の代理で、開花当日の降雨や湿度そのものではない。
3. S1の種間／種内とS2の在来／導入分布に対し、**同一元レコード・空間セル・写真品質の層**でVPDおよび降水指標を追加できたときだけ、新規の環境対照を評価。EXIF画像角度＝実際の下垂姿勢と置かない。
4. **系統的に独立なD/U反復**を持つ複数の核種系統で同じ比較を行う。現状T2真正核NA44 D3・Japan12 D5は探索段階。原データの生物学的有効反復数を増やす必要がある。

## 参照と再現

- EAzamiの公表元証拠：[Common9結果](https://github.com/zuizui0223/EAzami/blob/main/docs/chapter2/THREE_TRAIT_COMMON9_RESULT_V1.md)／[source JSON](https://github.com/zuizui0223/EAzami/blob/main/data/evidence/chapter2_three_trait_common9_result_v1.json)。GitHub Actions **34208671960**, artifact **10049165416** に`taxa.csv`・`univariate_rows.csv`を同梱。上記VPD平均・相関は**この同一17種元表**を用いた記述再集計であり、新たな生物学的独立実験ではない。
- Azami元結果：[生物学的な向き×VPDの種間表](https://github.com/zuizui0223/azami/blob/main/reproducibility/current_reference/axes/biological_axes_among_min5.csv)／[種内表](https://github.com/zuizui0223/azami/blob/main/reproducibility/current_reference/axes/biological_axes_within.csv)。
- EAzami主結果：[H1の複合ニッチ対応](https://github.com/zuizui0223/EAzami/blob/main/docs/chapter2/CURRENT_CLAIMS_AND_EVIDENCE_TRANSITION_REGIME_V2.md)／[BIO15の独立性未確証](https://github.com/zuizui0223/EAzami/blob/main/docs/chapter2/BIO15_SPECIFICITY_MATCHED_CONTRAST_V1.md)。
- 外部定義：[CHELSA V2.1 climatologies](https://www.chelsa-climate.org/datasets/chelsa_climatologies)のVPD単位Pa、[NOAAのVPD](https://psl.noaa.gov/forecasts/seasonal_vpd/description/)の`e_s-e_a`。
- この追補を追加してもEAzamiの凍結H1やV4の5環境モデルの数値は更新しない。VPDは**代替機構の鑑別**であり、事後に有意な予測子を探して主仮説を強める操作ではない。

# 結果3：気候ニッチ拡大に伴う種内表現型分化の方向・多様性・新規性
2026-11-26。彩度・頭花の画像内角度・投影輪郭縦横比の3軸のみ。ユーザーの指定により明度・色相・時間軸は除外。従来の8形質探索は履歴として残し、主解析に採用しない。

## 1. 気候ニッチ：分化が起こり得る環境差
背景補正済み年間5気候による在来・移入の占有D（arvense/vulgare）=.474/.442。移入密度の在来主要95%域外=18.8%/44.0%。操作的な現代の占有分布で、祖先気候や生理限界、経過時間を測っていない。
全適格セルの3形質分布との在来重複Dは移入類似 .894/.917、新規 .732/.845（等セル観測重み）。地域マッチ .714→.732/.864→.817、2種Holm=.559/.124。気候Dと形質Dは異なる空間の重複指標なので、差を直接的な効果量としない。
出典：[whole range](../whole_range_overlap/RESULTS_JA.md)、[record overlap](../record_overlap/RESULTS_JA.md)。

## 2. 方向的分化：どの形質へ偏るか
固定95%気候区分、地域マッチ1°セルの中央値の差（新規−類似）は C. arvense: 彩度−0.136、角度−7.708°、輪郭+.009239、C. vulgare: 彩度+4.203、角度+2.511°、輪郭−.000331。候補の角度と彩度は未補正P=.03906/.03638、6比較Holm=.21826。新規気候に伴う形質の方向的差は探索的。等セル観測重みでも arvense 角度−7.701°、vulgare 彩度+3.116。密度超過仕様のvulgare彩度+2.425。以前の希少度関連は地域内調整で弱まり、少雨や季節性が分化を引き起こすとの結論は採用しない。
出典：[matched centroids](../configuration/matched_centroids.csv)、[observation audit](../phenotype_observation_audit/RESULTS_JA.md)。

## 3. 既存の表現型区分の再配分 vs 在来未観測区分の増加
在来256参照セル、独立地理128在来検証セル、32抽出。観測重み、5分位×3＋各軸範囲外。95%主仕様の類似→新規の正の頻度変化のうち既存区分 E/(E+U) は87.5214%/82.9304%、固定ブロック置換帰無平均84.6012%/81.8107%。全32 Holmで既存区分が優先される支持なし。未観測区分の全頻度は在来検証13.5986/13.9893%、移入類似15.2247/13.9745%、新規11.5841/13.9810%、差−3.6406/+.0065 pp、Holm=1。全分布の頻度再配分TVはvulgareで帰無より大きい候補（raw P=.00268555、Holm=.0859375）。既存区分の寄与割合が高いことは在来区分占有率が高いことで説明し得る。未観測区分は新規突然変異と同義ではない。等価性も未証明。
出典：[3軸分解](../phenotype_three_axis/RESULTS_JA.md)、[統合仮説台帳と帰無](../RESULT3_THREE_AXIS_FINAL_JA.md)。

## 4. 新規追加：分布の多様化 vs 収縮
[Run 45781091217](https://github.com/zuizui0223/azami/actions/runs/45781091217) success. commit a6e8348573a110175fa2d3ec9b4adfa92074477e; artifact 15098068622. ZIP実体のSHA256 5970f107e0d5e061453678259b7c65eb1cebb528e7ae44c58930f53074a84dfa が公開digestと一致。dispersion_groups.csv 20行、dispersion_tests.csv 8行。全分散=セル間+セル内（最大誤差6.67e-16）、3軸別分散和=全分散（8.89e-16）、出力の差と比を各群の集計値から独立再計算、Holm再計算誤差0。全8検定Holm>=.05。全在来等セル観測平均・SDで各軸を固定標準化、3軸の分散の和trace（幾何学的体積ではない）。在来全分散3。独立標本は地理セルで、同じセル内の記録数に比例した重みではない。

### 全域・全適格セルの記述
| 種 | 群 | セル/観測 | 全分散 | セル間 | セル内 |
|---|---|---:|---:|---:|---:|
| arvense | 在来 | 1366/18505 | 3.000000 | 1.328019 | 1.671981 |
| | 移入全体 | 860/22269 | 2.607760 | .807287 | 1.800473 |
| | 移入類似 | 823/22060 | 2.630606 | .816239 | 1.814367 |
| | 移入新規 | 30/167 | 1.908714 | .619830 | 1.288884 |
| vulgare | 在来 | 1113/13115 | 3.000000 | 1.215312 | 1.784688 |
| | 移入全体 | 1139/25532 | 3.037481 | 1.026039 | 2.011442 |
| | 移入類似 | 931/22314 | 3.055581 | .970332 | 2.085249 |
| | 移入新規 | 190/3158 | 2.985649 | 1.272427 | 1.713223 |

移入全体の全分散は在来より−13.075%/+1.249%、セル間−39.211%/−15.574%、セル内+7.685%/+12.706%。全域差は地域構成が違い、検定していない。

### 地理マッチ（類似→新規）
| 種 | 定義 | 対/地域 | 指標 | 類似→新規 | 差% | 置換P | Holm全8 |
|---|---|---:|---|---:|---:|---:|---:|
| arvense | 95% | 30/9 | 全 | 3.378296→1.908714 | −43.501 | .08203125 | .65625 |
| | | | セル間 | 1.539241→.619830 | −59.731 | .296875 | 1 |
| vulgare | 95% | 146/14 | 全 | 2.900574→3.012473 | +3.858 | .573730 | 1 |
| | | | セル間 | 1.037998→1.318752 | +27.048 | .097290 | .681030 |
| arvense | 密度超過 | 183/14 | 全 | 2.486342→2.446593 | −1.599 | .924805 | 1 |
| | | | セル間 | .933115→.794143 | −14.893 | .938354 | 1 |
| vulgare | 密度超過 | 316/15 | 全 | 3.078864→3.206469 | +4.145 | .129089 | .774536 |
| | | | セル間 | .953731→.878087 | −7.931 | .851624 | 1 |

全地域ラベル反転、2側検定。arvense主仕様全分散差の点別ブロックbootstrap95%区間[-2.309465,-.287306]は0を含まないが、置換P=.082、Holm=.656。推論手順の違いを維持し、補正後の支持を主判断とする。95%主仕様のセル内分散−29.916%/−9.066%（記述のみ、別検定なし）。密度超過では+6.388%/+9.564%。arvense主仕様の全分散差−1.469583のうち輪郭分散−1.234682（約84.0%）、彩度−.230674、角度−.004227。画像輪郭が回転に依存し、定義変更で収縮幅が小さくなるため、実物の形態多様性の収縮を支持したとは言えない。vulgareは全分散+4%前後が2定義で同方向だが、セル間成分の符号は反転する。地理的な多様化を確証したとは言えない。各軸の標準化や分解は撮影誤差・地理的交絡を除去しない。

## 5. 採用する結果3の構成
A 気候占有差 → B 3形質の重複と方向 → C 既存/未観測の頻度変化と帰無 → D 多変量分散で分布中心の変化と多様化を区別。ポスターはA、B、Cを中心にし、Dを小さく併記するか本文・補足へ。観測値は「新しい気候への進出が形質空間の一様な拡大を伴わない」と整合する。vulgareの高彩度側への偏りは探索的、arvenseの上向き側への偏りも同様。自然選択により既存変異が優先された、分化・放散が急速に進んだ、進化速度や適応度が増したとは主張しない。形質多様性が不変という等価性も未検証。短期の種内進化を動機にすることは可能だが、導入年や歴史的祖先を測らない横断比較。系統分岐を伴うadaptive radiationと区別する。

### Results draft (English; exploratory)
The introduced ranges occupied climatic conditions outside the core native climatic niche (18.8% and 44.0% of background-corrected introduced occupancy density in C. arvense and C. vulgare, respectively), while their three-dimensional capitulum phenotype distributions broadly overlapped those of the native ranges. In geographically matched comparisons, novel-climate cells showed tendencies toward more upward-facing capitula in C. arvense and higher chroma in C. vulgare, but neither contrast survived multiple-testing correction. Most positive frequency differences between matched analogous- and novel-climate groups occurred in phenotype bins observed in the native reference (87.5% and 82.9%); these shares were not greater than expected under regional label exchange. Neither the frequency of native-unobserved phenotype bins nor multivariate dispersion showed a significant increase after correction. The descriptive results are compatible with occupancy of new climates using phenotypes broadly represented in the native range, without evidence of a general increase in phenotypic disparity. They do not establish rapid evolutionary change, genetic novelty, preferential selection on standing variation, or phenotypic equivalence.

### 主解析の終了
新規モデルを追加せず、上記の設計を固定。結果3のデータ・方法・図表候補を統一。main・既存発表資料は未変更。今後の実測の遺伝性・機能・画像精度評価は別課題であり、今回の有意性を求めた再探索は行わない。

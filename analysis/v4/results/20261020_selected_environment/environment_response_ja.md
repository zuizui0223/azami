# 頭花形質の環境勾配応答

問い：他の環境条件を同時に考慮したとき、どの環境勾配に沿って頭花形質が変わるか。

6群の全9候補から環境データのみでvifstep選定し、採択セットを固定して同時投入：温度季節性、降水、降水季節性、日射、VPD、風速、生育期降水量、NPP。観測値を平均せず、環境を種平均と種内偏差に分け、成分ごとに種ランダム切片を推定した。色全体は明度・彩度・円周色相の共同検定。各成分の利用可能な観測を使う。

過去の結果を閲覧済みの改訂観測解析であり、事前登録した独立検証ではない。48共同検定でHolm補正。成分の95%区間は点別である。

|形質|環境|種間P（Holm）|種内P（Holm）|
|---|---|---:|---:|
|色全体|温度季節性|0.6722|0.0004381|
|色全体|降水|2.977e-07|1|
|色全体|降水季節性|0.0656|0.8182|
|色全体|日射|1.022e-05|0.002617|
|色全体|VPD|0.0001368|0.004335|
|色全体|風速|1|4.035e-11|
|色全体|生育期降水量|6.69e-05|0.0001334|
|色全体|NPP|0.006973|0.4663|
|向き|温度季節性|1|1|
|向き|降水|0.02398|1|
|向き|降水季節性|1|1|
|向き|日射|1|0.2325|
|向き|VPD|1|0.09556|
|向き|風速|0.5096|1|
|向き|生育期降水量|1|1|
|向き|NPP|1|1|
|総苞片の輪郭突出|温度季節性|1|0.3428|
|総苞片の輪郭突出|降水|1|1|
|総苞片の輪郭突出|降水季節性|1|1|
|総苞片の輪郭突出|日射|1|1|
|総苞片の輪郭突出|VPD|0.7169|1|
|総苞片の輪郭突出|風速|1|1|
|総苞片の輪郭突出|生育期降水量|1|0.1274|
|総苞片の輪郭突出|NPP|1|0.113|

## 各成分の変化

全自生域母集団の環境1標準偏差当たりの変化。単独の成分検定で色全体の関連を判定しない。色相sin/cosは共同で読む。

|成分|観測数|種数|環境・尺度|変化［95%区間］|
|---|---:|---:|---|---|
|corolla_lab_lightness|186437|244|chelsa_bio04__between|+0.657 [-0.35, +1.66]|
|corolla_lab_lightness|186437|244|chelsa_bio04__within|-0.216 [-0.741, +0.309]|
|corolla_lab_lightness|186437|244|chelsa_bio12__between|-0.503 [-1.18, +0.169]|
|corolla_lab_lightness|186437|244|chelsa_bio12__within|-0.0903 [-0.5, +0.319]|
|corolla_lab_lightness|186437|244|chelsa_bio15__between|-1.57 [-2.32, -0.816]|
|corolla_lab_lightness|186437|244|chelsa_bio15__within|-0.383 [-0.841, +0.0757]|
|corolla_lab_lightness|186437|244|chelsa_rsds_mean__between|+2.03 [+0.502, +3.55]|
|corolla_lab_lightness|186437|244|chelsa_rsds_mean__within|-0.684 [-1.83, +0.462]|
|corolla_lab_lightness|186437|244|chelsa_vpd_mean__between|+0.193 [-0.886, +1.27]|
|corolla_lab_lightness|186437|244|chelsa_vpd_mean__within|+1.77 [+0.364, +3.17]|
|corolla_lab_lightness|186437|244|chelsa_sfcwind_mean__between|-0.00877 [-0.792, +0.774]|
|corolla_lab_lightness|186437|244|chelsa_sfcwind_mean__within|+0.0864 [-0.363, +0.536]|
|corolla_lab_lightness|186437|244|chelsa_gsp__between|+1.04 [-0.104, +2.19]|
|corolla_lab_lightness|186437|244|chelsa_gsp__within|-0.0126 [-0.0855, +0.0603]|
|corolla_lab_lightness|186437|244|chelsa_npp__between|+0.472 [-0.585, +1.53]|
|corolla_lab_lightness|186437|244|chelsa_npp__within|-0.252 [-1.24, +0.735]|
|corolla_lab_chroma|186437|244|chelsa_bio04__between|-0.542 [-1.7, +0.617]|
|corolla_lab_chroma|186437|244|chelsa_bio04__within|-0.661 [-1.5, +0.173]|
|corolla_lab_chroma|186437|244|chelsa_bio12__between|-1.25 [-1.74, -0.756]|
|corolla_lab_chroma|186437|244|chelsa_bio12__within|+0.119 [-0.265, +0.503]|
|corolla_lab_chroma|186437|244|chelsa_bio15__between|+0.961 [+0.0237, +1.9]|
|corolla_lab_chroma|186437|244|chelsa_bio15__within|+0.729 [+0.228, +1.23]|
|corolla_lab_chroma|186437|244|chelsa_rsds_mean__between|-4.42 [-6.11, -2.74]|
|corolla_lab_chroma|186437|244|chelsa_rsds_mean__within|-0.34 [-1.21, +0.531]|
|corolla_lab_chroma|186437|244|chelsa_vpd_mean__between|+0.143 [-1.08, +1.36]|
|corolla_lab_chroma|186437|244|chelsa_vpd_mean__within|-1.23 [-1.77, -0.691]|
|corolla_lab_chroma|186437|244|chelsa_sfcwind_mean__between|+0.0568 [-0.616, +0.73]|
|corolla_lab_chroma|186437|244|chelsa_sfcwind_mean__within|+0.173 [-0.158, +0.505]|
|corolla_lab_chroma|186437|244|chelsa_gsp__between|-0.197 [-1.16, +0.767]|
|corolla_lab_chroma|186437|244|chelsa_gsp__within|-0.132 [-0.182, -0.0824]|
|corolla_lab_chroma|186437|244|chelsa_npp__between|+1.8 [+0.661, +2.94]|
|corolla_lab_chroma|186437|244|chelsa_npp__within|+0.806 [-0.484, +2.1]|
|corolla_hue_sin|186388|244|chelsa_bio04__between|-0.0411 [-0.111, +0.0283]|
|corolla_hue_sin|186388|244|chelsa_bio04__within|+0.0846 [+0.0416, +0.128]|
|corolla_hue_sin|186388|244|chelsa_bio12__between|-0.0154 [-0.05, +0.0191]|
|corolla_hue_sin|186388|244|chelsa_bio12__within|-0.0184 [-0.0433, +0.0066]|
|corolla_hue_sin|186388|244|chelsa_bio15__between|-0.0548 [-0.108, -0.0012]|
|corolla_hue_sin|186388|244|chelsa_bio15__within|-0.034 [-0.121, +0.0532]|
|corolla_hue_sin|186388|244|chelsa_rsds_mean__between|+0.178 [+0.0776, +0.278]|
|corolla_hue_sin|186388|244|chelsa_rsds_mean__within|-0.043 [-0.0987, +0.0127]|
|corolla_hue_sin|186388|244|chelsa_vpd_mean__between|-0.145 [-0.218, -0.0725]|
|corolla_hue_sin|186388|244|chelsa_vpd_mean__within|+0.0393 [-0.0107, +0.0892]|
|corolla_hue_sin|186388|244|chelsa_sfcwind_mean__between|-0.000697 [-0.0459, +0.0445]|
|corolla_hue_sin|186388|244|chelsa_sfcwind_mean__within|+0.0119 [-0.00712, +0.0309]|
|corolla_hue_sin|186388|244|chelsa_gsp__between|-0.013 [-0.0789, +0.0529]|
|corolla_hue_sin|186388|244|chelsa_gsp__within|+0.00108 [-0.00308, +0.00524]|
|corolla_hue_sin|186388|244|chelsa_npp__between|-0.0856 [-0.157, -0.0143]|
|corolla_hue_sin|186388|244|chelsa_npp__within|+0.0266 [-0.00303, +0.0561]|
|corolla_hue_cos|186388|244|chelsa_bio04__between|-0.0188 [-0.0393, +0.00162]|
|corolla_hue_cos|186388|244|chelsa_bio04__within|+0.0076 [-0.0049, +0.0201]|
|corolla_hue_cos|186388|244|chelsa_bio12__between|+0.00117 [-0.0159, +0.0183]|
|corolla_hue_cos|186388|244|chelsa_bio12__within|+0.00171 [-0.00794, +0.0114]|
|corolla_hue_cos|186388|244|chelsa_bio15__between|+0.0176 [+0.000477, +0.0347]|
|corolla_hue_cos|186388|244|chelsa_bio15__within|-0.0016 [-0.016, +0.0128]|
|corolla_hue_cos|186388|244|chelsa_rsds_mean__between|+0.0289 [-0.000339, +0.0581]|
|corolla_hue_cos|186388|244|chelsa_rsds_mean__within|+0.0345 [+0.0146, +0.0545]|
|corolla_hue_cos|186388|244|chelsa_vpd_mean__between|-0.00252 [-0.0275, +0.0225]|
|corolla_hue_cos|186388|244|chelsa_vpd_mean__within|-0.0243 [-0.0478, -0.000817]|
|corolla_hue_cos|186388|244|chelsa_sfcwind_mean__between|-0.00676 [-0.025, +0.0115]|
|corolla_hue_cos|186388|244|chelsa_sfcwind_mean__within|+0.00907 [+0.00201, +0.0161]|
|corolla_hue_cos|186388|244|chelsa_gsp__between|-0.00903 [-0.023, +0.00491]|
|corolla_hue_cos|186388|244|chelsa_gsp__within|+0.00133 [-0.000285, +0.00295]|
|corolla_hue_cos|186388|244|chelsa_npp__between|-0.0199 [-0.0434, +0.00362]|
|corolla_hue_cos|186388|244|chelsa_npp__within|+0.00302 [-0.015, +0.021]|
|orientation_image_vertical_angle|170064|244|chelsa_bio04__between|-0.534 [-3.33, +2.27]|
|orientation_image_vertical_angle|170064|244|chelsa_bio04__within|+0.11 [-0.312, +0.531]|
|orientation_image_vertical_angle|170064|244|chelsa_bio12__between|+2.39 [+1.03, +3.76]|
|orientation_image_vertical_angle|170064|244|chelsa_bio12__within|-0.33 [-1.61, +0.948]|
|orientation_image_vertical_angle|170064|244|chelsa_bio15__between|+2.26 [-0.018, +4.55]|
|orientation_image_vertical_angle|170064|244|chelsa_bio15__within|-0.515 [-1.23, +0.198]|
|orientation_image_vertical_angle|170064|244|chelsa_rsds_mean__between|-0.39 [-3.86, +3.08]|
|orientation_image_vertical_angle|170064|244|chelsa_rsds_mean__within|-1.74 [-3, -0.479]|
|orientation_image_vertical_angle|170064|244|chelsa_vpd_mean__between|+0.0882 [-2.41, +2.58]|
|orientation_image_vertical_angle|170064|244|chelsa_vpd_mean__within|+2.76 [+0.971, +4.56]|
|orientation_image_vertical_angle|170064|244|chelsa_sfcwind_mean__between|+2.16 [+0.389, +3.92]|
|orientation_image_vertical_angle|170064|244|chelsa_sfcwind_mean__within|-0.124 [-0.565, +0.316]|
|orientation_image_vertical_angle|170064|244|chelsa_gsp__between|+0.306 [-0.814, +1.43]|
|orientation_image_vertical_angle|170064|244|chelsa_gsp__within|+0.0826 [-0.129, +0.294]|
|orientation_image_vertical_angle|170064|244|chelsa_npp__between|-0.708 [-3.32, +1.91]|
|orientation_image_vertical_angle|170064|244|chelsa_npp__within|-0.0317 [-0.979, +0.915]|
|bract_projection_p95|1303|111|chelsa_bio04__between|+0.00115 [-0.0114, +0.0137]|
|bract_projection_p95|1303|111|chelsa_bio04__within|+0.0156 [+0.00465, +0.0266]|
|bract_projection_p95|1303|111|chelsa_bio12__between|-0.00394 [-0.0142, +0.00632]|
|bract_projection_p95|1303|111|chelsa_bio12__within|-0.0023 [-0.00899, +0.00439]|
|bract_projection_p95|1303|111|chelsa_bio15__between|+0.00136 [-0.01, +0.0127]|
|bract_projection_p95|1303|111|chelsa_bio15__within|+0.0019 [-0.00417, +0.00796]|
|bract_projection_p95|1303|111|chelsa_rsds_mean__between|-0.0102 [-0.0261, +0.00558]|
|bract_projection_p95|1303|111|chelsa_rsds_mean__within|+0.00251 [-0.0127, +0.0177]|
|bract_projection_p95|1303|111|chelsa_vpd_mean__between|-0.0124 [-0.0257, +0.00079]|
|bract_projection_p95|1303|111|chelsa_vpd_mean__within|+0.00324 [-0.00915, +0.0156]|
|bract_projection_p95|1303|111|chelsa_sfcwind_mean__between|-0.00407 [-0.021, +0.0128]|
|bract_projection_p95|1303|111|chelsa_sfcwind_mean__within|+0.00731 [-0.000742, +0.0154]|
|bract_projection_p95|1303|111|chelsa_gsp__between|-0.00414 [-0.0113, +0.00304]|
|bract_projection_p95|1303|111|chelsa_gsp__within|-0.00199 [-0.00314, -0.000841]|
|bract_projection_p95|1303|111|chelsa_npp__between|+0.00652 [-0.00574, +0.0188]|
|bract_projection_p95|1303|111|chelsa_npp__within|+0.0189 [+0.00802, +0.0298]|
|bract_spread_fraction|1303|111|chelsa_bio04__between|+0.00089 [-0.0109, +0.0126]|
|bract_spread_fraction|1303|111|chelsa_bio04__within|+0.00518 [-0.00374, +0.0141]|
|bract_spread_fraction|1303|111|chelsa_bio12__between|-0.00953 [-0.0193, +0.000218]|
|bract_spread_fraction|1303|111|chelsa_bio12__within|+0.00201 [-0.0087, +0.0127]|
|bract_spread_fraction|1303|111|chelsa_bio15__between|+0.00492 [-0.00901, +0.0188]|
|bract_spread_fraction|1303|111|chelsa_bio15__within|-0.000987 [-0.00773, +0.00576]|
|bract_spread_fraction|1303|111|chelsa_rsds_mean__between|+1.19e-05 [-0.0144, +0.0145]|
|bract_spread_fraction|1303|111|chelsa_rsds_mean__within|-0.00218 [-0.0154, +0.011]|
|bract_spread_fraction|1303|111|chelsa_vpd_mean__between|-0.0209 [-0.0366, -0.00529]|
|bract_spread_fraction|1303|111|chelsa_vpd_mean__within|-0.00316 [-0.0221, +0.0158]|
|bract_spread_fraction|1303|111|chelsa_sfcwind_mean__between|-0.00123 [-0.0143, +0.0118]|
|bract_spread_fraction|1303|111|chelsa_sfcwind_mean__within|+0.00174 [-0.00428, +0.00776]|
|bract_spread_fraction|1303|111|chelsa_gsp__between|+5.43e-05 [-0.00833, +0.00844]|
|bract_spread_fraction|1303|111|chelsa_gsp__within|-0.00265 [-0.00552, +0.000214]|
|bract_spread_fraction|1303|111|chelsa_npp__between|+0.0154 [+0.00176, +0.0291]|
|bract_spread_fraction|1303|111|chelsa_npp__within|+0.00122 [-0.00879, +0.0112]|

## 環境の重なりと解釈の限界

相関・VIF・欠測数はJSONに全成分分を保存。高いVIFを理由に結果を見て追加で変数を削除しない。選定は母集団の環境値で行ったため、種間・種内や測定対象群では相関が残ることがある。相関が強い場合、個別係数は不安定であり、独立した環境作用を識別できたとは限らない。

種内の関連は可塑性の証明ではなく、種間の関連も適応の証明ではない。空間・系統・未測定環境を調整した結論ではない。写真の色は色素量ではなく、向きは画像基準角、総苞片は輪郭突出の代理指標である。

再現情報：Actions run 39416125522、commit b8e18fdb39c73326267649fd29dc4000912c5c9d。

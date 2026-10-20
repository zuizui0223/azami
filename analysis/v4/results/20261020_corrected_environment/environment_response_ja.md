# 頭花形質の環境勾配応答

問い：他の環境条件を同時に考慮したとき、どの環境勾配に沿って頭花形質が変わるか。

6群の全9候補から環境データのみでvifstep選定し、採択セットを固定して同時投入：温度季節性、降水、降水季節性、日射、VPD、風速、生育期降水量、NPP。観測値を平均せず、環境を種平均と種内偏差に分け、成分ごとに種ランダム切片を推定した。色全体は明度・彩度・円周色相の共同検定。各成分の利用可能な観測を使う。

過去の結果を閲覧済みの改訂観測解析であり、事前登録した独立検証ではない。64共同検定でHolm補正。成分の95%区間は点別である。

主解析は色全体・向き・頭花輪郭4項目。総苞片はSIに置く。補正対象は今回の全モジュール・環境・種内／種間の共同検定。

|形質|環境|種間P（Holm）|種内P（Holm）|
|---|---|---:|---:|
|色全体|温度季節性|1|0.009996|
|色全体|降水|0.1482|1|
|色全体|降水季節性|0.0498|1|
|色全体|日射|4.385e-06|0.00153|
|色全体|VPD|0.0003518|0.04544|
|色全体|風速|1|4.633e-11|
|色全体|生育期降水量|1|1|
|色全体|NPP|0.04343|1|
|向き|温度季節性|1|1|
|向き|降水|1|1|
|向き|降水季節性|1|1|
|向き|日射|1|0.3388|
|向き|VPD|1|0.1484|
|向き|風速|0.4653|1|
|向き|生育期降水量|1|1|
|向き|NPP|1|1|
|頭花輪郭4項目|温度季節性|1|4.889e-16|
|頭花輪郭4項目|降水|0.4287|1|
|頭花輪郭4項目|降水季節性|1|1|
|頭花輪郭4項目|日射|6.921e-11|0.04577|
|頭花輪郭4項目|VPD|9.963e-06|0.6011|
|頭花輪郭4項目|風速|1|1|
|頭花輪郭4項目|生育期降水量|0.01365|1|
|頭花輪郭4項目|NPP|0.00233|1|

## 各成分の変化

全自生域母集団の環境1標準偏差当たりの変化。単独の成分検定で色全体の関連を判定しない。色相sin/cosは共同で読む。

|成分|観測数|種数|環境・尺度|変化［95%区間］|
|---|---:|---:|---|---|
|corolla_lab_lightness|185564|244|chelsa_bio04__between|+0.671 [-0.344, +1.69]|
|corolla_lab_lightness|185564|244|chelsa_bio04__within|-0.244 [-0.769, +0.28]|
|corolla_lab_lightness|185564|244|chelsa_bio12__between|-1.49 [-2.5, -0.482]|
|corolla_lab_lightness|185564|244|chelsa_bio12__within|+0.0448 [-0.748, +0.838]|
|corolla_lab_lightness|185564|244|chelsa_bio15__between|-1.75 [-2.55, -0.949]|
|corolla_lab_lightness|185564|244|chelsa_bio15__within|-0.375 [-0.841, +0.0902]|
|corolla_lab_lightness|185564|244|chelsa_rsds_mean__between|+1.98 [+0.473, +3.49]|
|corolla_lab_lightness|185564|244|chelsa_rsds_mean__within|-0.719 [-1.91, +0.47]|
|corolla_lab_lightness|185564|244|chelsa_vpd_mean__between|-0.0342 [-1.12, +1.05]|
|corolla_lab_lightness|185564|244|chelsa_vpd_mean__within|+1.81 [+0.321, +3.29]|
|corolla_lab_lightness|185564|244|chelsa_sfcwind_mean__between|+0.19 [-0.542, +0.923]|
|corolla_lab_lightness|185564|244|chelsa_sfcwind_mean__within|+0.0845 [-0.372, +0.541]|
|corolla_lab_lightness|185564|244|chelsa_gsp__between|+1.47 [+0.18, +2.75]|
|corolla_lab_lightness|185564|244|chelsa_gsp__within|-0.225 [-0.999, +0.55]|
|corolla_lab_lightness|185564|244|chelsa_npp__between|-0.271 [-1.45, +0.907]|
|corolla_lab_lightness|185564|244|chelsa_npp__within|-0.187 [-1.2, +0.829]|
|corolla_lab_chroma|185564|244|chelsa_bio04__between|-0.602 [-1.76, +0.552]|
|corolla_lab_chroma|185564|244|chelsa_bio04__within|-0.727 [-1.56, +0.105]|
|corolla_lab_chroma|185564|244|chelsa_bio12__between|-0.362 [-1.25, +0.528]|
|corolla_lab_chroma|185564|244|chelsa_bio12__within|+0.503 [-0.0772, +1.08]|
|corolla_lab_chroma|185564|244|chelsa_bio15__between|+1.14 [+0.211, +2.08]|
|corolla_lab_chroma|185564|244|chelsa_bio15__within|+0.744 [+0.252, +1.24]|
|corolla_lab_chroma|185564|244|chelsa_rsds_mean__between|-4.48 [-6.13, -2.82]|
|corolla_lab_chroma|185564|244|chelsa_rsds_mean__within|-0.396 [-1.28, +0.485]|
|corolla_lab_chroma|185564|244|chelsa_vpd_mean__between|+0.349 [-0.869, +1.57]|
|corolla_lab_chroma|185564|244|chelsa_vpd_mean__within|-1.14 [-1.7, -0.579]|
|corolla_lab_chroma|185564|244|chelsa_sfcwind_mean__between|-0.095 [-0.787, +0.597]|
|corolla_lab_chroma|185564|244|chelsa_sfcwind_mean__within|+0.179 [-0.154, +0.513]|
|corolla_lab_chroma|185564|244|chelsa_gsp__between|-1.24 [-2.26, -0.229]|
|corolla_lab_chroma|185564|244|chelsa_gsp__within|-0.609 [-1.14, -0.0811]|
|corolla_lab_chroma|185564|244|chelsa_npp__between|+2.33 [+1.09, +3.56]|
|corolla_lab_chroma|185564|244|chelsa_npp__within|+0.977 [-0.313, +2.27]|
|corolla_hue_sin|185515|244|chelsa_bio04__between|-0.0371 [-0.106, +0.0322]|
|corolla_hue_sin|185515|244|chelsa_bio04__within|+0.0846 [+0.0404, +0.129]|
|corolla_hue_sin|185515|244|chelsa_bio12__between|-0.0598 [-0.112, -0.00737]|
|corolla_hue_sin|185515|244|chelsa_bio12__within|-0.0202 [-0.0642, +0.0238]|
|corolla_hue_sin|185515|244|chelsa_bio15__between|-0.064 [-0.118, -0.00988]|
|corolla_hue_sin|185515|244|chelsa_bio15__within|-0.0341 [-0.123, +0.0545]|
|corolla_hue_sin|185515|244|chelsa_rsds_mean__between|+0.181 [+0.083, +0.28]|
|corolla_hue_sin|185515|244|chelsa_rsds_mean__within|-0.0433 [-0.103, +0.0165]|
|corolla_hue_sin|185515|244|chelsa_vpd_mean__between|-0.155 [-0.227, -0.0829]|
|corolla_hue_sin|185515|244|chelsa_vpd_mean__within|+0.0391 [-0.0171, +0.0952]|
|corolla_hue_sin|185515|244|chelsa_sfcwind_mean__between|+0.00852 [-0.0356, +0.0526]|
|corolla_hue_sin|185515|244|chelsa_sfcwind_mean__within|+0.0119 [-0.00729, +0.031]|
|corolla_hue_sin|185515|244|chelsa_gsp__between|+0.0638 [-0.00514, +0.133]|
|corolla_hue_sin|185515|244|chelsa_gsp__within|+0.000916 [-0.0612, +0.063]|
|corolla_hue_sin|185515|244|chelsa_npp__between|-0.112 [-0.193, -0.0316]|
|corolla_hue_sin|185515|244|chelsa_npp__within|+0.0267 [-0.0155, +0.069]|
|corolla_hue_cos|185515|244|chelsa_bio04__between|-0.0185 [-0.0389, +0.00188]|
|corolla_hue_cos|185515|244|chelsa_bio04__within|+0.0093 [-0.00314, +0.0217]|
|corolla_hue_cos|185515|244|chelsa_bio12__between|+0.0104 [-0.0159, +0.0367]|
|corolla_hue_cos|185515|244|chelsa_bio12__within|-0.00765 [-0.0235, +0.00823]|
|corolla_hue_cos|185515|244|chelsa_bio15__between|+0.019 [+0.00122, +0.0368]|
|corolla_hue_cos|185515|244|chelsa_bio15__within|-0.00191 [-0.0161, +0.0123]|
|corolla_hue_cos|185515|244|chelsa_rsds_mean__between|+0.0294 [+0.000398, +0.0583]|
|corolla_hue_cos|185515|244|chelsa_rsds_mean__within|+0.0363 [+0.0158, +0.0567]|
|corolla_hue_cos|185515|244|chelsa_vpd_mean__between|-0.000126 [-0.0252, +0.025]|
|corolla_hue_cos|185515|244|chelsa_vpd_mean__within|-0.0267 [-0.0504, -0.00291]|
|corolla_hue_cos|185515|244|chelsa_sfcwind_mean__between|-0.00791 [-0.0259, +0.0101]|
|corolla_hue_cos|185515|244|chelsa_sfcwind_mean__within|+0.00899 [+0.00197, +0.016]|
|corolla_hue_cos|185515|244|chelsa_gsp__between|-0.0136 [-0.048, +0.0209]|
|corolla_hue_cos|185515|244|chelsa_gsp__within|+0.0152 [-0.00787, +0.0382]|
|corolla_hue_cos|185515|244|chelsa_npp__between|-0.0132 [-0.0402, +0.0137]|
|corolla_hue_cos|185515|244|chelsa_npp__within|-0.00123 [-0.0223, +0.0199]|
|orientation_image_vertical_angle|169263|244|chelsa_bio04__between|-0.442 [-3.23, +2.34]|
|orientation_image_vertical_angle|169263|244|chelsa_bio04__within|+0.0397 [-0.413, +0.493]|
|orientation_image_vertical_angle|169263|244|chelsa_bio12__between|+0.913 [-2.05, +3.88]|
|orientation_image_vertical_angle|169263|244|chelsa_bio12__within|+0.071 [-1.37, +1.51]|
|orientation_image_vertical_angle|169263|244|chelsa_bio15__between|+2 [-0.287, +4.28]|
|orientation_image_vertical_angle|169263|244|chelsa_bio15__within|-0.492 [-1.21, +0.226]|
|orientation_image_vertical_angle|169263|244|chelsa_rsds_mean__between|-0.27 [-3.62, +3.08]|
|orientation_image_vertical_angle|169263|244|chelsa_rsds_mean__within|-1.78 [-3.07, -0.49]|
|orientation_image_vertical_angle|169263|244|chelsa_vpd_mean__between|-0.279 [-2.67, +2.12]|
|orientation_image_vertical_angle|169263|244|chelsa_vpd_mean__within|+2.83 [+0.97, +4.7]|
|orientation_image_vertical_angle|169263|244|chelsa_sfcwind_mean__between|+2.33 [+0.554, +4.11]|
|orientation_image_vertical_angle|169263|244|chelsa_sfcwind_mean__within|-0.122 [-0.575, +0.33]|
|orientation_image_vertical_angle|169263|244|chelsa_gsp__between|+2.03 [-1.2, +5.26]|
|orientation_image_vertical_angle|169263|244|chelsa_gsp__within|-0.64 [-1.97, +0.689]|
|orientation_image_vertical_angle|169263|244|chelsa_npp__between|-1.48 [-4.2, +1.24]|
|orientation_image_vertical_angle|169263|244|chelsa_npp__within|+0.149 [-0.895, +1.19]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_bio04__between|+0.00276 [-0.0045, +0.01]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_bio04__within|+0.00552 [+0.0044, +0.00665]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_bio12__between|-0.000894 [-0.0096, +0.00781]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_bio12__within|-0.000134 [-0.00298, +0.00271]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_bio15__between|-0.00434 [-0.0117, +0.00304]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_bio15__within|+0.00242 [-0.000363, +0.0052]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_rsds_mean__between|+0.00419 [-0.00698, +0.0154]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_rsds_mean__within|-0.00391 [-0.00897, +0.00114]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_vpd_mean__between|-0.00253 [-0.00788, +0.00282]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_vpd_mean__within|+0.000435 [-0.00237, +0.00324]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_sfcwind_mean__between|-0.00319 [-0.00768, +0.0013]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_sfcwind_mean__within|-0.000152 [-0.00152, +0.00122]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_gsp__between|-0.00296 [-0.0136, +0.00768]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_gsp__within|+0.00145 [-0.00146, +0.00435]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_npp__between|+0.00779 [-0.00457, +0.0201]|
|capitulum_outline_aspect_ratio|175223|243|chelsa_npp__within|+0.00118 [-0.000809, +0.00317]|
|capitulum_outline_circularity|175223|243|chelsa_bio04__between|-0.00802 [-0.0142, -0.00183]|
|capitulum_outline_circularity|175223|243|chelsa_bio04__within|+0.0027 [+0.00145, +0.00396]|
|capitulum_outline_circularity|175223|243|chelsa_bio12__between|-0.00941 [-0.0154, -0.00346]|
|capitulum_outline_circularity|175223|243|chelsa_bio12__within|+0.00135 [-0.00215, +0.00485]|
|capitulum_outline_circularity|175223|243|chelsa_bio15__between|-0.000121 [-0.00483, +0.00459]|
|capitulum_outline_circularity|175223|243|chelsa_bio15__within|+0.00324 [+6.07e-05, +0.00642]|
|capitulum_outline_circularity|175223|243|chelsa_rsds_mean__between|-0.0137 [-0.0215, -0.00577]|
|capitulum_outline_circularity|175223|243|chelsa_rsds_mean__within|-0.00816 [-0.0148, -0.00155]|
|capitulum_outline_circularity|175223|243|chelsa_vpd_mean__between|+0.00171 [-0.00479, +0.0082]|
|capitulum_outline_circularity|175223|243|chelsa_vpd_mean__within|+0.00286 [-0.00094, +0.00666]|
|capitulum_outline_circularity|175223|243|chelsa_sfcwind_mean__between|+0.00369 [-0.000343, +0.00771]|
|capitulum_outline_circularity|175223|243|chelsa_sfcwind_mean__within|+0.000344 [-0.00177, +0.00246]|
|capitulum_outline_circularity|175223|243|chelsa_gsp__between|+0.00847 [+0.000407, +0.0165]|
|capitulum_outline_circularity|175223|243|chelsa_gsp__within|-0.000414 [-0.00407, +0.00324]|
|capitulum_outline_circularity|175223|243|chelsa_npp__between|-0.00581 [-0.0125, +0.000921]|
|capitulum_outline_circularity|175223|243|chelsa_npp__within|-0.000204 [-0.00233, +0.00192]|
|capitulum_outline_solidity|175223|243|chelsa_bio04__between|-0.00425 [-0.0079, -0.000606]|
|capitulum_outline_solidity|175223|243|chelsa_bio04__within|+0.000417 [-0.000114, +0.000948]|
|capitulum_outline_solidity|175223|243|chelsa_bio12__between|-0.00534 [-0.00869, -0.00198]|
|capitulum_outline_solidity|175223|243|chelsa_bio12__within|-6.31e-05 [-0.00145, +0.00132]|
|capitulum_outline_solidity|175223|243|chelsa_bio15__between|-0.000663 [-0.00368, +0.00235]|
|capitulum_outline_solidity|175223|243|chelsa_bio15__within|+0.0014 [-0.000462, +0.00327]|
|capitulum_outline_solidity|175223|243|chelsa_rsds_mean__between|+0.000222 [-0.00449, +0.00493]|
|capitulum_outline_solidity|175223|243|chelsa_rsds_mean__within|-0.00293 [-0.00652, +0.000653]|
|capitulum_outline_solidity|175223|243|chelsa_vpd_mean__between|-0.00211 [-0.00601, +0.0018]|
|capitulum_outline_solidity|175223|243|chelsa_vpd_mean__within|+0.000513 [-0.000996, +0.00202]|
|capitulum_outline_solidity|175223|243|chelsa_sfcwind_mean__between|+0.00247 [+0.000148, +0.0048]|
|capitulum_outline_solidity|175223|243|chelsa_sfcwind_mean__within|+0.00027 [-0.00108, +0.00162]|
|capitulum_outline_solidity|175223|243|chelsa_gsp__between|+0.00579 [+0.00172, +0.00985]|
|capitulum_outline_solidity|175223|243|chelsa_gsp__within|+0.000562 [-0.000977, +0.0021]|
|capitulum_outline_solidity|175223|243|chelsa_npp__between|-0.00377 [-0.00757, +3.17e-05]|
|capitulum_outline_solidity|175223|243|chelsa_npp__within|-8.44e-05 [-0.00126, +0.00109]|
|capitulum_width_profile_cv|175223|243|chelsa_bio04__between|+0.00154 [-0.00273, +0.00581]|
|capitulum_width_profile_cv|175223|243|chelsa_bio04__within|-0.000807 [-0.00181, +0.000198]|
|capitulum_width_profile_cv|175223|243|chelsa_bio12__between|+0.0068 [+0.00269, +0.0109]|
|capitulum_width_profile_cv|175223|243|chelsa_bio12__within|+0.000672 [-0.00131, +0.00265]|
|capitulum_width_profile_cv|175223|243|chelsa_bio15__between|+0.00063 [-0.0024, +0.00366]|
|capitulum_width_profile_cv|175223|243|chelsa_bio15__within|-0.00148 [-0.00353, +0.000574]|
|capitulum_width_profile_cv|175223|243|chelsa_rsds_mean__between|-0.00589 [-0.0112, -0.000608]|
|capitulum_width_profile_cv|175223|243|chelsa_rsds_mean__within|+0.00381 [+0.000307, +0.0073]|
|capitulum_width_profile_cv|175223|243|chelsa_vpd_mean__between|+0.00628 [+0.00232, +0.0102]|
|capitulum_width_profile_cv|175223|243|chelsa_vpd_mean__within|-0.00245 [-0.00397, -0.000929]|
|capitulum_width_profile_cv|175223|243|chelsa_sfcwind_mean__between|-0.00127 [-0.00425, +0.00171]|
|capitulum_width_profile_cv|175223|243|chelsa_sfcwind_mean__within|+2.81e-05 [-0.00142, +0.00148]|
|capitulum_width_profile_cv|175223|243|chelsa_gsp__between|-0.0109 [-0.0158, -0.00597]|
|capitulum_width_profile_cv|175223|243|chelsa_gsp__within|-0.000629 [-0.00346, +0.0022]|
|capitulum_width_profile_cv|175223|243|chelsa_npp__between|+0.00908 [+0.00502, +0.0131]|
|capitulum_width_profile_cv|175223|243|chelsa_npp__within|+0.000179 [-0.00217, +0.00253]|

## 環境の重なりと解釈の限界

相関・VIF・欠測数はJSONに全成分分を保存。高いVIFを理由に結果を見て追加で変数を削除しない。選定は母集団の環境値で行ったため、種間・種内や測定対象群では相関が残ることがある。相関が強い場合、個別係数は不安定であり、独立した環境作用を識別できたとは限らない。

種内の関連は可塑性の証明ではなく、種間の関連も適応の証明ではない。種内と種間の係数差を直接検定していないため、片方だけの有意性を尺度間の差の証拠にはしない。空間・系統・未測定環境を調整した結論ではない。写真の色は色素量ではなく、向きは画像基準角、総苞片は輪郭突出の代理指標である。

再現情報：Actions run 39422378058、commit 93e8d874bc2c5c09ba5692e2c77bad723194ee5f。

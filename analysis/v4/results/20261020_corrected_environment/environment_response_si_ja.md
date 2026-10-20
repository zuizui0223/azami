# SI：総苞片の輪郭突出・広がり

既存測定1,303観測の補足解析。主解析より観測規模・対象被覆が限定されるためSIに置く。この判断は有意性による選別ではなく、測定範囲と指標の制約に基づく。主解析と合わせた64検定のHolm補正を使用する。

輪郭突出・広がりは個々の総苞片の開出角や棘による防御効果そのものではない。気候との関連が明確でなくても、食害者・訪花者に対する防御機能の否定にはならない。

|環境|種間P（Holm）|種内P（Holm）|
|---|---:|---:|
|温度季節性|1|0.3876|
|降水|1|1|
|降水季節性|1|1|
|日射|1|1|
|VPD|0.5661|1|
|風速|1|1|
|生育期降水量|1|1|
|NPP|1|1|

|成分|観測数|種数|環境・尺度|変化［点別95%区間］|
|---|---:|---:|---|---|
|bract_projection_p95|1298|111|chelsa_bio04__between|+0.00168 [-0.0109, +0.0142]|
|bract_projection_p95|1298|111|chelsa_bio04__within|+0.0155 [+0.00515, +0.0258]|
|bract_projection_p95|1298|111|chelsa_bio12__between|-0.00847 [-0.0241, +0.00719]|
|bract_projection_p95|1298|111|chelsa_bio12__within|-0.000481 [-0.0118, +0.0108]|
|bract_projection_p95|1298|111|chelsa_bio15__between|+0.000245 [-0.0116, +0.0121]|
|bract_projection_p95|1298|111|chelsa_bio15__within|+0.00217 [-0.00422, +0.00856]|
|bract_projection_p95|1298|111|chelsa_rsds_mean__between|-0.00907 [-0.0252, +0.0071]|
|bract_projection_p95|1298|111|chelsa_rsds_mean__within|+0.00183 [-0.0137, +0.0173]|
|bract_projection_p95|1298|111|chelsa_vpd_mean__between|-0.0137 [-0.0277, +0.000359]|
|bract_projection_p95|1298|111|chelsa_vpd_mean__within|+0.004 [-0.0089, +0.0169]|
|bract_projection_p95|1298|111|chelsa_sfcwind_mean__between|-0.00376 [-0.0206, +0.013]|
|bract_projection_p95|1298|111|chelsa_sfcwind_mean__within|+0.00722 [-0.000761, +0.0152]|
|bract_projection_p95|1298|111|chelsa_gsp__between|+0.00626 [-0.0142, +0.0267]|
|bract_projection_p95|1298|111|chelsa_gsp__within|-0.00235 [-0.0199, +0.0152]|
|bract_projection_p95|1298|111|chelsa_npp__between|+0.00472 [-0.00892, +0.0183]|
|bract_projection_p95|1298|111|chelsa_npp__within|+0.0195 [+0.00467, +0.0344]|
|bract_spread_fraction|1298|111|chelsa_bio04__between|+0.00222 [-0.01, +0.0144]|
|bract_spread_fraction|1298|111|chelsa_bio04__within|+0.00564 [-0.00329, +0.0146]|
|bract_spread_fraction|1298|111|chelsa_bio12__between|-0.0171 [-0.0362, +0.00201]|
|bract_spread_fraction|1298|111|chelsa_bio12__within|+0.00141 [-0.0134, +0.0162]|
|bract_spread_fraction|1298|111|chelsa_bio15__between|+0.0029 [-0.0111, +0.0169]|
|bract_spread_fraction|1298|111|chelsa_bio15__within|-0.00118 [-0.00812, +0.00575]|
|bract_spread_fraction|1298|111|chelsa_rsds_mean__between|+0.00211 [-0.012, +0.0163]|
|bract_spread_fraction|1298|111|chelsa_rsds_mean__within|-0.00303 [-0.0171, +0.0111]|
|bract_spread_fraction|1298|111|chelsa_vpd_mean__between|-0.0232 [-0.0386, -0.00774]|
|bract_spread_fraction|1298|111|chelsa_vpd_mean__within|-0.00267 [-0.0227, +0.0174]|
|bract_spread_fraction|1298|111|chelsa_sfcwind_mean__between|-0.000875 [-0.0141, +0.0123]|
|bract_spread_fraction|1298|111|chelsa_sfcwind_mean__within|+0.00133 [-0.00469, +0.00736]|
|bract_spread_fraction|1298|111|chelsa_gsp__between|+0.0148 [-0.0115, +0.0412]|
|bract_spread_fraction|1298|111|chelsa_gsp__within|+0.00379 [-0.0139, +0.0215]|
|bract_spread_fraction|1298|111|chelsa_npp__between|+0.00931 [-0.00655, +0.0252]|
|bract_spread_fraction|1298|111|chelsa_npp__within|-0.000337 [-0.0117, +0.0111]|

全係数・診断・来歴はenvironment_response.jsonに保存。再測定はこの結果の成立条件にしない。
推定元：Actions run 39422378058、commit 93e8d874bc2c5c09ba5692e2c77bad723194ee5f。

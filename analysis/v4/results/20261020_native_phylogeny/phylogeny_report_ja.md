# V4 系統配置感度検証

Actions run: 39435174409; code: `c1266ec807f83b60e111a5d9125e7b14abc29786`.

自生域の観測値を保持し、既存の5環境要因の種内・種間成分を同時に投入した。
系統に沿う種効果と独立な種効果を分けた観測レベルのGaussian REMLを使用。
比較対象は、同じ樹掲載種・同じ観測だけを使う独立種効果モデル。
52本は既存の代替配置樹であり、独立した52回の実証ではない。
区間は分散成分を推定値に固定したモデルベース近似95%区間であり、以前の種クラスター区間とは異なる。全樹の区間の外包は統合信頼区間ではない。
個別成分の区間は多重比較補正前。色全体・輪郭全体の合同検定、系統的形質間共分散、適応・種分化の検証とは区別する。

## 樹の網羅率

|応答|元の観測数|樹掲載観測数|元の種数|樹掲載種数|
|---|---:|---:|---:|---:|
|orientation_image_vertical_angle|170066|169912|244|194|
|corolla_lab_lightness|186439|186265|244|194|
|corolla_lab_chroma|186439|186265|244|194|
|corolla_hue_sin|186390|186216|244|194|
|corolla_hue_cos|186390|186216|244|194|
|capitulum_outline_aspect_ratio|176048|175880|243|194|
|capitulum_outline_circularity|176048|175880|243|194|
|capitulum_outline_solidity|176048|175880|243|194|
|capitulum_width_profile_cv|176048|175880|243|194|

## 向きと降水

- chelsa_bio12__between: 降水100 mm当たり、独立種モデル 0.4880 [0.2220, 0.7539] 度。系統モデルの係数範囲 0.4722〜0.5078 度、区間の外包 [0.2051, 0.7748]。正の区間 52/52、負の区間 0/52。
- chelsa_bio12__within: 降水100 mm当たり、独立種モデル -0.0454 [-0.1131, 0.0223] 度。系統モデルの係数範囲 -0.0457〜-0.0454 度、区間の外包 [-0.1134, 0.0223]。正の区間 0/52、負の区間 0/52。

## 全成分

係数は表記した環境増分当たりの元の形質単位。

|応答・項|環境増分|系統係数範囲|区間の外包|正/負区間の樹数|
|---|---:|---:|---:|---:|
|capitulum_outline_aspect_ratio / chelsa_bio01__between|1|0.00046474〜0.00059714|[-0.00068045, 0.0017419]|0/0|
|capitulum_outline_aspect_ratio / chelsa_bio01__within|1|0.00025123〜0.00025207|[-9.4458e-05, 0.00059776]|0/0|
|capitulum_outline_aspect_ratio / chelsa_bio04__between|100|0.00020865〜0.00072773|[-0.0031446, 0.0040505]|0/0|
|capitulum_outline_aspect_ratio / chelsa_bio04__within|100|0.0024108〜0.0024123|[0.0018968, 0.0029263]|52/0|
|capitulum_outline_aspect_ratio / chelsa_bio12__between|100|-0.00027074〜-0.00015702|[-0.0010506, 0.00060494]|0/0|
|capitulum_outline_aspect_ratio / chelsa_bio12__within|100|0.00032447〜0.00032553|[9.5554e-05, 0.00055445]|52/0|
|capitulum_outline_aspect_ratio / chelsa_bio15__between|10|-0.0021218〜-0.0016971|[-0.0051507, 0.0013748]|0/0|
|capitulum_outline_aspect_ratio / chelsa_bio15__within|10|0.0011502〜0.0011525|[0.00052476, 0.001778]|52/0|
|capitulum_outline_aspect_ratio / chelsa_rsds_mean__between|1|-0.00015669〜0.00055397|[-0.0025352, 0.0028933]|0/0|
|capitulum_outline_aspect_ratio / chelsa_rsds_mean__within|1|-0.0011824〜-0.0011816|[-0.0017547, -0.00060921]|0/52|
|capitulum_outline_circularity / chelsa_bio01__between|1|-0.00010082〜0.00030239|[-0.001089, 0.00127]|0/0|
|capitulum_outline_circularity / chelsa_bio01__within|1|0.00020456〜0.00020596|[-6.3719e-05, 0.00047423]|0/0|
|capitulum_outline_circularity / chelsa_bio04__between|100|-0.0031349〜-0.0020777|[-0.0059585, 0.00073531]|0/2|
|capitulum_outline_circularity / chelsa_bio04__within|100|0.0012673〜0.0012689|[0.0008684, 0.0016678]|52/0|
|capitulum_outline_circularity / chelsa_bio12__between|100|-0.00072113〜-0.00048226|[-0.0013741, 0.00017251]|0/6|
|capitulum_outline_circularity / chelsa_bio12__within|100|0.000126〜0.00012753|[-5.1672e-05, 0.0003052]|0/0|
|capitulum_outline_circularity / chelsa_bio15__between|10|-0.00020091〜0.00075106|[-0.0027052, 0.0033186]|0/0|
|capitulum_outline_circularity / chelsa_bio15__within|10|0.0015816〜0.0015843|[0.0010962, 0.0020697]|52/0|
|capitulum_outline_circularity / chelsa_rsds_mean__between|1|-0.0038416〜-0.0030396|[-0.0059222, -0.00098415]|0/52|
|capitulum_outline_circularity / chelsa_rsds_mean__within|1|-0.0020535〜-0.0020516|[-0.0024977, -0.0016074]|0/52|
|capitulum_outline_solidity / chelsa_bio01__between|1|-0.00051054〜-0.00034995|[-0.0010703, 0.00020496]|0/0|
|capitulum_outline_solidity / chelsa_bio01__within|1|-5.1494e-05〜-5.0743e-05|[-0.00021083, 0.00010859]|0/0|
|capitulum_outline_solidity / chelsa_bio04__between|100|-0.0017403〜-0.0013733|[-0.003364, 0.00025034]|0/4|
|capitulum_outline_solidity / chelsa_bio04__within|100|2.9346e-05〜3.0163e-05|[-0.00020757, 0.00026707]|0/0|
|capitulum_outline_solidity / chelsa_bio12__between|100|-0.00017694〜-0.00010126|[-0.00055412, 0.00027652]|0/0|
|capitulum_outline_solidity / chelsa_bio12__within|100|2.2128e-05〜2.2888e-05|[-8.3384e-05, 0.0001284]|0/0|
|capitulum_outline_solidity / chelsa_bio15__between|10|-0.0002781〜0.0001174|[-0.0017105, 0.0015969]|0/0|
|capitulum_outline_solidity / chelsa_bio15__within|10|0.00066386〜0.00066542|[0.00037557, 0.00095371]|52/0|
|capitulum_outline_solidity / chelsa_rsds_mean__between|1|-0.00014352〜0.00025002|[-0.001317, 0.0014076]|0/0|
|capitulum_outline_solidity / chelsa_rsds_mean__within|1|-0.00075191〜-0.00075112|[-0.0010157, -0.00048731]|0/52|
|capitulum_width_profile_cv / chelsa_bio01__between|1|0.0012061〜0.0013855|[0.00057892, 0.0020155]|52/0|
|capitulum_width_profile_cv / chelsa_bio01__within|1|-0.00031591〜-0.00031526|[-0.00051766, -0.00011351]|0/52|
|capitulum_width_profile_cv / chelsa_bio04__between|100|0.00056731〜0.00094808|[-0.0012806, 0.0028133]|0/0|
|capitulum_width_profile_cv / chelsa_bio04__within|100|-0.00052801〜-0.00052701|[-0.00082799, -0.00022704]|0/52|
|capitulum_width_profile_cv / chelsa_bio12__between|100|-0.00021632〜-9.5983e-05|[-0.00065046, 0.00034257]|0/0|
|capitulum_width_profile_cv / chelsa_bio12__within|100|7.0971e-05〜7.1431e-05|[-6.2621e-05, 0.00020502]|0/0|
|capitulum_width_profile_cv / chelsa_bio15__between|10|-0.00075936〜-0.00045039|[-0.0024367, 0.0012415]|0/0|
|capitulum_width_profile_cv / chelsa_bio15__within|10|-0.00077754〜-0.00077561|[-0.0011426, -0.00041058]|0/52|
|capitulum_width_profile_cv / chelsa_rsds_mean__between|1|-0.0016854〜-0.0013111|[-0.0030177, -7.1771e-06]|0/52|
|capitulum_width_profile_cv / chelsa_rsds_mean__within|1|0.00094067〜0.00094121|[0.00060663, 0.0012752]|52/0|
|corolla_hue_cos / chelsa_bio01__between|1|-0.002256〜-0.0012169|[-0.0057795, 0.0023605]|0/0|
|corolla_hue_cos / chelsa_bio01__within|1|-0.0050831〜-0.0050788|[-0.0059, -0.004262]|0/52|
|corolla_hue_cos / chelsa_bio04__between|100|-0.01099〜-0.0078508|[-0.021327, 0.0023999]|0/3|
|corolla_hue_cos / chelsa_bio04__within|100|-0.0021103〜-0.0021052|[-0.003322, -0.00089349]|0/52|
|corolla_hue_cos / chelsa_bio12__between|100|-0.0015663〜-0.00062846|[-0.0038865, 0.0017229]|0/0|
|corolla_hue_cos / chelsa_bio12__within|100|0.00030176〜0.0003056|[-0.00023978, 0.00084716]|0/0|
|corolla_hue_cos / chelsa_bio15__between|10|0.0069677〜0.0094605|[-0.0023314, 0.018798]|2/0|
|corolla_hue_cos / chelsa_bio15__within|10|-0.0017349〜-0.00173|[-0.0032101, -0.00025481]|0/52|
|corolla_hue_cos / chelsa_rsds_mean__between|1|0.0067095〜0.00903|[-0.00074246, 0.016465]|47/0|
|corolla_hue_cos / chelsa_rsds_mean__within|1|0.0082046〜0.0082102|[0.0068511, 0.0095638]|52/0|
|corolla_hue_sin / chelsa_bio01__between|1|-0.036365〜-0.034523|[-0.047983, -0.022998]|0/52|
|corolla_hue_sin / chelsa_bio01__within|1|0.012846〜0.012848|[0.011428, 0.014265]|52/0|
|corolla_hue_sin / chelsa_bio04__between|100|-0.028911〜-0.02278|[-0.060805, 0.0093527]|0/0|
|corolla_hue_sin / chelsa_bio04__within|100|0.043465〜0.043468|[0.041362, 0.045571]|52/0|
|corolla_hue_sin / chelsa_bio12__between|100|-0.0049148〜-0.0036079|[-0.011957, 0.003448]|0/0|
|corolla_hue_sin / chelsa_bio12__within|100|-0.00077722〜-0.000775|[-0.0017174, 0.00016516]|0/0|
|corolla_hue_sin / chelsa_bio15__between|10|-0.024757〜-0.020445|[-0.054254, 0.0092747]|0/0|
|corolla_hue_sin / chelsa_bio15__within|10|-0.015274〜-0.01527|[-0.017835, -0.012709]|0/52|
|corolla_hue_sin / chelsa_rsds_mean__between|1|0.046232〜0.051181|[0.021669, 0.075889]|52/0|
|corolla_hue_sin / chelsa_rsds_mean__within|1|-0.016749〜-0.016745|[-0.019098, -0.014396]|0/52|
|corolla_lab_chroma / chelsa_bio01__between|1|0.22919〜0.2491|[0.045309, 0.43404]|52/0|
|corolla_lab_chroma / chelsa_bio01__within|1|-0.17005〜-0.17|[-0.19915, -0.1409]|0/52|
|corolla_lab_chroma / chelsa_bio04__between|100|-0.3348〜-0.2687|[-0.84705, 0.25064]|0/0|
|corolla_lab_chroma / chelsa_bio04__within|100|-0.52916〜-0.5291|[-0.57233, -0.48594]|0/52|
|corolla_lab_chroma / chelsa_bio12__between|100|-0.16465〜-0.14964|[-0.28033, -0.034024]|0/52|
|corolla_lab_chroma / chelsa_bio12__within|100|0.068558〜0.068584|[0.049266, 0.087876]|52/0|
|corolla_lab_chroma / chelsa_bio15__between|10|0.44508〜0.4993|[-0.029761, 0.97588]|2/0|
|corolla_lab_chroma / chelsa_bio15__within|10|0.28208〜0.28216|[0.22953, 0.33471]|52/0|
|corolla_lab_chroma / chelsa_rsds_mean__between|1|-1.3918〜-1.3143|[-1.7777, -0.92362]|0/52|
|corolla_lab_chroma / chelsa_rsds_mean__within|1|-0.16839〜-0.16833|[-0.21661, -0.12012]|0/52|
|corolla_lab_lightness / chelsa_bio01__between|1|-0.028324〜0.016427|[-0.21028, 0.19777]|0/0|
|corolla_lab_lightness / chelsa_bio01__within|1|0.32178〜0.32188|[0.29906, 0.3446]|52/0|
|corolla_lab_lightness / chelsa_bio04__between|100|0.29492〜0.40552|[-0.2092, 0.91161]|0/0|
|corolla_lab_lightness / chelsa_bio04__within|100|0.17467〜0.17476|[0.14097, 0.20846]|52/0|
|corolla_lab_lightness / chelsa_bio12__between|100|-0.061676〜-0.032792|[-0.17332, 0.080832]|0/0|
|corolla_lab_lightness / chelsa_bio12__within|100|-0.024711〜-0.024637|[-0.039775, -0.0095738]|0/52|
|corolla_lab_lightness / chelsa_bio15__between|10|-0.68671〜-0.62151|[-1.1533, -0.15705]|0/52|
|corolla_lab_lightness / chelsa_bio15__within|10|-0.11848〜-0.1184|[-0.15951, -0.077369]|0/52|
|corolla_lab_lightness / chelsa_rsds_mean__between|1|0.52672〜0.65863|[0.13741, 1.0464]|52/0|
|corolla_lab_lightness / chelsa_rsds_mean__within|1|-0.17854〜-0.17847|[-0.21619, -0.14082]|0/52|
|orientation_image_vertical_angle / chelsa_bio01__between|1|-0.20355〜-0.15784|[-0.60742, 0.24419]|0/0|
|orientation_image_vertical_angle / chelsa_bio01__within|1|0.51658〜0.51698|[0.41415, 0.61941]|52/0|
|orientation_image_vertical_angle / chelsa_bio04__between|100|0.017119〜0.27153|[-1.1446, 1.4021]|0/0|
|orientation_image_vertical_angle / chelsa_bio04__within|100|0.46937〜0.46979|[0.31756, 0.62159]|52/0|
|orientation_image_vertical_angle / chelsa_bio12__between|100|0.47222〜0.50784|[0.20507, 0.77483]|52/0|
|orientation_image_vertical_angle / chelsa_bio12__within|100|-0.045665〜-0.04541|[-0.11336, 0.022289]|0/0|
|orientation_image_vertical_angle / chelsa_bio15__between|10|0.76792〜0.9174|[-0.29014, 1.9809]|0/0|
|orientation_image_vertical_angle / chelsa_bio15__within|10|-0.15337〜-0.15265|[-0.33934, 0.033324]|0/0|
|orientation_image_vertical_angle / chelsa_rsds_mean__between|1|0.048171〜0.27016|[-0.7919, 1.0925]|0/0|
|orientation_image_vertical_angle / chelsa_rsds_mean__within|1|-0.43256〜-0.43228|[-0.60113, -0.26371]|0/52|

失敗記録: []

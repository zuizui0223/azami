# GEB V4：全球標高と *Cirsium arvense* 花冠画像明度の解釈 — 改稿提案

**Status: REVIEW ONLY.** 2026-11-07。元の凍結原稿 `GEB_RESEARCH_ARTICLE.md`、S1 80 Holm、S2 128 Holm、H3c帰無には一切変更しない。以下は方法・結果・考察への差し込み用案。新しい実測効果、因果解釈、有意性を主張しない。

## 認証済みの根拠：モデルが実際に推定したもの

- 全球自生分類群S1：撮影場所のWorldClim v2.1地表標高はCHELSA v2.1 BIO18（最暖四半期月平均降水量）・BIO05（最暖月最高気温）・RSDS（地表への下向き短波放射）と同時に調整。元4変数の最大VIF **4.9559**。気候×画像頭花向きの元のHolm陽性はBIO18の種間係数（β=+0.045447, Holm P=.00920635）。標高自体の独立したS1効果が確認されたという意味ではない。
- *C. arvense* S2：**移入域**花冠画像CIELAB L* × 標高：標準化偏回帰係数β=**−0.173731874**, pointwise 95%CI **−0.250914 to −0.096550**, original 128-test family Holm P=**0.00131233**, **26,171**適格写真。**在来域**β=**+0.056266732**, 95%CI **−0.018662 to +0.131195**, Holm P=1, **22,936**写真。回帰は両域でそれぞれ標準化（環境変数・明度とも層内のSD基準）、時期調整、観測者×5°系セルの写真重複除去、セル等重み・セルクラスタ誤差。
- 移入域10個の20°緯度×経度空間ブロックを399回再標本化した標高偏回帰の経験95%区間 **[−0.22949, −0.05578]**；最大寄与の5ブロックを個別に抜いたとき符号維持5/5。在来域24ブロックの区間 **[−0.05208,+0.12724]**。これは空間サンプル依存性の感度診断であり独立地理反復・空間交絡の除去を意味しない。
- **S2の標高×移入域交互作用は元凍結論文では未実施。** 最近のX-only気候5変数での交互作用80件には標高を含めておらず、これを標高の範囲差検定と誤記してはいけない。既存の移入域有意vs在来域非有意の差は勾配差のP値ではない。

## 標高という変数の位置づけ

WorldClim地表標高は**地理的位置の物理的属性**であり、気候や標高に伴う形質への個別の生理的刺激ではない。温度・気圧・大気透過に一定の物理的関係がある一方、降水・雲量・日射・積雪・表層地質・土壌・植生・管理・送粉者相は地理域ごとに別々の勾配を取り得る（Körner 2007）。同じ標高1,000 mの観測値が高緯度・熱帯・温帯島嶼で同じ温度・季節長・照明環境を意味しない。

BIO18/BIO05/RSDSを同時に調整した **conditional elevation coefficient** は、これらの広域気候平年値で線形的に説明されない標高方向の花冠画像明度との*残余関連*を要約する。標高が気候の形成原因でもある以上、同時調整後の係数は「標高が与える総効果」ではない。残余には気候平年値の解像度より細かな環境、短期の開花期温度・積雪・UV、土壌や生息地、地域／集団史、撮影条件、未測定の空間構造が混在し得る。「気候から独立した標高の適応効果」「高標高により色素増加が生じた」とは記述しない。

CIELAB L* は画像領域での測定であって、分光測定の反射率、色素含量、頭花内部温度、UV照射量ではない。RSDSは気候平年の総短波放射であり、撮影時の光やUVだけを表すものではない。暗く写る関係は（適切なQC後でも）写真の露出・日照・標本間形質・撮影者の偏りなどで一部説明され得る。

## 論文に直接使用できる改稿案（English, proposed）

### Methods 2.2 (replace/add after environmental source sentence)

> We extracted elevation (metres above sea level) from WorldClim v2.1 separately from the CHELSA v2.1 climatic normals. Across the global study extent, elevation represents a location-specific topographic gradient rather than a direct measurement of reproductive-season temperature, ultraviolet exposure or a common climatic treatment. Its regression coefficient is interpreted conditionally on warm-quarter precipitation (BIO18), maximum warm-month temperature (BIO05) and downwelling shortwave radiation (RSDS), without identifying an independent causal effect of elevation. The original S2 standardisation was performed within each species-by-range stratum; accordingly, standardised slopes cannot be interpreted as effects per fixed metre or compared between strata as if they shared identical environmental variation.

### Results 3.2 (replace only the altitude sentence, preserving all original numbers)

> In introduced *Cirsium arvense*, corolla-image lightness declined with elevation after adjustment for BIO18, BIO05 and RSDS (within-stratum standardised β = −0.173732, pointwise 95% CI −0.250914 to −0.096550; Holm-adjusted P = 0.00131233; 26,171 eligible photographs). The negative association persisted under 20° spatial-block resampling (399 draws; percentile interval −0.22949 to −0.05578), although only ten such blocks contributed to this introduced-range analysis. In the native range the estimated slope was +0.056267 (95% CI −0.018662 to +0.131195; 22,936 photographs), but the separate regressions and stratum-specific standardisation do not establish a native–introduced difference in slopes.

### Discussion: replace the current elevation-temperature paragraph with the following two paragraphs

> The negative association between elevation and photographic corolla lightness in introduced *C. arvense* is compatible with a thermal explanation, but is not specific to it. Elevation combines physically constrained gradients, including atmospheric pressure and generally declining air temperature, with locally contingent variation in precipitation, cloud cover, radiation, snow, soils and vegetation (Körner, 2007). Moreover, equivalent elevations across latitudes and regions do not imply equivalent climates. Because the observed coefficient was conditional on long-term BIO18, BIO05 and RSDS, it summarises residual environmental and geographic correspondence rather than the total effect of climbing an elevational gradient or a direct effect of elevation itself. Its persistence under spatial-block resampling makes the observation less dependent on the largest sampled geographic blocks, but does not exclude residual spatial confounding, population composition or photography-related variation.

> Darker reproductive structures can absorb more solar energy under some conditions: experimental work in *Plantago lanceolata* links developmental temperature, spike reflectance and reproductive-tissue heating (Lacey & Herr, 2005), while elevational variation in floral pigmentation in *Silene vulgaris* covaries with microclimatic conditions (Berardi et al., 2016). These studies make a temperature-mediated interpretation plausible, not demonstrated for *Cirsium*. Our CIELAB lightness values are image-derived and do not measure spectral reflectance, ultraviolet absorption, pigments or capitulum temperature. The native-range estimate was weakly positive and the introduced-range estimate negative, but neither separate significance decisions nor slopes standardised within distinct sampling distributions test a range-by-elevation interaction. Testing whether the association reflects local thermal conditions, other correlates of elevation, or differences in sampled populations would require comparable elevational and climatic support across regions, a shared-scale range interaction and direct floral reflectance and temperature measurements. We therefore interpret elevation as a geographically contingent environmental proxy associated with variation in photographed corolla lightness, without inferring thermal adaptation or an evolutionary change after introduction.

### Discussion short alternative (if only 80–100 words are available)

> Elevation was treated as a geographically contingent environmental proxy rather than an independent climatic driver. The lower photographic corolla lightness of introduced *C. arvense* at higher elevations persisted after adjustment for BIO18, BIO05 and RSDS and under spatial-block resampling, consistent with but not diagnostic of temperature-mediated colour variation (Körner, 2007; Lacey & Herr, 2005; Berardi et al., 2016). Elevational gradients also reflect local radiation, snow, soils, population structure and sampling conditions, which differ across latitudes and ranges. Image lightness does not directly measure reflectance or floral temperature. Separate native and introduced coefficients do not demonstrate different responses to elevation.

## 追加するなら最小限の検証（未実施、今回の本文の前提条件ではない）

1. 標高の元の移入域標本で、地理的大区分（生態地理区・山地／島域）別または地域内中心化した標高偏差の符号と支持数を可視化し、*地域間の高低差*と*同じ地域内の上下差*を分ける。
2. *C. arvense*在来・移入域の共通の**標高と気候支持域**で、両域共通のX尺度による `L* ~ elevation * introduced + BIO18 + BIO05 + RSDS + date + spatial terms` を推定。支持不足なら推定不能と表示。
3. 花冠の分光反射（可視・UV）・温度を追加測定するまでは熱的機能・UV保護を示したと書かない。

## References checked

- Körner, C. (2007). The use of ‘altitude’ in ecological research. *Trends in Ecology & Evolution*, 22:569–574. https://doi.org/10.1016/j.tree.2007.09.006
- Lacey, E. P. & Herr, D. (2005). Phenotypic plasticity, parental effects, and parental care in plants? I. An examination of spike reflectance in *Plantago lanceolata*. *American Journal of Botany*, 92:920–930. https://doi.org/10.3732/ajb.92.6.920
- Berardi, A. E., Fields, P. D., Abbate, J. L. & Taylor, D. R. (2016). Elevational divergence and clinal variation in floral color and leaf chemistry in *Silene vulgaris*. *American Journal of Botany*, 103:1508–1523. https://doi.org/10.3732/ajb.1600106
- Lacey, E. P. (2026). Temperature and the evolution of flower color: A review. *American Journal of Botany*, 113:e70106. https://doi.org/10.1002/ajb2.70106
- McVicar, T. R. & Körner, C. (2013). On the use of elevation, altitude, and height in the ecological and climatological literature. *Oecologia*, 171:335–337. https://doi.org/10.1007/s00442-012-2416-7

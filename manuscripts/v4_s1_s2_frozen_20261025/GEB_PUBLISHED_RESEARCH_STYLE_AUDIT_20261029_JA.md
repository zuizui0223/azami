# GEB掲載Research Articleの構成をV4へ適用するための比較監査（2026-10-29）

**対象:** *Global Ecology and Biogeography* (GEB) に実際に掲載された、広域植物形質研究・在来／外来植物研究。出版物の文章を複製せず、**研究の問い・実測結果の提示順・主張の程度・図表の機能**のみを参考にする。V4原資料の係数・標本数・Holm多重検定群は変更しない。

## 1. 比較した採択論文

| GEB掲載論文・原典 | Abstract / Introduction / Results / Discussion の特徴 | V4で採用する編集原則 |
|---|---|---|
| Kuppler et al. (2020), [*Global gradients in intraspecific variation in vegetative and floral traits are partially associated with climate and species richness*](https://doi.org/10.1111/geb.13077) | 種内形質変異がなぜ重要か→温度・降水など地理的勾配→標本単位と仮説を明示。Aimは「どんな共通則か」を先に問う | 「種平均・カテゴリは地理的な連続形質変異を隠す」というギャップを、写真の数より先に置く。分類群間と分類群内を混同しない |
| Hanz et al. (2022), [*Climatic and biogeographical drivers of functional diversity in the flora of the Canary Islands*](https://doi.org/10.1111/geb.13507) | 系統・地域史、気候と形質の生態学的質問を一貫してつなぐ。導入・在来の地理文脈の解釈を分け、Resultsには観測値を置く | 日本・北米の複数のアザミ放散を比較群の選定理由に置く。ただしV4は独立適応放散や系統的な反復進化を**検定していない** |
| Grenié et al. (2025), [*Gaps in Global Alien Plant Trait Data and How to Fill Them*](https://doi.org/10.1111/geb.70131) | 全球形質データの欠測・地理的偏りという研究上の具体的な欠点を先に説明。データ構築と再現性が研究成果そのものを支える。Introductionから実測単位が明確 | 公開写真による反復した連続形質計測が種代表値に比べ何を増やすかを明確にする。同時に源泉観測数≠撮影枚数≠計測可能な頭花数と明記。写真の標本・計測偏りは隠さない |
| Riera et al. (2026), [*From Purpose to Persistence: How Intentional Introductions by Humans Influence Ecological Niches of Non-Native Plants*](https://doi.org/10.1111/geb.70245) | 導入を環境障壁の通過という生態学的な問いへ直接つなぐ。Aimと主要な比較・標本規模を冒頭で明確化 | S2の在来・導入域比較は「異なる分布史で形質―環境関係が維持／変化するか」という自然実験**的**対照。環境勾配の傾き差自体は未検定であり、導入後進化の証明ではない |
| Wiethase et al. (2024), [*Spatio-temporal integrated Bayesian species distribution models reveal lack of broad relationships between traits and range shifts*](https://doi.org/10.1111/geb.13819) | 候補形質に関する弱い・支持されない結果を中心的発見に偽装せず、サンプリングと解析の欠点を考慮した結果を提示 | V4の画像形質間の r≈0.2 前後は探索的な記述としてSIへ移し、Holm通過した環境関連を本文の主役にする |

## 2. GEB編集上の方針：単なる「Part 1/2/3」の書き方から脱却

1. **Introduction:** 広域の形態的多様性は、地域史と気候にどう関係するか→種代表値・カテゴリーの限界→なぜ *Cirsium*（日本・北米の地域的放散と豊富な公開写真）→具体的な降水・色彩勾配仮説→導入による新しい地域史を使う比較の問い→本研究の目的。
2. **Methods:** 全写真への計測成功を示唆しない。YOLO頭花局在化→連続角度・色・2D輪郭→品質判定、WCVP在来／導入、S1/S2別母集団、4環境変数、分類群間・内推定、多重比較数を再現可能に。
3. **Results:** ①観測可能範囲と計測数→②**S1 BIO18×角度（分類群間、Holm P0.00920635）**→③**S2 導入域の花色×日射／標高（Holm通過3件）**→④弱い形質間関連はSIへの短い参照。解釈や否定命題をResultsで繰り返さない。
4. **Discussion:** 観測結果から開始し、次に雨による花粉の濡れ／流失回避という具体的な候補を先行実験と照合。在来／導入は既存変異、可塑性、進化、検出力の競合理論を**実際の在来係数を見て**解釈。限界は最後に集約し、実証した因果として書かない。
5. **Main displays:** Holm支持のBIO18角度・日射／標高色が主図。測定工程・全球分布地図はデータ裏付けと権利確認が必要な未完了素材。28ペアのRV、形質間相関、8組合せはSIのまま。旧Figure 2の3弱相関は**Main Figureではない**。

## 3. 変更しない科学的境界

- 固定V4 S1：313,350元在来観測、角度168,907件／178分類群、分類群間β=+0.045446843、Holm P=0.0092063506（80検定中1件）。
- 固定V4 S2：297,622元野生記録（在来141,014、導入156,608）、2種4層、128検定中3件（導入 *C. arvense* 標高×明度、導入 *C. vulgare* 日射×彩度／明度）。
- native *C. vulgare* 日射×彩度β=+0.136883、introduced β=+0.134504。Holmの判定差は**効果の域間差を意味しない**。
- オリジナルのSHA256 `6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205`、S1 80/S2 128、多重検定補正、全28ペア・SIでの保全、画像角度は重力校正ではない点を維持。

## 4. 投稿規程

GEBの公式[Author Guidelines](https://onlinelibrary.wiley.com/page/journal/14668238/homepage/forauthors.html)はResearch Articleに **Aim / Location / Time period / Major taxa studied / Methods / Results / Main conclusions** の7項目・300語以内の構造化要旨を要求し、おおむね本文5000語、文献50、図表6–8点を目安とする。研究の**一般的な生態学的意義**がAbstractと主要図で読み取れる必要があり、データ・コード利用可能性と匿名レビューを要求する。引用は著者–年、参考文献は著者名順。実際の画像権利と匿名コード公開は未完了なので、**投稿準備完了とはしない**。

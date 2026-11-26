# 第3段階：気候密度のカバー率と形質変化（2026-11-26）
Run [45755106891](https://github.com/zuizui0223/azami/actions/runs/45755106891), commit 2b93941a8e45989baaf6a483cc33ef5480ee00f5, artifact 15087113735. ZIP SHA256 8bda36e924b4119984d19d6696c7f00d358b2086c7909d9a0251b77be506fc16 verified against downloaded ZIP. 240 mass summaries, 20 climate contrasts, 720 slopes. No missing p values; independent Holm recomputation agrees to <1e-14.

## 固定気候質量の格子粗さ感度（記述）
| Species | Trait | Representation | Factor 1 | 2 | 4 | 8 |
|---|---|---|---:|---:|---:|---:|
| C. vulgare | chroma | median | .0312 | 1.0514 | 2.6702 | 1.2441 |
| C. vulgare | chroma | record_equalcell | .3643 | 1.2512 | 2.7906 | 1.1307 |
| C. arvense | angle | median | -4.9331 | -8.9011 | -11.4147 | -6.3513 |
| C. arvense | angle | record_equalcell | -5.1060 | -8.9023 | -11.7816 | -6.5877 |
| C. arvense | outline | median | .0219 | .0290 | .0431 | .0174 |

C. vulgare excess mass coverage .15885/.35620/.57532/.83602; corresponding ESS 87.90/78.35/60.99/51.72. C. arvense excess mass coverage .21255/.46734/.84071/.95315; ESS 78.15/43.60/23.43/35.59. Higher coverage by coarsening does not provide more independent observations: assumption is within-coarse-bin representativeness. C. arvense lower angle tail rises consistently, C. vulgare higher chroma tail rises consistently. These are descriptive and reuse the same underlying measurements. Outline is rotation dependent.

## 5 climate axes and matched phenotype contrasts
720 rows (two species, two novelty definitions, two representations, two adjustment specifications, three traits, three summaries, five climate variables). Eight Holm<.05 in clustered t inference; zero in approximate block-score sign inference (minimum Holm .137329). Strongest t results include C. vulgare low outline frequency versus mean temperature, C. arvense high outline frequency versus annual rain; these are image indices. C. vulgare high chroma frequency versus annual rain in density-surplus record/date-count adjusted model beta+.022514 per native climate SD, Holm .049460 but block-sign Holm 1. C. arvense high chroma frequency versus BIO15 in 95% median climate-only beta -.494159, Holm .049240 but block-sign Holm .351562 (30 pairs/9 blocks). No robust individual climate-axis association across inference methods. Not evidence for dryness driving chroma or causally identified selection. Keep the full comparison visible, do not choose specifications by p.

## Decision
Two exploratory candidates (C. vulgare chroma; C. arvense angle) persist descriptively across grid choices. Next wave scans all eight archived measurable constructs from the source using fixed original pairs, native-reference ranks/circular hue, and one distributional omnibus per species/endpoint. The new scan is exploratory, not independent confirmation; do not use it to retroactively promote a favorable test. No genetic or fitness evidence.
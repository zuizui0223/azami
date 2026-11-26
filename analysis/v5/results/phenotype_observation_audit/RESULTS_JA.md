# 第8段階：既報の気候希少度–彩度関連の観測監査（2026-11-26）
Run [45770028974](https://github.com/zuizui0223/azami/actions/runs/45770028974), commit 495a1adf400d4d274ebb2ad2220f2c280fc8b12c, artifact 15093926331. ZIP SHA-256 2262c617b07efc27d7348f657efc51240863c96c61734df18f5b4e37bb6b52cb verified against downloaded bytes. Independent recomputation of Holm for sensitivity (within-specification 4 planned comparisons) and coverage (12 planned comparisons), both t and sign families, agreed within 3.1e-16. Original four coefficients reproduced to 1e-8 by the execution gate.

## C. vulgare
| Estimand | Specification | n cells | beta [pointwise 95% CI] | Holm t | Holm score-sign |
|---|---|---:|---|---:|---:|
| Introduced rarity per .1, C* | original | 1023 | .761 [.337,1.186] | .003467 | .0651 |
| | within region | 1023 | .344 [-.601,1.290] | 1 | 1 |
| | high observation support | 580 | .552 [-.012,1.116] | .164122 | .332 |
| | high support within region | 580 | .233 [-.151,.616] | .658983 | .8406 |
| Range×rarity, SD/SD | original | 2014 | .175 [.111,.239] | .000010 | .002 |
| | within block×range | 2014 | .046 [-.097,.190] | 1 | 1 |
| | high observation support | 1032 | .114 [-.003,.231] | .164122 | .1276 |
| | high support within block×range | 1032 | -.018 [-.176,.140] | .818903 | .8492 |

Original Holm here is from the new targeted 4-comparison audit, not the prior broader 8/6-comparison families; the latter remain the primary historical values. Region centering attenuates slopes by ~55% and ~73%; sign inference does not support within-region or high-support results. Different estimands: pooled geographic covariance is not an independent within-region response; attenuation does not prove that photography is the cause. Region-level source sorting, environment differences, demographic structure and measurement differences remain possible.

32 sensitivity rows: 24 estimated, 8 insufficient_support (all models including the proxy vector; all had enough cells/blocks, so inspect design rank; not proof of bias or robustness, no p-value interpretation for these 8). All complete-proxy baseline/within-region rows retained original cell counts. Need a rank diagnostic and a transparently reduced proxy model only if required to characterize data quality, not a search for significant adjusted effects.

## Measurement availability
80 descriptive coverage rows and 12 within-block rarity slopes for measured-row availability, usable C*, and 3-trait completeness; all Holm=1 for both inference methods. In common-background introduced C. vulgare, usable C* cell-equal fraction .62449 inside / .63424 outside native 95% core (959/190 cells); three complete .54541/.53537. Within-region C* availability slope per rarity .1 +.02747 [−.01515,.07009], corrected p=1. No detected large systematic linear coverage trend, not proof of equivalent coverage or calibrated colour. Native rare climate cells are scarce with lower descriptive coverage. Measurements cannot be validated without independent image review.

## Decision
Do not present the previously strong pooled association as robust independent evidence of adaptation. In particular, regional confounding substantially weakens the association. Keep phenotypic redistribution as an exploratory candidate and standing genetic variation (SGV) as a motivating mechanism, not established. Native prediction route remains stopped. No evidence of elapsed-time trait evolution, fitness or allele frequency changes.
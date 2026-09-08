# Hypervolume synthesis notes

The hypervolume layer is intentionally downstream of the primary ecological estimators. It is included for two pre-specified purposes only: (1) auditing whether successfully measurable observations occupy a restricted abiotic region relative to the full source cohort, and (2) summarizing multivariate phenotype breadth and centroid structure at the taxon level.

The primary implementation target is a shrinkage-covariance ellipsoid rather than a kernel-density hypervolume. This keeps the definition reproducible in moderate dimension, makes rarefaction and support thresholds explicit, and avoids choosing smoothing parameters after seeing ecological results.

Colour representation must respect circular hue and closed composition. Gross shape is limited to the four qualified endpoints. Fine-architecture and surface endpoints remain excluded because they did not pass the frozen measurement qualification gate.

Any source-versus-measurable environmental contraction is interpreted as potential image-assessability selection, not as a biological response. Phenotype-environment hypervolume correspondence remains associational and cannot establish plasticity, local adaptation, genetic differentiation, or mechanism.

# S1 discrete transitions implementation plan

**Goal:** Apply the verified EAzami discrete ancestry and transition-environment workflow to orientation, corolla colour and capitulum shape on the S1-linked nuclear panel.

**Architecture:** Freeze source-backed categorical states separately from continuous photographic measurements. Reuse EAzami's actual symmetric Mk likelihood, edge joint probabilities and count-preserving state-map null; preserve unresolved and polymorphic source states. Extend to unordered multistate traits only after their codebook and observations are verified.

**Tech stack:** Python, Bio.Phylo, NumPy/SciPy, pandas; pinned EAzami source c9870e91604ed1514344ddc8f765eae5237dbf27.

- [ ] Verify original code, source state crosswalk, hypothesis contract, branch-length and null definitions; reproduce a frozen EAzami reference.
- [ ] Create a 24-species evidence ledger for each requested categorical trait, exact source concepts and uncertainty. Annotation templates are not observations; do not replace capitulum shape with corolla throat or phyllary posture.
- [ ] Join present S1 common-cohort BIO1/BIO15 and four V4 environmental axes. Keep old occurrence/photo cohorts separate.
- [ ] Run source-backed orientation Mk ancestry, edge direction probabilities and minimum-step bounds on both nuclear topologies; compare low-recovery-tip exclusion where applicable. Label edge endpoint-change probability separately from within-edge jump counts.
- [ ] Extend validated equal-rates Mk calculations to source-backed colour and shape categories. Do not force a directional order onto nominal states or recode undefined categories to obtain results.
- [ ] Prespecify orientation BIO15-up/BIO1-down projection as an EAzami-motivated transfer hypothesis, not an independent confirmation if source species overlap. Test other axes/traits with explicit multiplicity families and count-preserving nulls, refitting rate per state map.
- [ ] Validate likelihood/posterior marginal consistency and small-tree exhaustive examples; export ancestor probabilities, branch joint matrices, all state maps or reproducible Monte Carlo settings, sensitivity and figures.
- [ ] Update T2/T3 status. Previous continuous analyses are supplementary; their completion does not complete this discrete goal. Any remaining evidence gaps stay explicit.

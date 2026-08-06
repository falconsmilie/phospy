# ADR: Stochastic Reproducibility and Adaptive Prediction Seed Governance

## Document Control

- **ADR ID:** ADR-0017
- **Title:** Stochastic Reproducibility and Adaptive Prediction Seed Governance
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines governance for stochastic behavior in scientific methods,
including explicit seed requirements, policy metadata provenance, stable
child-RNG derivation, and replayability expectations. Adaptive kinase
prediction remains one governed stochastic path.

This governance also applies to stochastic dataset preprocessing methods (for
example MinProb missing-data imputation): seed and policy assumptions must be
explicitly recorded for scientific replayability.

This governance also applies to activity-method stochastic tests such as
ssGSEA-style substrate-set permutation p-values.

## Status

Accepted.

## Context and Problem Statement

Adaptive prediction includes stochastic sampling. Activity methods may also
include stochastic test procedures, for example seeded substrate-set
permutation p-values. Results can change when seed and sampling/permutation
policy change. The code already enforces explicit seeds in adaptive mode and
records policy metadata/provenance fields. This must be ADR-governed as
scientific-output policy.

For stochastic tests indexed by scientific identity, iteration order must not
determine random streams. Sorting a loop is not sufficient governance because
adding a previously unrelated test can still advance a shared generator and
change existing named results. Stable semantic identifiers must derive child
RNG streams from the user seed.

## Decision Drivers

1. Scientific replayability and traceability.
2. Explicit stochastic governance over convenience defaults.
3. Clear distinction between stable and parity-oriented sampling policies.
4. Provenance completeness for policy-level audit.

## Decision

1. Stochastic scientific behavior must be explicit.
2. Adaptive kinase prediction requires an explicit seed.
3. Adaptive reproducible mode must record seed strategy in provenance.
4. Stable deterministic per-kinase sampling and R-parity/global-stream
   sampling are distinct policies.
5. Policy IDs, policy names, versions, parameters, and seed strategy must be
   provenance-visible.
6. Missing-seed errors are intentional and must not be downgraded to warnings.
7. Replayability is more important than convenience for stochastic scientific
   modes.
8. Deterministic scientific imputers that do not use RNG still require
   reproducibility regression coverage for repeated output, diagnostics, and
   provenance.
9. Stable stochastic scientific tests must derive deterministic child RNG
   streams from stable semantic identifiers and the user seed. A single
   sequential RNG shared by identity-indexed tests is not an acceptable stable
   seed strategy.
10. The ssGSEA substrate enrichment permutation policy derives each
    profile/kinase/method permutation stream from method ID, method version,
    profile ID, kinase name, stream name, and caller-supplied seed. Public
    provenance reports the active seed policy
    `stable_by_method_profile_kinase` and its version; the v1 hash encoding
    retains a private compatibility salt solely to preserve existing
    deterministic streams.

## Consequences

1. Users must provide `random_state` for adaptive reproducible prediction.
2. Changing the sampling policy may legitimately change results.
3. Policy version changes must be treated as scientific-output changes.
4. Parity mode exists for parity testing, not as the default recommendation.
5. Golden/provenance tests must cover seed and policy metadata.
6. Stochastic imputation paths must serialize seed, method parameters, and
   imputation summary diagnostics.
7. Deterministic imputation paths must keep tie-breaking deterministic and must
   be covered by replay-oriented tests even when no seed is recorded.
8. Reordering kinases or profiles must not change named stochastic activity
   results.
9. Adding an unrelated kinase must not change existing named permutation
   p-values. Multiple-testing adjustments may legitimately change when the
   tested family changes and must be interpreted separately from RNG-stream
   identity.

## Affected Modules

- `src/phospy/api/configs/prediction.py`
- `src/phospy/science/prediction/policies.py`
- `src/phospy/science/prediction/sampling_runtime.py`
- `src/phospy/science/prediction/sampling_core.py`
- `src/phospy/science/prediction/execution.py`
- `src/phospy/science/activities/methods/ssgsea_substrate_enrichment.py`
- `src/phospy/science/activities/scientific_policies.py`
- `src/phospy/workflows/kinase/provenance.py`
- `tests/unit/test_prediction_adaptive_sampling.py`
- `tests/unit/test_activity_science.py`
- `tests/parity/test_adaptive_replay_parity.py`
- `tests/fixtures/public_workflow_reference/kinase_public_predmat_provenance_golden.json`

## Scope Boundaries

This ADR governs stochastic reproducibility and seed-policy provenance for
scientific methods. It does not define general test suite structure (ADR-0014)
or public namespace governance (ADR-0001). It does not require order-invariant
streams for intentionally parity-oriented global-stream modes when those modes
are explicitly named, versioned, and recorded as such.

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. Adaptive stochastic paths reject missing seeds.
2. Selected sampling policy and seed strategy remain provenance-visible.
3. Policy metadata (`id`, `name`, `version`, `parameters`) remains serialized in
   provenance.
4. Stable and parity sampling policies remain explicitly distinct.
5. Policy or version changes are reviewed as scientific-output changes.
6. Identity-indexed stable stochastic tests use child RNG streams derived from
   semantic identifiers and the user seed.
7. Reversed mapping order, profile reordering, unrelated-test insertion,
   repeated runs, seed divergence, and input serialization round trips are
   covered for stochastic activity permutations.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR

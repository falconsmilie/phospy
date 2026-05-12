# ADR: Stochastic Reproducibility and Adaptive Prediction Seed Governance

## Document Control

- **ADR ID:** ADR-0017
- **Title:** Stochastic Reproducibility and Adaptive Prediction Seed Governance
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines governance for stochastic behavior in adaptive kinase
prediction, including explicit seed requirements, policy metadata provenance,
and replayability expectations.

This governance also applies to stochastic dataset preprocessing methods (for
example MinProb missing-data imputation): seed and policy assumptions must be
explicitly recorded for scientific replayability.

## Status

Accepted.

## Context and Problem Statement

Adaptive prediction includes stochastic sampling. Results can change when seed
and sampling policy change. The code already enforces explicit seeds in adaptive
mode and records policy metadata/provenance fields. This must be ADR-governed
as scientific-output policy.

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

## Consequences

1. Users must provide `random_state` for adaptive reproducible prediction.
2. Changing the sampling policy may legitimately change results.
3. Policy version changes must be treated as scientific-output changes.
4. Parity mode exists for parity testing, not as the default recommendation.
5. Golden/provenance tests must cover seed and policy metadata.
6. Stochastic imputation paths must serialize seed, method parameters, and
   imputation summary diagnostics.

## Affected Modules

- `src/phospy/api/configs/prediction.py`
- `src/phospy/prediction/policies.py`
- `src/phospy/prediction/sampling_runtime.py`
- `src/phospy/prediction/sampling_core.py`
- `src/phospy/prediction/execution.py`
- `src/phospy/workflows/kinase/provenance.py`
- `tests/unit/test_prediction_adaptive_sampling.py`
- `tests/parity/test_adaptive_replay_parity.py`
- `tests/fixtures/public_workflow_reference/kinase_public_predmat_provenance_golden.json`

## Scope Boundaries

This ADR governs stochastic reproducibility and seed-policy provenance for
adaptive prediction. It does not define general test suite structure
(ADR-0014) or public namespace governance (ADR-0001).

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. Adaptive stochastic paths reject missing seeds.
2. Selected sampling policy and seed strategy remain provenance-visible.
3. Policy metadata (`id`, `name`, `version`, `parameters`) remains serialized in
   provenance.
4. Stable and parity sampling policies remain explicitly distinct.
5. Policy or version changes are reviewed as scientific-output changes.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR

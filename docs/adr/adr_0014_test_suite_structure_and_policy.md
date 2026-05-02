# ADR: Test Suite Structure and Policy for PhosPy

## Document Control

- **ADR ID:** ADR-014
- **Title:** Test Suite Structure and Policy for PhosPy
- **Status:** Accepted
- **Date:** 2026-05-02
- **Decision Type:** Architecture Decision Record

## Abstract

This ADR defines governance for test categories, golden/provenance regression
policy, property-based invariants, and performance contracts.

## Status

Accepted.

This ADR supersedes earlier narrower test-structure-only guidance.

## Context and Problem Statement

The codebase now includes:

- golden fixtures for representative public workflow outputs
- provenance golden contracts
- parity tests with distinct purpose from golden tests
- dedicated performance contract tests under `tests/performance/`

These are architecture governance checks and release policy checks. They are
not ordinary local default unit-test behavior and are not manual-only checks.

## Decision Drivers

1. Keep architectural and scientific regression signals explicit.
2. Separate parity, golden, and performance purposes.
3. Preserve deterministic replay-critical metadata in provenance tests.
4. Make release blocking policy explicit for scientific production quality.

## Decision

### Test Category Structure

PhosPy test categories remain:

- `tests/unit`
- `tests/integration`
- `tests/parity`

`parity` remains opt-in at local default pytest invocation.

### Golden and Provenance Regression Governance

1. Golden tests are architecture-governance tests, not ordinary snapshot tests.
2. Golden fixtures should lock stable reproducibility-critical outputs.
3. Golden fixtures should not lock incidental runtime/environment details.
4. Provenance golden tests must cover:
   - workflow name
   - workflow config
   - reference metadata
   - reference fingerprints
   - output fingerprints
   - policy IDs
   - policy names
   - policy versions
   - seed strategy where stochastic behavior is used
   - score-preconditioning policy where applicable
5. Fixture updates require intentional review because they may represent
   scientific-output changes.
6. Parity tests and golden tests have different purposes and must not be merged
   conceptually.

### Property-Based Testing Governance

1. Property-based tests should protect scientific invariants.
2. Property-based tests should not assert incidental implementation order unless
   order is the invariant.
3. Property-based tests complement golden and parity tests rather than replacing
   them.

### Performance Contract Governance

1. tests/performance are release-gate checks.
2. tests/performance remain excluded from the default local
   unit/integration run.
3. tests/performance are not manual-only checks.
4. `tests/performance` run in a dedicated CI/release job or explicit
   release-validation command.
5. Failures in `tests/performance` block release until fixed, waived, or the
   contract is intentionally updated.
6. Contract updates require updating both `docs/performance.md` and related
   test expectations.
7. Scale guardrails are part of supported-behavior governance, not incidental
   implementation limits.

### Why Release-Gate Is the Correct Policy

1. Performance tests can be environment-sensitive for every local/default run.
2. Manual-only policy is too weak for production-quality scientific software.
3. Release-gate policy is the required middle ground between the two extremes.

## Consequences

### Positive Consequences

- Golden/provenance fixtures become explicit contract assets.
- Performance regressions cannot silently ship.
- Property-based tests stay focused on real scientific invariants.
- Review criteria for fixture/contract updates become explicit.

### Negative Consequences

- Release pipelines must run dedicated performance contract jobs.
- Fixture updates require explicit scientific review discipline.

## Affected Modules

- `tests/fixtures/public_workflow_reference/`
- `tests/integration/test_kinase_workflow_integration.py`
- `tests/integration/test_signalome_workflow_integration.py`
- `tests/unit/test_scientific_invariants.py`
- `tests/parity/`
- `tests/performance/test_performance_contracts.py`
- `tests/support/performance_contracts.py`
- `docs/performance.md`
- `src/phospy/signalomes/clustering/scale_guards.py`
- `pyproject.toml`

## Scope Boundaries

This ADR governs testing and release-gate policy. It does not define public
API namespace ownership (ADR-001), validation ownership (ADR-007), or internal
module-splitting governance (ADR-010).

## Validation and Review Criteria

Future changes must satisfy all of the following:

1. Is this test in the correct category for its purpose?
2. Does the change preserve golden vs parity separation?
3. Are provenance golden fields still reproducibility/audit critical?
4. If performance contracts changed, were both docs and thresholds updated?
5. Is release-gate semantics preserved (not default local, not manual-only)?

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR*. GitHub repository.
https://github.com/PYangLab/PhosR

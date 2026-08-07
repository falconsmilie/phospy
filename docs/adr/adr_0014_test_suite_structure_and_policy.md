# ADR: Test Suite Structure and Policy for PhosPy

## Document Control

- **ADR ID:** ADR-0014
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
- public-boundary adversarial gates for exported API signatures, provenance
  binding, DataFrame ownership, and recursive JSON immutability
- external-consumer public API contract tests under `tests/contract/`
- parity tests with distinct purpose from golden tests
- dedicated performance contract tests under `tests/performance/`

These are architecture governance checks and release policy checks. They are
not ordinary local default unit-test behavior and are not manual-only checks.

## Decision Drivers

1. Keep architectural and scientific regression signals explicit.
2. Separate parity, golden, and performance purposes.
3. Preserve deterministic replay-critical metadata in provenance tests.
4. Make release blocking policy explicit for scientific production quality.
5. Prevent public-boundary regressions from hiding behind broad coverage volume.

## Decision

### Test Category Structure

PhosPy test categories remain:

- `tests/unit`
- `tests/integration`
- `tests/parity`
- `tests/architecture`
- `tests/contract`
- `tests/release`
- `tests/golden`
- `tests/performance`

`tests/contract` is an explicit external-consumer contract suite. It remains
outside the default local pytest invocation and is selected by `make
test-contract`, which is release-blocking. Contract tests must exercise public
construction and workflow entrypoints rather than validation or construction
internals.

`parity`, `tests/contract`, `tests/release`, `tests/golden`, and
`tests/performance` remain outside the default local pytest invocation and are
selected by explicit release targets.
PhosPy-owned release-validation science contracts that are not external parity
may live under `tests/science` and carry `release_gate` when they protect
adverse scientific cases.

### Public-Boundary Adversarial Governance

1. `tests/architecture/test_public_boundary_integrity.py` inventories
   `phospy.__all__`, `phospy.api.__all__`, and public result/evidence JSON-like
   model fields.
2. `tests/unit/test_public_boundary_adversarial.py` runs compact runtime probes
   for public signatures, dataset provenance binding, DataFrame ownership, and
   JSON immutability.
3. The installed-distribution verifier
   (`scripts/verify_installed_distributions.py`) installs the built wheel and
   sdist in isolated environments outside the checkout, imports from the
   installed package location, verifies bundled reference resources and
   SHA-256 values, and exercises representative public dataset, differential,
   kinase, and supported public-boundary contracts without importing
   repository tests or fixtures.
4. Installed public-boundary and resource failures block release.
5. Source-tree external-consumer contracts in `tests/contract/` prove that a
   separate-style consumer can build an analysis-ready dataset, run supported
   workflows, receive typed results, and query result tables by `site_key` using
   only supported imports.
6. Architecture tests statically reject unsupported imports in
   `tests/contract/`, including `phospy.science`, `phospy.validation`, private
   `phospy` module segments, and workflow implementation modules such as
   validators, interpreters, and executors.

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
   - missing-data imputation summaries (method, parameters, affected counts,
     rows removed, final missing count, and intensity-scale assumptions)
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

1. Bounded tests/performance contracts are release-gate checks.
2. tests/performance remain excluded from the default local
   unit/integration run.
3. Bounded tests/performance contracts are not manual-only checks.
4. `tests/performance` run in a dedicated CI/release job or explicit
   release-validation target (`make test-performance`).
5. Failures in bounded `tests/performance` contracts block release until fixed,
   waived, or the contract is intentionally updated.
6. Contract updates require updating both `docs/performance.md` and related
   test expectations.
7. Scale guardrails are part of supported-behavior governance, not incidental
   implementation limits.
8. Large machine-dependent workloads that are disproportionate for routine CI
   belong under `benchmarks/` as opt-in local benchmarks, not under `tests/`.

### Release Validation Command Policy

1. A documented release-gate command must run:
   - normal unit tests
   - integration tests
   - external-consumer public API contract tests through `make test-contract`
   - release-gated reproducibility/golden tests through `make test-release-gates`
   - threshold-bearing parity tests through
     `pytest tests/parity -m "parity and not parity_diagnostic" -s`
   - performance contract tests
   - strict documentation build through `make docs-build`
   - archive-level wheel/sdist metadata and packaged-reference checks
   - installed wheel/sdist verification outside the checkout
2. Fast local defaults may keep parity, release/golden, and performance suites
   out of the default marker selection.
3. Scientific parity failures are release-blocking.
4. Missing optional dependency failures in release-gate execution must include
   clear setup guidance in maintainer docs.
5. `make release-check` is the authoritative aggregate release command.
6. A collection-only selector audit must compare actual pytest node IDs and
   effective markers against the authoritative release targets so
   release-blocking tests cannot be missed silently.
7. CI must run the release-science selectors on every supported Python version
   for the non-parity suite, external-consumer contract suite,
   threshold-bearing parity suite, release/golden gates, and bounded performance
   contracts unless a performance waiver is documented.
8. CI must include a maintained minimum-dependency lane that uses a dedicated
   lower-bound constraint file rather than the current pinned CI constraints,
   runs `pip check`, and executes the non-parity plus release/golden selectors
   that do not require external scientific tools.
9. CI must not run opt-in local benchmarks that are explicitly documented as
   machine-dependent and non-release-blocking.
10. CI and publication workflows must run installed wheel/sdist verification
    against the single built artifact set on Python 3.11 and 3.12.
11. CI must run `mkdocs build --strict` or an equivalent strict documentation
    build so documentation warnings, including broken internal links, fail the
    build.
12. Release-policy tests must audit the effective command and import closure of
    release Make targets and CI/publication workflow commands so the prohibited
    50,000 x 48 workload cannot become release-reachable through renaming,
    simple arithmetic, helper modules, configuration objects, or marker changes.

### Why Release-Gate Is the Correct Policy

1. Performance tests can be environment-sensitive for every local/default run.
2. Manual-only policy is too weak for production-quality scientific software.
3. Release-gate policy is the required middle ground between the two extremes.

## Consequences

### Positive Consequences

- Golden/provenance fixtures become explicit contract assets.
- Bounded performance regressions cannot silently ship.
- Property-based tests stay focused on real scientific invariants.
- Review criteria for fixture/contract updates become explicit.
- Machine-dependent scale observations remain available locally without
  blocking releases on hosted-runner variability.

### Negative Consequences

- Release pipelines must run dedicated bounded performance contract jobs.
- Fixture updates require explicit scientific review discipline.

## Affected Modules

- `tests/fixtures/public_workflow_reference/`
- `tests/integration/test_kinase_workflow_integration.py`
- `tests/integration/test_signalome_workflow_integration.py`
- `tests/unit/test_scientific_invariants.py`
- `tests/parity/`
- `tests/contract/`
- `tests/architecture/test_public_boundary_integrity.py`
- `tests/architecture/test_public_consumer_contract_imports.py`
- `tests/unit/test_public_boundary_adversarial.py`
- `tests/performance/test_performance_contracts.py`
- `benchmarks/measure_release_scale_builder_differential.py`
- `scripts/verify_installed_distributions.py`
- `tests/science/test_differential_adverse_design_contracts.py`
- `tests/science/test_evidence_resolution_regression_fixtures.py`
- `tests/science/test_kinase_sparse_support_regression_fixtures.py`
- `tests/science/test_signalome_safety_regression_fixtures.py`
- `tests/support/performance_contracts.py`
- `tools/testing/release_selector_coverage.py`
- `docs/testing/public_boundary_invariant_checklist.md`
- `docs/performance.md`
- `docs/maintenance.md`
- `docs/testing/pytest_markers.md`
- `Makefile` (`release-check`, `test-release-gates`, `test-performance`,
  `test-contract`, `docs-build`, `verify-installed-distributions`)
- `.github/workflows/publish.yml` (release-gate enforcement)
- `.github/workflows/ci.yml` (supported-version release-science matrix and
  minimum-dependency lane plus strict documentation gate)
- `constraints/minimum.txt`
- `mkdocs.yml`
- `src/phospy/science/signalomes/clustering/scale_guards.py`
- `pyproject.toml`

## Scope Boundaries

This ADR governs testing and release-gate policy. It does not define public
API namespace ownership (ADR-0001), validation ownership (ADR-0007), or
internal module-splitting governance (ADR-0010).

## Validation and Review Criteria

Future changes must satisfy all the following:

1. Is this test in the correct category for its purpose?
2. Does the change preserve golden vs parity separation?
3. Are provenance golden fields still reproducibility/audit critical?
4. If bounded performance contracts changed, were both docs and thresholds
   updated?
5. Is release-gate semantics preserved (not default local, not manual-only)?
6. For stochastic preprocessing methods, are seed and imputation provenance
   fields explicitly tested?
7. Are supported public boundaries still covered by source and installed-artifact
   adversarial probes?
8. Does the selector coverage audit still prove that every release-blocking
   collected node is selected by at least one authoritative release target?
9. Do supported Python and minimum-dependency CI lanes still reflect the
   declared support policy?
10. Are opt-in local benchmarks kept outside pytest collection, release-check,
    and CI by effective command/source reachability rather than only by file or
    target names?
11. Does installed wheel/sdist verification still install both artifacts
    outside the checkout, use isolated Python execution, verify bundled
    resources, exercise the supported public boundary, and avoid repository
    test helpers?
12. Do `tests/contract/` tests still use public construction/workflow routes
    only, avoid unsupported imports, receive typed results, and cover `site_key`
    table lookup behavior?
13. Does CI still run strict documentation builds so broken internal links fail
    before release?

## Amendment: Optional Release-Scale Benchmark Ownership (2026-07-29)

The 50,000-site x 48-sample end-to-end workload is no longer a release gate or
CI responsibility. It is retained as an opt-in local benchmark because its
runtime and memory cost are disproportionate for routine hosted CI.

The workload now belongs under `benchmarks/` and is invoked explicitly through
`make benchmark-release-scale`. It must not live under `tests/`, must not be
selected by `make test-performance`, must not be included in `make
release-check`, and must not be run by GitHub Actions. Its runtime and process
memory observations are informational local capacity data rather than release
budgets.

Release-policy validation audits effective 50,000 x 48 dimensions and command
reachability rather than only the historical benchmark target or filename. The
audit covers simple constant arithmetic, generic aliases, positional and keyword
dimensions, dictionary/config objects, local helper imports, Make target
closures, and CI/publication workflow command reachability. Bounded 50,000 x 12
and 50,000 x 24 contracts remain allowed.

## Amendment: Python 3.11 Minimum Support Policy (2026-08-03)

PhosPy's active supported interpreter matrix is Python 3.11 and Python 3.12.
Python 3.10 is no longer a supported CI, publication, minimum-dependency, or
installed-distribution verification target. The minimum-dependency lane remains
required, but it now runs on Python 3.11 as the lowest supported interpreter.

## References

Yang, P., Patrick, E., Humphrey, S. J., Ghazanfar, S., James, D. E., Jothi,
R., & Yang, J. Y. H. (2019). Kinase activity inference from quantitative
phosphoproteomics data using multiple linear models. *Bioinformatics, 35*(14),
i349-i356.

YangLab. (n.d.). *PhosR* [Computer software]. GitHub.
https://github.com/PYangLab/PhosR

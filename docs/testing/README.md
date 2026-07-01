# Testing Audit and Consolidation Docs

This folder is the entry point for the testing audit system used to inventory tests, protect high-value contracts, and guide safe test-suite consolidation.

## Purpose

The testing audit system exists to:

- define what test quality and coverage expectations are,
- map current tests to behavioral contracts,
- identify what must not regress during consolidation,
- and separate machine-generated evidence from human review decisions.

## Release-Gate Policy

Default `pytest` excludes parity tests for developer speed through the
`pyproject.toml` addopts setting. Public releases must run
`make test-release-gate`; that target is the release-blocking command and
includes the default non-parity suite, provenance/golden checks, reference
manifest gates, threshold-bearing parity tests, and performance contracts.

`parity_diagnostic` tests are informational unless they are intentionally
promoted into the release selector.

## Read This First

1. [Test Audit Rubric](./test_audit_rubric.md): scoring and decision criteria used across the audit.
2. [Pytest Markers](./pytest_markers.md): marker taxonomy used by the inventory and downstream grouping.

## Core Audit Artifacts

- [Test Inventory (Markdown)](./test_inventory.md): human-readable inventory of discovered tests.
- [Test Inventory (CSV)](./test_inventory.csv): structured inventory for filtering, scripting, and diffing.
- [Protected Test Contracts](./protected_test_contracts.md): explicit list of contracts that must remain covered through consolidation.
- [Consolidation Coverage Map](./consolidation_coverage_map.md): maps contracts and existing tests to consolidation targets/gaps.

## Generated Analysis Reports

These reports are generated from static analysis scripts and provide focused slices of risk/opportunity:

- [Validator Test Pattern Report](./validator_test_pattern_report.md)
- [DataFrame Ownership Test Report](./dataframe_ownership_test_report.md)
- [Diagnostic Payload Test Report](./diagnostic_payload_test_report.md)
- [Orchestration Test Candidates](./orchestration_test_candidates.md)

## Human Review Documents

- [Diagnostic Payload Review](./diagnostic_payload_review.md): reviewer decisions for diagnostic candidates; rewrite/narrow work should follow this decision document, not only the generated candidate report.
- [Orchestration Test Review](./orchestration_test_review.md): reviewer triage and decisions on candidate consolidation paths.

## Relationship Between Documents

1. The [rubric](./test_audit_rubric.md) defines evaluation standards.
2. The [inventory](./test_inventory.md) and [CSV](./test_inventory.csv) capture what exists.
3. [Protected contracts](./protected_test_contracts.md) identify behavior that cannot be dropped.
4. The [coverage map](./consolidation_coverage_map.md) connects existing tests/contracts to planned consolidation.
5. Generated reports provide domain-specific candidate evidence (input to review), not final rewrite authority.
6. Human review documents record final consolidation decisions where judgment is required (for diagnostics, use [diagnostic_payload_review.md](./diagnostic_payload_review.md) as the rewrite/narrow decision source).
7. Coverage/protected-contract documents remain the guardrails that final changes must satisfy.

## Static Audit Tools

`tools/testing/` is the home of static testing-audit scripts:

- `tools/testing/generate_test_inventory.py`
- `tools/testing/find_validator_test_patterns.py`
- `tools/testing/find_dataframe_ownership_tests.py`
- `tools/testing/find_diagnostic_assertion_clusters.py`
- `tools/testing/find_orchestration_test_candidates.py`

### Local Freshness Check

To verify generated audit docs are current (and match CI), run:

1. `python tools/testing/generate_test_inventory.py`
2. `python tools/testing/find_validator_test_patterns.py`
3. `python tools/testing/find_dataframe_ownership_tests.py`
4. `python tools/testing/find_diagnostic_assertion_clusters.py`
5. `python tools/testing/find_orchestration_test_candidates.py`
6. `git diff --exit-code -- docs/testing/test_inventory.csv docs/testing/test_inventory.md docs/testing/validator_test_pattern_report.md docs/testing/dataframe_ownership_test_report.md docs/testing/diagnostic_payload_test_report.md docs/testing/orchestration_test_candidates.md`

Only generated audit artifacts are checked. Human review documents (for example
`diagnostic_payload_review.md` and `orchestration_test_review.md`) are intentionally
excluded from the freshness gate.

## Recommended Review Order for Consolidation Work

1. [test_audit_rubric.md](./test_audit_rubric.md)
2. [pytest_markers.md](./pytest_markers.md)
3. [test_inventory.md](./test_inventory.md) and [test_inventory.csv](./test_inventory.csv)
4. [protected_test_contracts.md](./protected_test_contracts.md)
5. [consolidation_coverage_map.md](./consolidation_coverage_map.md)
6. Generated reports:
   - [validator_test_pattern_report.md](./validator_test_pattern_report.md)
   - [dataframe_ownership_test_report.md](./dataframe_ownership_test_report.md)
   - [diagnostic_payload_test_report.md](./diagnostic_payload_test_report.md)
   - [orchestration_test_candidates.md](./orchestration_test_candidates.md)
7. Human review decision documents:
   - [diagnostic_payload_review.md](./diagnostic_payload_review.md)
   - [orchestration_test_review.md](./orchestration_test_review.md)
8. Re-check [protected_test_contracts.md](./protected_test_contracts.md) and [consolidation_coverage_map.md](./consolidation_coverage_map.md) before merging consolidation changes.

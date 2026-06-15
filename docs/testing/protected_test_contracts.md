# Protected Test Contracts

This registry identifies high-risk test groups that must not be removed casually during
test consolidation work.

Use this file before proposing rewrites, consolidations, or deletions.

## Policy

- Protected does not mean immutable; it means changes require explicit contract-preserving justification.
- Any consolidation must preserve the protected risk/contract listed below.
- If a rewrite is proposed, reviewers should verify equivalent protection remains with
  clear before/after coverage mapping.

## Protected Groups

| Test path/group | Protected risk/contract | Consolidation allowed? | What must remain true after rewrite |
| --- | --- | --- | --- |
| `tests/parity/**` | Python outputs must remain aligned with accepted R/PhosR reference behavior on declared parity seams. | Limited, with explicit parity sign-off. | Reference-backed parity drift detection still covers all release-bearing seams and declared tolerances. |
| `tests/performance/**` | Runtime/memory envelopes for critical scientific workflows must not regress silently. | Limited, threshold tuning only with evidence. | Performance contracts still guard key hot paths, with representative workload coverage and deterministic thresholds. |
| Provenance/golden contract tests (for example `*provenance*`, `*bundle*`, `*snapshot*`, `*baseline*`, `test_quantitative_meaning_output_audit.py`) | Reproducibility-critical payloads/manifests/provenance fields must remain stable and machine-auditable. | Yes, if goldens/contracts are preserved or intentionally versioned. | Contract-critical fields, schema, and fingerprint/provenance semantics remain enforced end-to-end. |
| `tests/unit/test_frame_ownership_policy.py` (DataFrame ownership/copy-policy) | Public exports must be mutation-safe; owned/borrowed frame semantics and copy budgets must not regress. | Yes, but high-risk checks must remain targeted. | No caller-mutation leaks; public accessor snapshots remain defensive; borrow/owned alias guarantees stay explicit. |
| Workflow diagnostic boundary tests (for example `tests/unit/test_kinase_workflow_diagnostics.py`, `tests/unit/test_signalome_workflow_diagnostics.py`, `tests/unit/test_workflow_boundary_error.py`) | Public diagnostics/error payload contracts and provenance-bearing boundary context must remain stable enough for downstream consumers. | Yes, if public contract fields remain covered. | Required diagnostic/error fields, structured boundary context, and provenance/summary contract surfaces remain validated. |
| Maintainer script smoke coverage (for example `tests/unit/test_maintainer_scripts_smoke.py` over `scripts/active/**`, `scripts/support/**`, `scripts/run_pyright.py`) | Release infrastructure for parity fixture generation, provenance goldens, and public-workflow reference artefacts must remain import-safe and path-stable. | Yes, with explicit output-path and side-effect checks preserved. | Maintainer scripts still parse/import, avoid generation side effects on import, and keep documented output locations aligned to protected fixture/golden paths. |
| Validation-boundary ownership tests (for example `tests/unit/test_validator_boundaries.py`, `tests/unit/test_validation_ownership.py`, `tests/unit/test_domain_boundaries.py`) | Validation/interpreter/executor responsibilities must stay in the correct layer; boundary errors must remain explicit and deterministic. | Yes, with strict ownership-preservation review. | Validators continue owning validation rules; interpreters/executors do not absorb unrelated validation/science responsibilities. |
| Public workflow contract tests (for example `tests/integration/test_dataset_builder_integration.py`, `tests/integration/test_kinase_workflow_integration.py`, `tests/integration/test_signalome_workflow_integration.py`, `tests/unit/test_public_contract_*.py`) | Public API/workflow boundaries and published output contracts must remain backward-safe. | Yes, with contract-equivalence proof. | Public request/result surface and documented output contracts remain stable, with compatible error semantics. |
| Scientific invariant tests (for example `tests/unit/test_scientific_invariants.py`, `tests/unit/test_scientific_policies.py`, `tests/unit/test_activity_science.py`, `tests/unit/test_kinase_science.py`, `tests/unit/test_signalome_science*.py`, `tests/unit/test_prediction_science.py`) | Scientifically meaningful invariants must remain true even if internal implementation changes. | Yes, if invariant coverage is preserved or strengthened. | Invariant assertions remain direct and scientifically interpretable, not replaced by incidental implementation checks. |

## Explicit High-Risk DataFrame Ownership Paths

The following tests are high-risk and should not be removed without a direct equivalent:

- `tests/unit/test_frame_ownership_policy.py::test_builder_dataframe_copy_churn_regression_budget`
- `tests/unit/test_frame_ownership_policy.py::test_internal_borrowed_dataset_access_aliases_owned_frames_without_copy`
- `tests/unit/test_frame_ownership_policy.py::test_internal_borrowed_prediction_and_scoring_access_aliases_owned_frames`
- `tests/unit/test_frame_ownership_policy.py::test_signalome_validator_read_path_does_not_mutate_internal_frames`
- `tests/unit/test_frame_ownership_policy.py::test_signalome_interpreter_read_path_does_not_mutate_dataset_frames`
- `tests/unit/test_frame_ownership_policy.py::test_owned_construction_frames_can_be_mutated_after_owned_transfer`
- `tests/unit/test_frame_ownership_policy.py::test_safe_public_export_does_not_change_owned_provenance_state`
- `tests/unit/test_frame_ownership_policy.py::test_kinase_activity_result_series_properties_are_defensive_snapshots`

Regeneration note (TST-AUDIT-007, 2026-05-09): `docs/testing/dataframe_ownership_test_report.md`
still reports repeated ownership clusters after TST-OWN-001. Those remaining clusters are expected
because they preserve targeted high-risk protections (borrow-alias, mutation-safety, copy-budget,
and provenance-safe export behavior) rather than incidental duplication.

## Release Sensitivity

- `tests/parity/**`: release-sensitive, required for parity sign-off.
- `tests/performance/**`: release-sensitive, required for performance sign-off.

## Domain Ownership Boundary Guardrail

This registry is intentionally aligned to domain ownership boundaries:

- Validators own validation rules and boundary diagnostics.
- Interpreters own interpretation/planning concerns.
- Executors/workflows own orchestration and stage execution.

Test refactors must preserve these separations and must not relocate responsibilities
across layers without an explicit architecture decision.

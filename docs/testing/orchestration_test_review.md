# Orchestration Test Review (TST-AUDIT-004)

This review classifies private-member/orchestration candidate tests from
`docs/testing/orchestration_test_candidates.md` and records whether each case
was kept or rewritten.

Legend:
- `keep`: private seam is protecting a real contract (public/scientific/perf/provenance).
- `rewrite`: private coupling was incidental and replaced with public/interpreted behavior.

## Decisions By Candidate File

| File | Decision | Protected Contract / Reason |
| --- | --- | --- |
| `tests/integration/test_quantitative_meaning_output_audit.py` | keep | Provenance/bundle payload compatibility boundary (`_shared` serializers are the contract seam). |
| `tests/integration/test_signalome_bundle_integration.py` | keep | Bundle compatibility/provenance contract; `_from_owned` is fixture construction for payload variants. |
| `tests/integration/test_signalome_publishing_integration.py` | keep | Publishing/bundle output contract; private constructor used for controlled result-shape variants. |
| `tests/parity/test_kinase_workflow_parity.py` | keep | Parity-sensitive regression setup; direct mutation supports reproducibility checks. |
| `tests/performance/test_performance_contracts.py` | keep | Performance guard rails and scale/latency contracts. |
| `tests/performance/test_signalome_clustering_benchmark.py` | keep | Benchmark contract at clustering backend seam. |
| `tests/unit/test_boundary_constructor_no_repair.py` | keep | Boundary constructor ownership/no-repair contract. |
| `tests/unit/test_dataset_preprocessing_subsystem.py` | keep | Internal stage handoff aliasing/copy-budget contract. |
| `tests/unit/test_dataset_run_provenance.py` | keep | Provenance environment serialization fallback contract. |
| `tests/unit/test_dataset_site_sequence_resolution.py` | keep | Bundle schema compatibility boundary (`_shared` payload codec seam). |
| `tests/unit/test_dataset_transformation_state_establishment.py` | keep | Transformation-state provenance and bundle compatibility contract. |
| `tests/unit/test_frame_ownership_policy.py` | keep | High-risk ownership/aliasing/perf/copy-budget contracts (explicitly protected). |
| `tests/unit/test_kinase_workflow_components.py` | rewrite + keep | Rewrote incidental overlap helper use; kept `_from_owned` alias assertions for zero-copy ownership contract. |
| `tests/unit/test_kinase_workflow_diagnostics.py` | keep | Boundary diagnostics and activity seam behavior (scientific failure-mode reporting). |
| `tests/unit/test_prediction_sequence_validation.py` | rewrite | Replaced private scoring stage call with public workflow scoring result assertion. |
| `tests/unit/test_processing_state_bundle_payload.py` | keep | Bundle payload compatibility boundary (`_shared` serialization seam). |
| `tests/unit/test_public_contract_import_routes.py` | keep | Import-route compatibility contract for public compatibility adapters. |
| `tests/unit/test_runtime_contract_guards.py` | keep | Defensive internal contract guard not reachable via public request path. |
| `tests/unit/test_signalome_bundle_compatibility.py` | keep | Signalome bundle compatibility payload contract seam. |
| `tests/unit/test_signalome_workflow_diagnostics.py` | keep | Boundary/seam diagnostics for scientific and workflow failure modes. |
| `tests/unit/test_site_identifier_canonicalization.py` | rewrite | Replaced private overlap/index helper usage with interpreter/public contract assertions. |
| `tests/unit/test_validation_ownership.py` | keep | Validation ownership architecture contract (single-owner policy). |
| `tests/unit/test_validator_boundaries.py` | keep | Boundary-layer validation behavior (targeted setup for missing protein mapping cases). |

## Rewrites Performed

1. `tests/unit/test_prediction_sequence_validation.py`
   - Replaced `KinaseWorkflowExecutor()._run_scoring_stage(...)` with
     `KinaseWorkflow().run(request).scoring_result`.
2. `tests/unit/test_site_identifier_canonicalization.py`
   - Replaced private interpreter helper calls with `KinaseWorkflowInterpreter().run(...)`
     and public/interpreted overlap/index assertions.
3. `tests/unit/test_kinase_workflow_components.py`
   - Replaced private overlap helper call with overlap computed from interpreted/public
     tables.
4. `tests/unit/test_runtime_contract_guards.py`
   - Added explicit rationale comment for retained private-seam guard.
5. `tests/unit/test_kinase_workflow_components.py`
   - Added explicit rationale comment for retained private alias/copy-budget assertions.

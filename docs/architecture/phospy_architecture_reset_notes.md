# PhosPy Architecture Reset Notes

- Last reconciled: 2026-04-20
- Document role: reset-era architecture note with explicit historical markers

> Audience: maintainers and contributors reviewing architecture history and authority order.
> New users should start at [Getting started](../getting-started/index.md).

## Document Status

- **Active contract source**: accepted ADRs (`docs/adr/adr_0001` through
  `adr_0015`) plus implemented public contracts under `src/phospy/`.
- **This file status**: **Historical / superseded in part**. Keep as context for
  reset rationale; do not treat older prescriptive sections as higher authority
  than accepted ADRs.

## Source-of-Truth Order (Governance)

1. Direct maintainer instruction
2. Accepted ADRs
3. This reset note (historical context and still-valid principles only)
4. Existing code, only when not conflicting with higher sources

## Active Architecture Snapshot (Current Contract)

### Package boundary

- Supported package: `src/phospy/`
- Historical reference material: archived project snapshots in git history
- Legacy package structure is not a migration target (ADR-012).

### Public product shape

- One dataset model: `AnalysisReadyPhosphoDataset`
- One builder story: `AnalysisReadyDatasetBuilder.run(DatasetBuildRequest(...))`
- Two public workflows:
  - `KinaseWorkflow.run(KinaseWorkflowRequest(...))`
  - `SignalomeWorkflow.run(SignalomeWorkflowRequest(...))`

### Dataset contract (supersedes reset-era draft wording)

- Dataset model uses explicit `intensity_scale_state` and `processing_state`
  (ADR-006), not a broad `transformation_state` label.
- Required site metadata columns are `gene_symbol` and `site`; `site_sequence`
  is optional at final dataset boundary and validated when present (ADR-003
  update note; `src/phospy/validation/datasets/analysis_ready.py`).

### Workflow architecture

- Validator -> Interpreter -> Executor staging is active and implemented for
  builder, kinase workflow, and signalome workflow (ADR-002).
- Public workflows expose `run(request)` as the contract entrypoint.

### Result model contract

- Public results are nested typed containers (ADR-005):
  - `KinaseWorkflowResult` with nested scoring/prediction/activity outputs
  - `SignalomeWorkflowResult` with nested kinase lineage + signalome outputs
- No compatibility flattening or mirror accessor layer is part of the active
  contract.

## Scientific Contract Snapshot (Current Lane)

- Kinase scoring publishes authoritative `profile_scores` and `combined_scores`;
  `motif_scores` and `weights` are optional diagnostic outputs.
- Authoritative downstream matrix is `combined_scores` first with
  `profile_scores` fallback.
- Prediction and signalome consume the same resolved downstream score lane.
- `expanded_signalome` is an official supported signalome output and is
  populated in the supported executor lane (optional by type for compatibility).
- Detailed science status and follow-ons live in
  `docs/architecture/legacy_science_gap_audit.md`, including governance
  details for repaired kinase ranking comparison surfaces and active
  ranking parity gates in the supported L6 lane.

## Historical / Superseded Reset Content

The following reset-era ideas are now historical and must not be treated as
live contract guidance:

- Reset-era dataset wording that required `site_sequence` universally.
- Reset-era dataset wording that used `intensity_scale` as the final contract.
- "Immediate next steps" checklists that predate accepted ADRs and completed
  implementation milestones.
- Package-shape sketches presented as future intent rather than current truth.

These items are retained only to explain rewrite-phase thinking.

## Still-Valid Reset Principles

The following principles remain valid and are now enforced by ADRs/code:

- Fresh-start rewrite: reuse science where useful, do not preserve legacy
  architecture by default.
- Keep public API small and workflow-oriented.
- Prefer explicit boundaries and typed contracts.
- Reject wrapper-heavy or alias-heavy compatibility drift.

## ADR Cross-Links for Active Contract Decisions

- Public API contract: `docs/adr/adr_0001_public_api_contract.md`
- Workflow staging: `docs/adr/adr_0002_internal_workflow_architecture.md`
- Dataset boundary: `docs/adr/adr_0003-analysis_ready_dataset_and_preprocessing_boundary.md`
- Result model rules: `docs/adr/adr_0005__result_model_design.md`
- Intensity-scale and processing-state contract:
  `docs/adr/adr_0006_transformation_state_and_transformer_contract.md`
- Builder public contract: `docs/adr/adr_0011_builder_public_api_contract.md`
- Fresh-start rewrite strategy: `docs/adr/adr_0012_rewrite_roadmap_and_fresh_start_plan.md`
- Scientific parity policy: `docs/adr/adr_0013_scientific_parity_strategy_and_parity_testing_policy.md`

## Conflict Handling Rule

If this note and an accepted ADR diverge, ADR text governs. If code diverges
from accepted ADRs, document the conflict explicitly rather than normalizing it
through ambiguous wording.

## Where Next

- Current architecture navigation: [Architecture index](index.md)
- Decision records: [ADR Index](../adr/index.md)
- Public usage contract: [API Guide](../api.md)

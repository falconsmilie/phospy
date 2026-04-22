# Core Concepts

PhosPy is easiest to use when you keep five core boundaries in mind.

## 1. Analysis-Ready Dataset Boundary

`AnalysisReadyPhosphoDataset` is the strict workflow input boundary.

- Workflows consume this model, not raw files or loose frames.
- Structural/content validation is enforced at construction.

Deep detail: [Validation Guide](../validation.md).

## 2. Builder Boundary

`AnalysisReadyDatasetBuilder` is the public ingestion path.

- It accepts file paths or DataFrames.
- It applies supported preprocessing policies.
- It returns a strict `AnalysisReadyPhosphoDataset`.

Deep detail: [API Guide](../api.md#builder-contract).

## 3. Workflow Requests and Results

Each workflow is one request DTO in, one typed result DTO out.

- `KinaseWorkflowRequest -> KinaseWorkflowResult`
- `SignalomeWorkflowRequest -> SignalomeWorkflowResult`

Results are nested by stage (for example `result.scoring_result.profile_scores`).

Deep detail: [API Guide](../api.md#result-contract-nested-stage-outputs).

## 4. References Are Explicit Contract Inputs

Workflows depend on explicit reference resolution (`ReferencePreset` or `ReferenceBundle`).

- Bundled runtime references are intentionally narrow in the current release.
- Non-bundled lanes require explicit caller-provided references.

Deep detail: [API Guide](../api.md#reference-resolution).

## 5. Supported Science vs Governance Claims

PhosPy distinguishes:

- what is implemented
- what is supported
- what is parity-gated
- what is contract-changed relative to legacy behavior

Deep detail: [Parity to PhosR](../parity.md), [Legacy science gap audit](../architecture/legacy_science_gap_audit.md).

## Next Step

Choose a deeper path in [Choose your path](../learning-paths/choose-your-path.md).

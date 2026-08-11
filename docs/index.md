# PhosPy documentation

PhosPy is a Python package for selected phosphoproteomics workflows.

Use these docs when you want to build an analysis-ready phosphosite dataset,
choose the right workflow, understand the request fields, and interpret the
result tables without digging through implementation modules.

## Start here

1. [Install PhosPy](installation.md).
2. [Prepare a dataset](api/dataset-build-workflow.md).
3. [Run your first analysis](quickstart.md).
4. Choose a workflow page:
   [differential analysis](api/differential-analysis.md),
   [enrichment](api/enrichment.md),
   [kinase analysis](api/kinase.md), or
   [signalome analysis](api/signalome.md).

## Workflow map

| Goal | Workflow page | Public entry point |
| --- | --- | --- |
| Test explicit condition contrasts | [Differential analysis](api/differential-analysis.md) | `DifferentialAnalysisWorkflow.run(DifferentialAnalysisRequest)` |
| Run offline ORA on selected identifiers | [Enrichment](api/enrichment.md) | `EnrichmentWorkflow.run(EnrichmentWorkflowRequest)` |
| Score and predict kinase support | [Kinase analysis](api/kinase.md) | `KinaseWorkflow.run(KinaseWorkflowRequest)` |
| Summarize kinase outputs into modules and network-style tables | [Signalome analysis](api/signalome.md) | `SignalomeWorkflow.run(SignalomeWorkflowRequest)` |

The dataset builder is shared by the workflows. It returns an
`AnalysisReadyPhosphoDataset` with `site_key` row identity, required
`site_sequence` metadata, and auditable preprocessing state.

Duplicate `display_id` values remain valid when the corresponding `site_key`
values differ. Duplicate rows that resolve to the same `site_key` are a
scientific ambiguity and fail by default unless you deliberately choose and
audit a non-error duplicate-site policy.

## Scientific scope

PhosPy is PhosR-inspired, not a full PhosR replacement. Scientific scope and
parity/open-gap status live in [Scientific Coverage](scientific-coverage.md).
For practical caveats, start with
[Scientific interpretation and limitations](scientific-interpretation.md).

PhosPy does not provide HTTP endpoints or a hosted service. The supported user
interface is the Python API.

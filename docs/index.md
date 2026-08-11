# Welcome to PhosPy

PhosPy provides a focused Python API for phosphoproteomics analysis. These
guides are designed to help you move from an intensity table to an interpretable
result without having to study PhosPy's internal architecture.

## Start Here

1. [Install PhosPy](installation.md).
2. [Prepare an analysis-ready dataset](api/dataset-build-workflow.md).
3. [Run a first analysis](quickstart.md).
4. Choose the workflow that matches your scientific question.

## Choose a Workflow

| Scientific question | Guide | Public entry point |
| --- | --- | --- |
| Which phosphosites differ between named conditions? | [Differential Analysis](api/differential-analysis.md) | `DifferentialAnalysisWorkflow.run(...)` |
| Are selected identifiers over-represented in my local sets? | [Enrichment](api/enrichment.md) | `EnrichmentWorkflow.run(...)` |
| Which kinases have the strongest support for each site? | [Kinase Analysis](api/kinase.md) | `KinaseWorkflow.run(...)` |
| How can kinase-supported sites be summarized into modules and associations? | [Signalome Analysis](api/signalome.md) | `SignalomeWorkflow.run(...)` |

Each workflow guide contains its complete request model, a runnable example,
its response format, interpretation guidance, and common problems.

## Understand the Dataset Boundary

The dataset builder returns an `AnalysisReadyPhosphoDataset`. Rows use a unique
`site_key`; `display_id` remains a readable label and may repeat when protein
context differs.

Duplicate `display_id` values remain valid when the corresponding `site_key`
values differ. Duplicate rows that resolve to the same `site_key` are a
scientific ambiguity and fail by default unless you deliberately choose and
audit a non-error duplicate-site policy.

## Keep the Scientific Limits in View

PhosPy is inspired by selected PhosR workflows, but it does not claim full PhosR
equivalence. Start with [Scientific Interpretation and
Limitations](scientific-interpretation.md) for practical guidance. The detailed
support matrix is available in [Scientific Coverage](scientific-coverage.md).

PhosPy does not provide HTTP endpoints or a hosted service. The supported user
interface is the Python API.

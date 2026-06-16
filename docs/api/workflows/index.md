# Workflow API Pages

Each public PhosPy workflow has a dedicated API page. Most workflows start with
a built `AnalysisReadyPhosphoDataset`; enrichment is independent and starts from
caller-supplied identifiers, sets, and background.

| Workflow | Public entrypoint | Request | Result |
| --- | --- | --- | --- |
| [Differential analysis](differential-analysis.md) | `DifferentialAnalysisWorkflow` | `DifferentialAnalysisRequest` | `DifferentialAnalysisResult` |
| [Enrichment](enrichment.md) | `EnrichmentWorkflow` | `EnrichmentWorkflowRequest` | `EnrichmentWorkflowResult` |
| [Kinase](kinase.md) | `KinaseWorkflow` | `KinaseWorkflowRequest` | `KinaseWorkflowResult` |
| [Signalome](signalome.md) | `SignalomeWorkflow` | `SignalomeWorkflowRequest` | `SignalomeWorkflowResult` |

`DifferentialAnalysisWorkflow`, `KinaseWorkflow`, and `SignalomeWorkflow` are
also available from top-level `phospy`. `EnrichmentWorkflow` is public from
`phospy.api` and `phospy.workflows`.

The dataset builder is documented separately because it is the strict
analysis-ready dataset boundary rather than a member of `phospy.api.workflows`.
See [Dataset Build API](../dataset-build-workflow.md).

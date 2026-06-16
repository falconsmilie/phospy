# Workflow API Pages

Each public PhosPy workflow has a dedicated API page. Start with a built
`AnalysisReadyPhosphoDataset`, then choose the workflow that matches your
analysis question.

| Workflow | Public entrypoint | Request | Result |
| --- | --- | --- | --- |
| [Differential analysis](differential-analysis.md) | `DifferentialAnalysisWorkflow` | `DifferentialAnalysisRequest` | `DifferentialAnalysisResult` |
| [Enrichment](enrichment.md) | `EnrichmentWorkflow` | `EnrichmentWorkflowRequest` | `EnrichmentWorkflowResult` |
| [Kinase](kinase.md) | `KinaseWorkflow` | `KinaseWorkflowRequest` | `KinaseWorkflowResult` |
| [Signalome](signalome.md) | `SignalomeWorkflow` | `SignalomeWorkflowRequest` | `SignalomeWorkflowResult` |

The dataset builder is documented separately because it is the strict
analysis-ready dataset boundary rather than a member of `phospy.api.workflows`.
See [Dataset Build API](../dataset-build-workflow.md).

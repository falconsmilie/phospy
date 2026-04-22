# Workflow Guides

PhosPy has one supported workflow chain:

1. Build `AnalysisReadyPhosphoDataset`
2. Run `KinaseWorkflow`
3. Run `SignalomeWorkflow` (optional, downstream of kinase result)

## Core Workflow Pages

- Beginner walkthrough: [Quickstart: first workflow](../getting-started/quickstart-first-workflow.md)
- API shape and request/result contracts: [API Guide](../api.md)
- Validation rules at each boundary: [Validation Guide](../validation.md)
- CLI workflow execution: [CLI Guide](../cli.md)
- Persisting results to bundles: [Output Bundles](../output_bundles.md)

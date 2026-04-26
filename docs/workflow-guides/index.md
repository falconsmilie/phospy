# Workflow Overview

PhosPy supports one public workflow chain:

1. build `AnalysisReadyPhosphoDataset`
2. run `KinaseWorkflow`
3. optionally run `SignalomeWorkflow`

For most users, the quickstart already covers this path:

- [Quickstart](../getting-started/quickstart-first-workflow.md)
- [CLI Guide](../cli.md)
- [API Guide](../api.md)
- [Performance Contracts](../performance.md)
- [Troubleshooting](../getting-started/troubleshooting-first-run.md)

Performance note:

- quantile normalisation and signalome clustering are the two most common dense/runtime-sensitive steps,
- kinase scoring cost depends strongly on dataset/reference overlap and retained diagnostic scoring tables,
- see [Performance Contracts](../performance.md) for current thresholds and CI regression policy.

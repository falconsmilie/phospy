# What Is PhosPy?

PhosPy is a Python package for phosphoproteomics workflows with a strict public contract.

It is designed for users who want:

- a typed, validated analysis-ready dataset boundary
- a supported kinase workflow (`KinaseWorkflow`)
- a supported signalome workflow (`SignalomeWorkflow`)
- explicit scientific and contract documentation for what is supported now

## Who Is It For?

PhosPy is useful if you are:

- scientifically strong and want a predictable Python workflow surface
- technically strong and want clear API/validation boundaries for phosphoproteomics analysis
- maintaining a reproducible workflow and need explicit contract and parity documentation

## Supported Product Shape

The supported public shape is intentionally focused:

- Build analysis-ready data:
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`
- Run kinase workflow:
  `KinaseWorkflow().run(KinaseWorkflowRequest(...))`
- Run signalome workflow:
  `SignalomeWorkflow().run(SignalomeWorkflowRequest(...))`

For full contract details, see [API Guide](../api.md).

## What PhosPy Is Not

- It is not a broad utility toolbox with many loosely defined entrypoints.
- It does not claim blanket parity for all historical PhosR lanes.
- It does not hide architectural/scientific constraints; those are documented explicitly.

For parity and scientific confidence tiers, see [Parity to PhosR](../parity.md).

## Next Step

Continue with [Quickstart: first workflow](quickstart-first-workflow.md).

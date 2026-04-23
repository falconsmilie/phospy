# What Is PhosPy?

PhosPy is a Python package for phosphoproteomics workflows with a strict public contract.

In plain terms, PhosPy helps you move through one supported workflow chain:

1. start with quantified phosphosite data
2. build a validated analysis-ready dataset
3. run kinase scoring and prediction
4. optionally group results into signalome outputs

PhosPy does **not** expose a web server or REST API. You use it through Python
classes or the `phospy` CLI.

## Who Is It For?

PhosPy is useful if you are:

- a scientist who wants a reproducible, documented workflow surface
- an engineer who wants clear typed boundaries and validation rules
- a maintainer who needs explicit support, parity, and contract documentation

## Supported Product Shape

The supported public shape is intentionally focused:

- Build analysis-ready data:
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`
- Run kinase workflow:
  `KinaseWorkflow().run(KinaseWorkflowRequest(...))`
- Run signalome workflow:
  `SignalomeWorkflow().run(SignalomeWorkflowRequest(...))`

For full contract details, see [API Guide](../api.md).

## What You Need for a First Success

For the beginner-friendly supported lane, keep the setup small:

- a `phospho` matrix with site IDs like `TSC2;S939;`
- matching `site_metadata` with `gene_symbol` and `site`
- `organism=Organism.RAT` if you want bundled references through `ReferencePreset.AUTO`
- `protein_id` only when you plan to run signalome

## What PhosPy Is Not

- It is not a broad utility toolbox with many loosely defined entrypoints.
- It does not claim blanket parity for all historical PhosR lanes.
- It does not hide scientific or architectural limits; those are documented explicitly.

For parity and scientific confidence tiers, see [Parity to PhosR](../parity.md).

## Next Step

Continue with [Quickstart: first workflow](quickstart-first-workflow.md).

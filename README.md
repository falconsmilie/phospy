# PhosPy

PhosPy is a Python package for phosphoproteomics workflows with an explicit, typed public contract.

It provides:

- one analysis-ready dataset boundary
- one dataset builder story
- one kinase workflow story
- one signalome workflow story

## Start Here

Install:

```bash
pip install .
```

Run the official examples:

```bash
python examples/dataset_builder_demo.py
python examples/kinase_workflow_demo.py
python examples/signalome_workflow_demo.py
```

Then follow the guided docs path:

1. [Docs Home](docs/index.md)
2. [What is PhosPy?](docs/getting-started/what-is-phospy.md)
3. [Quickstart: first workflow](docs/getting-started/quickstart-first-workflow.md)
4. [Core concepts](docs/concepts/core-concepts.md)

## Minimal API Example

```python
from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import DatasetBuildRequest, KinaseWorkflowRequest

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho_df,
        site_metadata=site_metadata_df,
    )
)

result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
    )
)
```

## Import Contract

`phospy.api` is the canonical namespace where public API types are defined and
organized in source.

Both namespaces are public, but they have different roles:

- top-level `phospy` is a curated convenience surface for the four main product
  entrypoints only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`.
- import requests, configs, results, enums/references, and errors from
  `phospy.api`.

## Supported Public Product Shape

- Dataset builder:
  `AnalysisReadyDatasetBuilder().run(DatasetBuildRequest(...))`
- Dataset boundary:
  `AnalysisReadyPhosphoDataset`
- Kinase workflow:
  `KinaseWorkflow().run(KinaseWorkflowRequest(...))`
- Signalome workflow:
  `SignalomeWorkflow().run(SignalomeWorkflowRequest(...))`

For strict contract details, use [docs/api.md](docs/api.md) and [docs/validation.md](docs/validation.md).

## Documentation Map

- Beginner onboarding: [docs/getting-started](docs/getting-started/index.md)
- Tutorials and user guides: [docs/user-guides](docs/user-guides/index.md)
- Workflow guides: [docs/workflow-guides](docs/workflow-guides/index.md)
- API/reference: [docs/reference](docs/reference/index.md)
- Validation/contracts: [docs/contracts](docs/contracts/index.md)
- Architecture/ADRs: [docs/architecture](docs/architecture/index.md), [docs/adr](docs/adr/index.md)
- Scientific/parity/governance: [docs/science](docs/science/index.md)
- Contributor/maintainer docs: [docs/contributor](docs/contributor/index.md)

## Project Boundary and Scientific Scope

- Supported package boundary: `src/phospy/`
- Historical migration archive: `legacy_archive/phospy_legacy/`
- Scientific confidence and parity claims are tiered and explicit; see:
  - [docs/parity.md](docs/parity.md)
  - [docs/architecture/legacy_science_gap_audit.md](docs/architecture/legacy_science_gap_audit.md)

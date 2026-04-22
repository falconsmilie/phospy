# PhosPy

PhosPy is a Python package for phosphoproteomics workflows with a strict public
contract around one dataset boundary and two workflow entrypoints
(`KinaseWorkflow`, `SignalomeWorkflow`).

## Install

For normal package use:

```bash
pip install phospy
```

If you plan to read/write `.parquet` inputs or outputs, install the optional
parquet extra:

```bash
pip install "phospy[parquet]"
```

For local contributor installs from a clone:

```bash
pip install -e ".[dev]"
```

For editable installs with parquet support:

```bash
pip install -e ".[dev,parquet]"
```

## First Run (Supported Happy Path)

Bundled runtime references in this release are rat-only. The recommended first
workflow is:

1. build an analysis-ready dataset with `organism=Organism.RAT`
2. run kinase with `references=ReferencePreset.AUTO`
3. optionally run signalome on the kinase result

Minimum input tables:

- `phospho`: numeric site-by-sample matrix with canonical site IDs as row index
- `site_metadata`: rows aligned to `phospho.index`, with `gene_symbol` and `site`
- add `protein_id` if you plan to run `SignalomeWorkflow`

For human/mouse lanes, pass an explicit `ReferenceBundle` instead of `AUTO`.

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow, SignalomeWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
    SignalomeWorkflowRequest,
)

phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70],
        "sample_b": [1.10, 0.80],
        "sample_c": [0.95, 0.75],
    },
    index=["TSC2;S939;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "_______MSGRPRTTSFAESCKPVQQPSAFG",
        ],
        "protein_id": ["TSC2", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferencePreset.AUTO,
    )
)
pred_mat = kinase_result.prediction_result.pred_mat

# Optional signalome step (requires site_metadata.protein_id)
signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

## Import Contract

`phospy.api` is the canonical namespace where public API types are defined and
organized in source.

Both namespaces are public, with different roles:

- top-level `phospy` is a curated convenience surface for only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`
- requests, configs, results, enums/references, and errors are imported from
  `phospy.api`

## Where To Go Next

- Release framing for current contract: [Release Notes 1.5.0](docs/release_notes/1.5.0.md)
- Guided onboarding: [Quickstart: first workflow](docs/getting-started/quickstart-first-workflow.md)
- Full contract details: [API Guide](docs/api.md), [Validation Guide](docs/validation.md)
- Runnable demos (after you understand the first-run flow):
  - `python examples/dataset_builder_demo.py`
  - `python examples/kinase_workflow_demo.py`
  - `python examples/signalome_workflow_demo.py`

# PhosPy

PhosPy is a Python package for phosphoproteomics workflows with a strict public
contract around one dataset boundary and two workflow entrypoints
(`KinaseWorkflow`, `SignalomeWorkflow`).

PhosPy does **not** expose HTTP endpoints or a web service. The supported user
interfaces are:

- the Python API
- the `phospy` command-line interface

## Install

PhosPy requires Python 3.10 or newer.

For normal package use:

```bash
pip install phospy
```

If you plan to read or write `.parquet` inputs or outputs, install the optional
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

Before you start, make sure your input tables match this minimum shape:

- `phospho`: numeric site-by-sample matrix with site IDs as the row index
- `site_metadata`: rows aligned to `phospho.index`, with `gene_symbol` and `site`
- site IDs should look like `TSC2;S939;`
- add `protein_id` if you plan to run `SignalomeWorkflow`

For human or mouse lanes, pass an explicit `ReferenceBundle` instead of `AUTO`.

### Minimal Python example

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

If you copy the example as-is, you should end up with:

- a strict `AnalysisReadyPhosphoDataset`
- `dataset.phospho.shape == (2, 3)`
- a non-empty kinase prediction matrix in `pred_mat`
- a signalome result only when `protein_id` is present for every interpreted site

### Minimal CLI example

```bash
phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --outdir ./out
```

That command writes a dataset directory, a kinase directory, and a short summary
of written file paths.

## Import Contract

`phospy.api` is the canonical namespace where public API types are defined and organised in source.

Both namespaces are public, with different roles:

- top-level `phospy` is a curated convenience surface for only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`
- requests, configs, results, enums or references, and errors are imported from
  `phospy.api`

## CLI vs Python API

The `phospy` CLI is intentionally narrow and file-first. It supports the public
happy-path flow (`dataset-build`, `kinase`, `signalome`) with selected
high-value runtime knobs.

Use the CLI for reproducible command-line execution from files. Use the Python
API (`phospy.api`) when you need the full request and config surface, including
`DatasetPreprocessingConfig`, explicit `ReferenceBundle` injection, DataFrame
inputs, or advanced scoring and signalome configuration.

## Citation

If you use PhosPy in scientific work, cite this software release using
[`CITATION.cff`](CITATION.cff) and also cite the upstream PhosR project and
publications described in [`NOTICE.md`](NOTICE.md).

## Where To Go Next

- Guided onboarding: [Quickstart: first workflow](docs/getting-started/quickstart-first-workflow.md)
- First failure recovery: [Troubleshooting: first-run and supported-lane failures](docs/getting-started/troubleshooting-first-run.md)
- CLI scope and command usage: [CLI Guide](docs/cli.md)
- Full contract details: [API Guide](docs/api.md), [Validation Guide](docs/validation.md)
- Release framing for the shipped 1.5.0 contract: [Release Notes 1.5.0](docs/release_notes/1.5.0.md)
- Runnable demos after you understand the first-run flow:
  - `python examples/dataset_builder_demo.py`
  - `python examples/kinase_workflow_demo.py`
  - `python examples/signalome_workflow_demo.py`

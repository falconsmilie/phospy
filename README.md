# PhosPy

PhosPy is a Python package for phosphoproteomics workflows.

It supports one clear public workflow:

1. build an analysis-ready dataset
2. run kinase scoring and prediction
3. optionally run signalome analysis

PhosPy does **not** expose HTTP endpoints or a web service. The supported user
interfaces are:

- the Python API
- the `phospy` command-line interface

## Install

PhosPy requires Python 3.10 or newer.

```bash
pip install phospy
```

If you need `.parquet` input or output support:

```bash
pip install "phospy[parquet]"
```

For local development from a clone:

```bash
pip install -e ".[dev]"
pip install -e ".[dev,parquet]"  # optional parquet support
```

## First Run

The recommended beginner lane is deliberately small:

1. build a dataset with `organism=Organism.RAT`
2. run kinase with `references=ReferencePreset.AUTO`
3. run signalome only when `site_metadata.protein_id` is present

Bundled runtime references in this release are rat-only. For human or mouse
work, pass an explicit `ReferenceBundle` in Python instead of `AUTO`.

### Minimum input shape

- `phospho`: numeric site-by-sample matrix
- `site_metadata`: rows aligned to `phospho.index`
- required `site_metadata` columns: `gene_symbol`, `site`
- site IDs should look like `TSC2;S939;`
- add `protein_id` if you plan to run signalome

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
        activity_config=None,  # keep the 2-site first-run example in the supported lane
    )
)

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
```

If you copy the example as-is, you should get:

- a strict `AnalysisReadyPhosphoDataset`
- `dataset.phospho.shape == (2, 3)`
- a non-empty kinase prediction matrix
- a signalome result only when `protein_id` is present for every interpreted site
- compact module summaries via `signalome_result.signalome_modules.table`
- optional provenance sidecars via `signalome_result.site_membership` and `signalome_result.protein_site_context`

### Minimal CLI example

```bash
phospy kinase   --phospho ./input/phospho.csv   --site-metadata ./input/site_metadata.csv   --organism rat   --reference auto   --outdir ./out
```

That command writes a dataset directory, a kinase directory, and a short summary
of written file paths.

## Import Contract

`phospy.api` is the canonical namespace where public API types are defined and organised in source.

Both namespaces are public, with different roles:

- top-level `phospy` is a curated convenience surface for only:
  `AnalysisReadyDatasetBuilder`, `AnalysisReadyPhosphoDataset`,
  `KinaseWorkflow`, `SignalomeWorkflow`
- requests, configs, results, enums, references, and errors are imported from
  `phospy.api`

## CLI vs Python API

Use the CLI for file-based runs in the supported public lane.

Use the Python API (`phospy.api`) when you need:

- DataFrame inputs
- explicit `ReferenceBundle` injection
- dataset preprocessing control
- advanced scoring or signalome configuration

## Documentation

Read these in order if you are new:

1. [Quickstart](docs/getting-started/quickstart-first-workflow.md)
2. [Troubleshooting](docs/getting-started/troubleshooting-first-run.md)
3. [CLI Guide](docs/cli.md) or [API Guide](docs/api.md)

## Citation

If you use PhosPy in scientific work, cite this software release using
[`CITATION.cff`](CITATION.cff) and also cite the upstream PhosR project and
publications described in [`NOTICE.md`](NOTICE.md).

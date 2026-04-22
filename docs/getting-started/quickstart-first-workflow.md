# Quickstart: First Workflow

This quickstart is the supported first-run lane:

1. install the package
2. build an analysis-ready dataset
3. run kinase with bundled references via `ReferencePreset.AUTO`
4. optionally run signalome

## 1. Install

Normal package install:

```bash
pip install phospy
```

If you need parquet file support (`.parquet` input/output), install the optional
parquet extra:

```bash
pip install "phospy[parquet]"
```

If you are developing from a local clone instead:

```bash
pip install -e ".[dev]"
```

For editable installs with parquet support:

```bash
pip install -e ".[dev,parquet]"
```

## 2. Know Required Data and Reference Scope

You need:

- `phospho`: numeric site-by-sample matrix with canonical site IDs in the index
- `site_metadata`: row-aligned to `phospho.index`, with `gene_symbol` and `site`
- `protein_id` in `site_metadata` if you plan to run signalome

Reference behavior:

- `ReferencePreset.AUTO` resolves bundled references from `dataset.organism`.
- Bundled runtime references in this release are rat-only.
- For human/mouse lanes, provide an explicit `ReferenceBundle`.

## 3. Build Dataset + Run Kinase (Recommended First Run)

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, KinaseWorkflow
from phospy.api import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
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
print(pred_mat.round(4))
```

## 4. Optional: Run Signalome

Signalome requires `dataset.site_metadata.protein_id` with non-empty values for
all sites.

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
print(signalome_result.module_assignments.table.head())
```

## 5. Then Use Examples and Deeper Docs

Runnable scripts:

- `python examples/dataset_builder_demo.py`
- `python examples/kinase_workflow_demo.py`
- `python examples/signalome_workflow_demo.py`

Learn next:

- Concepts: [Core concepts](../concepts/core-concepts.md)
- Practical usage: [Tutorials and user guides](../user-guides/index.md)
- Contract details: [API Guide](../api.md), [Validation Guide](../validation.md)

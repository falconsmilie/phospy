# Quickstart: First Workflow

This is the recommended first-run lane:

1. install the package
2. build an analysis-ready, missing-value-free dataset
3. run kinase with bundled references via `ReferencePreset.AUTO`
4. optionally run signalome

## 1. Install

PhosPy requires Python 3.10 or newer.

```bash
pip install phospy
```

If you need `.parquet` input or output:

```bash
pip install "phospy[parquet]"
```

## 2. Prepare two input tables

For a simple first run, you need:

- `phospho`: numeric site-by-sample matrix
- `site_metadata`: one row per phosphosite, aligned to `phospho.index`

Required `site_metadata` columns:

- `gene_symbol`
- `site`

Useful checks:

- row IDs should look like `TSC2;S939;`
- `site_metadata.index` should exactly match `phospho.index`
- add `protein_id` only if you plan to run signalome

Reference rules for this quickstart:

- `ReferencePreset.AUTO` resolves from `dataset.organism`
- bundled runtime references are rat-only in this release
- for human or mouse work, provide an explicit `ReferenceBundle`

## 3. Copy-paste Python example

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
        activity_config=None,  # keep the 2-site first-run example in the supported lane
    )
)

pred_mat = kinase_result.prediction_result.pred_mat
print(pred_mat.round(4))
```

What success looks like:

- `dataset.phospho.shape == (2, 3)`
- `dataset.organism.value == "rat"`
- `pred_mat` is present and non-empty

## 4. Optional: run signalome

Signalome requires explicit, non-empty `dataset.site_metadata.protein_id` for
all interpreted sites.

```python
from phospy import SignalomeWorkflow
from phospy.api import SignalomeWorkflowRequest

signalome_result = SignalomeWorkflow().run(
    SignalomeWorkflowRequest(kinase_result=kinase_result)
)
print(signalome_result.module_assignments.table.head())
```

What success looks like:

- `signalome_result.module_assignments.table` is non-empty
- `signalome_result.signalome_modules.table` is non-empty
- `signalome_result.expanded_signalome` is populated in the supported signalome lane

## 5. CLI version of the same lane

```bash
phospy kinase   --phospho ./input/phospho.csv   --site-metadata ./input/site_metadata.csv   --organism rat   --reference auto   --outdir ./out
```

The command prints a short summary of the files it wrote under `./out`.

## 6. If it fails

Use [Troubleshooting: first run](troubleshooting-first-run.md). It covers the
most common beginner mistakes first.

## 7. Next pages

- [CLI Guide](../cli.md)
- [API Guide](../api.md)
- [Validation Guide](../validation.md)
- `python examples/dataset_builder_demo.py`
- `python examples/kinase_workflow_demo.py`
- `python examples/signalome_workflow_demo.py`

# Quickstart: First Workflow

This quickstart is the supported first-run lane:

1. install the package
2. build an analysis-ready, missing-value-free dataset
3. run kinase with bundled references via `ReferencePreset.AUTO`
4. optionally run signalome

## 1. Install

PhosPy requires Python 3.10 or newer.

Normal package install:

```bash
pip install phospy
```

If you need parquet file support (`.parquet` input or output), install the optional
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

## 2. Know the Minimum Data You Need

For the smallest supported Python example, prepare these two tables:

- `phospho`: numeric site-by-sample matrix
- `site_metadata`: one row per phosphosite, aligned to `phospho.index`

Required `site_metadata` columns:

- `gene_symbol`
- `site`

Useful rule of thumb:

- row IDs should look like `TSC2;S939;`
- `site_metadata.index` should line up exactly with `phospho.index`
- add `protein_id` only if you want to run signalome

Signalome protein-identity prerequisite:

- signalome requires explicit, non-empty `site_metadata.protein_id`
- gene-symbol site-ID prefixes (for example `"<gene_symbol>;<site>;"`) are not a
  substitute for protein identity
- builder flexibility at ingestion does not weaken this downstream contract

Reference behaviour:

- `ReferencePreset.AUTO` resolves bundled references from `dataset.organism`
- bundled runtime references in this release are rat-only
- for human or mouse lanes, provide an explicit `ReferenceBundle`

## 3. Copy-Paste Python Example

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

What success looks like:

- `dataset.phospho.shape == (2, 3)`
- `dataset.organism.value == "rat"`
- `pred_mat` is present and non-empty

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

What success looks like here:

- `signalome_result.module_assignments.table` is non-empty
- `signalome_result.signalome_modules.table` is non-empty
- `signalome_result.expanded_signalome` is populated in the supported signalome lane

## 5. CLI Version of the Same Happy Path

If you prefer files instead of in-memory DataFrames, the matching CLI lane is:

```bash
phospy kinase \
  --phospho ./input/phospho.csv \
  --site-metadata ./input/site_metadata.csv \
  --organism rat \
  --reference auto \
  --outdir ./out
```

The command prints a short summary of the files it wrote under `./out`.

## 6. If It Fails Early

Use [Troubleshooting: first-run and supported-lane failures](troubleshooting-first-run.md) before reading the full validation contract. It is organised by the error you saw and covers the most common supported-lane mistakes.

## 7. Then Use Examples and Deeper Docs

Runnable scripts:

- `python examples/dataset_builder_demo.py`
- `python examples/kinase_workflow_demo.py`
- `python examples/signalome_workflow_demo.py`

Learn next:

- Concepts: [Core concepts](../concepts/core-concepts.md)
- Practical usage: [Tutorials and user guides](../user-guides/index.md)
- First-run recovery: [Troubleshooting: first-run and supported-lane failures](troubleshooting-first-run.md)
- Contract details: [API Guide](../api.md), [Validation Guide](../validation.md)

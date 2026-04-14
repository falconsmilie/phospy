# PhosPy

PhosPy is a focused Python library for a small, supported set of phosphoproteomics workflows inspired by `PhosR`.

The supported public entry points are:

- preprocessing: `PhosphoDataset`
- kinase scoring and prediction: `SimpleKinaseWorkflow`
- signalome analysis: `SignalomeWorkflow`

PhosPy is intentionally narrow. It is not a full `PhosR` replacement.

## Install

PhosPy supports Python 3.10 and newer.

```bash
pip install phospy
```

For parquet output:

```bash
pip install "phospy[parquet]"
```

## Choose a Starting Point

### 1) Preprocessing with `PhosphoDataset`

```python
from phospy.datasets import PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig, SiteMatrixPolicy

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(
    config=CorePreprocessingConfig(
        max_unmatched_fraction=0.1,
        site_matrix_policy=SiteMatrixPolicy(
            missing_data_policy="require_min_observed_values",
            minimum_observed_values=2,
        ),
    )
)

analysis_ready = dataset.preprocessing.run_analysis_ready(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)
```

### 2) Kinase Scoring with `SimpleKinaseWorkflow`

```python
from phospy.api import PredictionRunConfig, SimpleKinaseWorkflow

result = SimpleKinaseWorkflow().run(
    phospho="examples/data/simple_workflow/phospho.tsv",
    total="examples/data/simple_workflow/total.tsv",
    species="rat",
    prediction_config=PredictionRunConfig(
        min_substrates=1,
        min_motif_size=1,
        ensemble_size=2,
        top=3,
        inclusion=2,
        n_iterations=2,
        random_state=7,
    ),
)

pred_mat = result.pred_mat_result.to_frame()
weighted_activity = result.kinase_activity_result.weighted_activity
```

This path handles preprocessing, analysis-ready adaptation, bundled reference selection, prediction, and activity analysis.

### 3) Signalome Analysis with `SignalomeWorkflow`

```python
from phospy.api import SignalomeRunConfig, SignalomeWorkflow

signalome = SignalomeWorkflow().run_from_analysis_ready(
    dataset=result.analysis_ready_dataset,
    scoring_result=result.scoring_result,
    prediction_result=result.prediction_result,
    kinases_of_interest=list(result.pred_mat_result.kinase_names[:2]),
    config=SignalomeRunConfig(signalome_cutoff=0.5),
)
```

Runnable demos:

- [`examples/simple_workflow_demo.py`](examples/simple_workflow_demo.py)
- [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py)
- [`examples/signalome_map_demo.py`](examples/signalome_map_demo.py)

## File Inputs

PhosPy works with:

- total input as TSV
- phospho input as TSV
- `predMat` as CSV, with the first column used as the phosphosite index

For required columns and common validation rules, see [`docs/validation.md`](docs/validation.md).

## CLI

PhosPy ships with a small CLI for file-based preprocessing and optional kinase activity analysis from an existing `predMat`.

```bash
phospy \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --pred-mat examples/data/predMat.csv \
  --phospho-encoding utf-16le \
  --max-unmatched-fraction 0.1 \
  --outdir examples/output
```

## Read Next

- [`docs/api.md`](docs/api.md) for the supported Python API
- [`docs/validation.md`](docs/validation.md) for input rules and common failures
- [`docs/parity.md`](docs/parity.md) for parity scope and `svm_mode`
- [`docs/fixtures.md`](docs/fixtures.md) for fixture and trace rebuild commands
- [`docs/architecture/package-layout.md`](docs/architecture/package-layout.md) for the contributor-facing package layout

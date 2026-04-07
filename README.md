# PhosPy

PhosPy is a focused Python library for selected phosphoproteomics workflows inspired by `PhosR`.

It is designed for users who want to:

- preprocess total and phospho tables
- analyse kinase activity from an existing `predMat`
- run the native Python kinase workflow

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

The examples below use repository paths such as `examples/data/...`. If you installed from PyPI, use your own local file paths.

## Choose the Right Entry Point

### `PhosphoDataset`

Use `PhosphoDataset` when you want validated total and phospho inputs plus the standard preprocessing flow.

```python
from phospy import PhosphoDataset
from phospy.writers import CoreOutputWriter

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

CoreOutputWriter().write(core, outdir="examples/output", format="csv")

site_matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

### `KinaseActivityAnalyzer`

Use `KinaseActivityAnalyzer` when you already have a phosphosite matrix and a `predMat`.

```python
from phospy import KinaseActivityAnalyzer, PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

analyzer = KinaseActivityAnalyzer()
result = analyzer.run(
    pred_mat=analyzer.load_pred_mat("examples/data/predMat.csv"),
    phospho_matrix=core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)

ksea_scores = result.ksea_scores
```

The bundled example data is tiny, so it uses `min_substrates=1` and `top_n_substrates=1`.

### `PhosRPipeline`

Use `PhosRPipeline` when you want file loading, preprocessing, optional kinase analysis, and output publishing in one call.

```python
from phospy import PhosRPipeline

pipeline = PhosRPipeline.from_files(
    total_path="examples/data/total.tsv",
    phospho_path="examples/data/phospho.tsv",
    pred_mat_path="examples/data/predMat.csv",
    phospho_encoding="utf-16le",
    max_unmatched_fraction=0.1,
    kinase_activity_threshold=0.6,
    kinase_activity_min_substrates=1,
    kinase_activity_top_n_substrates=1,
)
outputs = pipeline.run(outdir="examples/output")
```

When `outdir` is set, the pipeline writes the core outputs, any kinase-analysis outputs, and `run_manifest.json`.

### `KinaseWorkflow`

Use `KinaseWorkflow` for the native Python end-to-end scoring and prediction workflow.

A runnable example lives in [`examples/native_workflow_demo.py`](examples/native_workflow_demo.py).

From a repository checkout:

```bash
make native-workflow-demo
```

## File Inputs

PhosPy works with:

- total input as TSV
- phospho input as TSV
- `predMat` as CSV, with the first column used as the phosphosite index

For the default table schema and method-level validation rules, see [`docs/api.md`](docs/api.md).

## CLI

PhosPy also ships with a small CLI for the file-based preprocessing path and optional `predMat` analysis.

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

- [`docs/api.md`](docs/api.md) for the public Python API and CLI options
- [`docs/validation.md`](docs/validation.md) for the validation checklist
- [`docs/parity.md`](docs/parity.md) for the PhosR parity scope
- [`docs/fixtures.md`](docs/fixtures.md) for fixture and trace directories
- [`CONTRIBUTING.md`](.github/CONTRIBUTING.md) for local development
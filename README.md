# PhosPy

`PhosPy` is a small Python library for selected phosphoproteomics workflows inspired by the R `PhosR` package.

Use it when you want to:

- preprocess total and phospho tables
- analyse kinase activity from a `predMat`
- run the native Python kinase workflow

It is intentionally narrow. It does not aim to reproduce all of `PhosR`.

## Install

PhosPy supports Python 3.10 and newer.

```bash
pip install phospy
```

The file paths below use `examples/data/...`, so they assume a repository checkout. If you installed from PyPI, use your own input-file paths instead.

## Pick the Method You Need

### 1. Core preprocessing from total and phospho tables

Use `PhosphoDataset` when you want validated inputs and the standard preprocessing flow. `PhosphoDataset` owns mutable pandas tables and isolates itself from caller-managed inputs at construction time. Prefer `dataset.total_df_copy`, `dataset.phospho_df_copy`, and `dataset.copy_inputs()` for inspection, reporting, and other read-oriented work. Use `dataset.total_df_live` and `dataset.phospho_df_live` only when you intentionally want shared workspace state.

Result bundles returned by the public analyzers, pipeline, and workflow are detached snapshot-style outputs. Mutating those result tables affects only the returned result bundle and does not feed changes back into dataset workspace state.

```python
from phospy import PhosphoDataset
from phospy.writers import CoreOutputWriter

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

writer = CoreOutputWriter()
writer.write(core, outdir="examples/output", format="csv")

site_matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

`dataset.preprocessing.run(...)` returns a `CoreProcessingResult` with:

- `total_unique`
- `total_filtered`
- `phospho_filtered`
- `phospho_corrected`
- `site_matrix`

Use `CoreOutputWriter` from `phospy.writers` when you want to write the core outputs to disk.

### 2. Kinase activity analysis from an existing `predMat`

Use `KinaseActivityAnalyzer` when you already have a phosphosite matrix and a `predMat`.

```python
from phospy import KinaseActivityAnalyzer, PhosphoDataset
from phospy.io import load_pred_mat

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

analyzer = KinaseActivityAnalyzer()
result = analyzer.analyze(
    pred_mat=load_pred_mat("examples/data/predMat.csv"),
    phospho_matrix=core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)
analyzer.write_outputs(result, outdir="examples/output")

ksea_scores = result.ksea_scores
target_counts = result.target_counts
```

The bundled example uses `min_substrates=1` and `top_n_substrates=1` because the example data is very small.

### 3. One-shot pipeline from files

Use `PhosRPipeline` when you want file loading, preprocessing, optional kinase analysis, and output writing in one call.

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

This writes the core CSV outputs and, when `pred_mat_path` is provided, the downstream kinase-analysis tables using the explicit `kinase_activity_*` settings carried on the pipeline request.
A pipeline run also writes `run_manifest.json`.

### 4. Native Python kinase workflow

Use `KinaseWorkflow` for the native end-to-end prediction workflow.

A complete runnable example lives in [`examples/native_workflow_demo.py`](examples/native_workflow_demo.py).

From a repository checkout, run:

```bash
make native-workflow-demo
```

## Input Files

For file-based workflows:

- total input is read as TSV
- phospho input is read as TSV
- `predMat` is read as CSV, using the first column as the phosphosite index

For the default schema, the expected columns are documented in [`docs/api.md`](docs/api.md).

## CLI

After installation, you can run the CLI on your own files.

```bash
phospy \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --pred-mat examples/data/predMat.csv \
  --phospho-encoding utf-16le \
  --max-unmatched-fraction 0.1 \
  --outdir examples/output
```

The CLI covers the file-based preprocessing path and optional `predMat` analysis.

## Where to Read Next

- [`docs/api.md`](docs/api.md) for method signatures, parameters, validation, and examples
- [`docs/validation.md`](docs/validation.md) for the validation quick guide
- [`docs/parity.md`](docs/parity.md) for parity to the R `PhosR` package
- [`docs/fixtures.md`](docs/fixtures.md) for fixture and trace layout
- [`CONTRIBUTING.md`](CONTRIBUTING.md) for local development and test commands

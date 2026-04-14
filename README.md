# PhosPy

PhosPy is a focused Python library for a small, supported set of phosphoproteomics workflows inspired by `PhosR`.

It is built for a narrow set of jobs:

- preprocess total and phospho tables
- analyse kinase activity from an existing `predMat`
- generate a `predMat` from phosphosite inputs
- run the native Python kinase workflow
- build downstream signalome outputs

PhosPy is intentionally narrow. It is **not** a full `PhosR` replacement.

## Install

PhosPy supports Python 3.10 and newer.

```bash
pip install phospy
```

For parquet output:

```bash
pip install "phospy[parquet]"
```

Repository examples use paths such as `examples/data/...`. If you installed from PyPI, point those examples at your own files instead.

## Choose the Right Starting Point

### `SimpleKinaseWorkflow`

Use this when you already have a phospho table, know the species, and want the shortest supported end-to-end path.

```python
from phospy.api import SimpleKinaseWorkflow

result = SimpleKinaseWorkflow().run(
    phospho="study_phospho.tsv",
    total="study_total.tsv",
    species="rat",
    reference="auto",
)

pred_mat = result.pred_mat_result.to_frame()
weighted_activity = result.kinase_activity_result.weighted_activity
```

This path handles:

- preprocessing
- the analysis-ready phosphosite boundary
- bundled reference resolution
- `predMat` generation
- kinase activity analysis

Bundled reference support is intentionally small:

- supported species: `rat`
- supported references: `auto`, `l6`, `l6_native`
- `auto` currently resolves to `l6_native`

Runnable example: [`examples/simple_workflow_demo.py`](examples/simple_workflow_demo.py)

### `PhosphoDataset`

Use this when you want validated total and phospho inputs plus the standard preprocessing flow.

```python
from phospy.datasets import PhosphoDataset
from phospy.io.writers import CoreOutputWriter
from phospy.preprocessing import CorePreprocessingConfig

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)

CoreOutputWriter().write(core, outdir="examples/output", format="csv")

site_matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

### `KinaseActivityAnalyzer`

Use this when you already have a phosphosite matrix and a `predMat`.

```python
from phospy.activities import KinaseActivityAnalyzer
from phospy.datasets import PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)

analyzer = KinaseActivityAnalyzer()
result = analyzer.run(
    pred_mat=analyzer.load_pred_mat("examples/data/predMat.csv"),
    phospho_matrix=core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)
```

The bundled example data is tiny, so it uses `min_substrates=1` and `top_n_substrates=1`.

### `PhosRPipeline`

Use this when you want file loading, preprocessing, optional kinase analysis, and output publishing in one call.

```python
from phospy.api import DatasetLoadOptions, KinaseActivityConfig
from phospy.pipeline import PhosRPipeline
from phospy.preprocessing import CorePreprocessingConfig

pipeline = PhosRPipeline.from_files(
    total_path="examples/data/total.tsv",
    phospho_path="examples/data/phospho.tsv",
    pred_mat_path="examples/data/predMat.csv",
    dataset_options=DatasetLoadOptions(phospho_encoding="utf-16le"),
    preprocessing_config=CorePreprocessingConfig(max_unmatched_fraction=0.1),
    activity_config=KinaseActivityConfig(
        threshold=0.6,
        min_substrates=1,
        top_n_substrates=1,
    ),
)
outputs = pipeline.run(outdir="examples/output")
```

When `outdir` is set, the pipeline writes core outputs, optional kinase-analysis outputs, and `run_manifest.json`.

### `PredMatWorkflow`, `KinaseWorkflow`, and `SignalomeWorkflow`

Use these when you need direct control over workflow-shaped inputs such as:

- `site_sequences`
- `substrate_map`
- `motif_sequences`
- `ReferenceBundle`
- intermediate scoring and prediction outputs
- signalome construction

Runnable examples:

- [`examples/native_workflow_demo.py`](examples/native_workflow_demo.py)
- [`examples/predmat_workflow_demo.py`](examples/predmat_workflow_demo.py)
- [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py)

## File Inputs

PhosPy works with:

- total input as TSV
- phospho input as TSV
- `predMat` as CSV, with the first column used as the phosphosite index

For required columns and common validation rules, see [`docs/validation.md`](docs/validation.md).

## CLI

PhosPy also ships with a small CLI for file-based preprocessing and optional `predMat` analysis.

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

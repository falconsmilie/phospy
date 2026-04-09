# PhosPy

PhosPy is a focused Python library for selected phosphoproteomics workflows inspired by `PhosR`.

It is designed for users who want to:

- preprocess total and phospho tables
- analyse kinase activity from an existing `predMat`
- generate a `predMat` from phosphosite inputs
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

The examples below use repository paths such as `examples/data/...`. If you installed from PyPI, use your own local file
paths.

## Choose the Right Entry Point

### `PhosphoDataset`

Use `PhosphoDataset` when you want validated total and phospho inputs plus the standard preprocessing flow.

```python
from pathlib import Path

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

Use `PhosRPipeline` when you want file loading, preprocessing, optional kinase analysis, and output publishing in one
call.

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

### `PredMatWorkflow`

Use `PredMatWorkflow` when your goal is to generate a `predMat` from phosphosite inputs and export it as CSV.

A runnable end-to-end example lives in [`examples/predmat_workflow_demo.py`](examples/predmat_workflow_demo.py).

```python
import json
from pathlib import Path

import pandas as pd

from phospy import PredMatWorkflow

phospho_matrix = pd.read_csv("predmat_phospho_matrix.csv", index_col=0)
site_sequences = json.loads(Path("predmat_site_sequences.json").read_text())
substrate_map = json.loads(Path("predmat_substrate_map.json").read_text())
motif_sequences = json.loads(Path("predmat_motif_sequences.json").read_text())

workflow = PredMatWorkflow(flank_size=2, svm_mode="default")
result = workflow.run(
    phospho_matrix=phospho_matrix,
    substrate_map=substrate_map,
    site_sequences=site_sequences,
    motif_sequences=motif_sequences,
    min_substrates=2,
    min_motif_size=2,
    ensemble_size=3,
    top=4,
    score_threshold=0.75,
    inclusion=3,
    n_iterations=2,
    random_state=17,
)

pred_mat = result.pred_mat_result.to_frame(copy=False)
result.pred_mat_result.to_csv("predMat.csv")
```

Use `svm_mode="default"` for the recommended stable native path. Use `svm_mode="r_parity"` when you want the supported parity-oriented learner, sampling, and final-scoring preset for parity-sensitive comparisons.

When prediction thresholds are too strict and no kinase candidates qualify, PhosPy raises `NoCandidateKinasesError` instead of returning an invalid empty `predMat`.

### `KinaseWorkflow`

Use `KinaseWorkflow` for the fuller native Python scoring and prediction workflow when you want the intermediate
profile, motif, and combined scoring outputs as well as the final prediction matrix.

A runnable example lives in [`examples/native_workflow_demo.py`](examples/native_workflow_demo.py).

From a repository checkout:

```bash
make native-workflow-demo
```

Use `svm_mode="default"` for the recommended stable native path.
Use `svm_mode="r_parity"` when you want the supported parity-oriented learner, sampling, and final-scoring preset.

### `SignalomeWorkflow`

Use `SignalomeWorkflow` when you already have scoring and prediction outputs and want the next downstream layer:
signalome modules plus map-ready and network-ready derived outputs.

A runnable end-to-end example lives in [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py).

```python
from phospy import PredMatWorkflow, SignalomeWorkflow

pred_mat_result = PredMatWorkflow(flank_size=2, svm_mode="default").run(...)
signalome_result = SignalomeWorkflow().run(
    scoring_result=pred_mat_result.scoring_result,
    prediction_result=pred_mat_result.prediction_result,
    expression_matrix=phospho_matrix,
    kinases_of_interest=["KINASE_A", "KINASE_B"],
    signalome_cutoff=0.5,
)

map_data = signalome_result.to_map_data()
network_data = signalome_result.to_network_data()
```

Use `signalome_result.to_csv(...)`, `map_data.to_csv(...)`, and `network_data.to_csv(...)` when you want exportable
tables. The same `PredMatWorkflow` call can use `svm_mode="r_parity"` when you want the parity-oriented prediction preset before constructing downstream signalome outputs.

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
- [`docs/parity.md`](docs/parity.md) for the PhosR parity scope, release thresholds, and prediction-mode intent
- [`docs/adr/0002-r-parity-public-preset.md`](docs/adr/0002-r-parity-public-preset.md) for the explicit public support decision on `r_parity`
- [`docs/fixtures.md`](docs/fixtures.md) for fixture and trace directories
- [`CONTRIBUTING.md`](.github/CONTRIBUTING.md) for local development
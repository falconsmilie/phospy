# PhosPy

PhosPy is a Python package for a small, supported set of phosphoproteomics workflows inspired by `PhosR`.

It is built for three main jobs:

- **Preprocessing** raw total and phospho tables with `PhosphoDataset`
- **Kinase scoring and prediction** with `SimpleKinaseWorkflow`
- **Signalome analysis** with `SignalomeWorkflow`

PhosPy is intentionally narrow. It is not a full `PhosR` replacement.

## Install

PhosPy supports Python 3.10 and newer.

```bash
pip install phospy
```

Optional parquet support:

```bash
pip install "phospy[parquet]"
```

## Start Here

### Use `PhosphoDataset` when you want preprocessing only

```python
from phospy.datasets import PhosphoDataset
from phospy.preprocessing import CorePreprocessingConfig

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)

analysis_ready = dataset.preprocessing.run_analysis_ready(
    config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
)
```

You get a validated, analysis-ready dataset that you can pass into downstream workflows.

### Use `SimpleKinaseWorkflow` for the common end-to-end path

```python
from phospy.api import PredictionRunConfig, SimpleKinaseWorkflow

with SimpleKinaseWorkflow().run(
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
) as result:
    pred_mat = result.pred_mat_result.to_owned_frame()
    weighted_activity = result.kinase_activity_result.weighted_activity
```

This path handles preprocessing, analysis-ready adaptation, bundled reference selection, prediction, and kinase activity analysis.

### Workflow result access patterns

`SimpleKinaseWorkflow.run(...)` returns `SimpleKinaseWorkflowResult`.

- Use `result.pred_mat_result` for the canonical predMat table contract.
- Use `result.prediction_result` when you need full prediction payload details (`pred_matrix`, `substrate_list`, optional traces).
- Use `result.scoring_result` for `profile_scores`, `combined_scores`, and `weights`.
- Use `result.kinase_activity_result` for activity summaries such as `weighted_activity` and `ksea_scores`.

`predMat` is already part of `prediction_result` (`prediction_result.pred_mat_result`), so there is no separate predMat workflow.

Result ownership naming is explicit:
- `to_<resource>()` returns a detached safe copy
- `to_owned_<resource>()` returns cheap shared owned state
- `to_mutable_<resource>_unsafe()` returns explicit mutable shared state

`phospy.internal.kinase_workflows.KinaseWorkflow` is an internal orchestration helper; the supported public lane is `phospy.api.SimpleKinaseWorkflow`.

### Use `SignalomeWorkflow` when you already have workflow outputs

The next example assumes you already have a `result` from `SimpleKinaseWorkflow.run(...)`.

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

`run_from_analysis_ready(...)` expects an `AnalysisReadyPhosphoDataset` (typically `result.analysis_ready_dataset` from `SimpleKinaseWorkflow.run(...)`).

## Minimum Input Expectations

PhosPy works with:

- **total** input as TSV
- **phospho** input as TSV
- **`predMat`** as CSV, with the first column used as the phosphosite index

Default required columns:

| Input | Required columns |
| --- | --- |
| Total table | `genes`, `group1` to `group6` |
| Phospho table | `uid`, `gene_names`, `gene_p_site`, `localization_prob`, `centralized_sequence`, `p_group1` to `p_group6` |

For full validation rules and common failure cases, see [`docs/validation.md`](docs/validation.md).

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

Use `phospy --help` to see every option.

## Runnable Examples

- [`examples/simple_workflow_demo.py`](examples/simple_workflow_demo.py)
- [`examples/signalome_workflow_demo.py`](examples/signalome_workflow_demo.py)
- [`examples/signalome_map_demo.py`](examples/signalome_map_demo.py)
- [`examples/kinase_activity_analyzer_demo.py`](examples/kinase_activity_analyzer_demo.py)

## Contributor Hygiene Check

CI enforces repository hygiene and fails if ignored/generated artefacts are tracked, including:

- `__pycache__/`
- `*.egg-info/`
- `benchmarks/reports/` output

If CI fails on this check, untrack the files before committing (for example with `git rm --cached <path>`).

## Read Next

### Beginner guides

- [`docs/api.md`](docs/api.md) for the supported Python API
- [`docs/validation.md`](docs/validation.md) for input rules and common failures

### Advanced and contributor docs

- [`docs/parity.md`](docs/parity.md) for parity scope and `svm_mode`
- [`docs/fixtures.md`](docs/fixtures.md) for fixture and trace rebuild commands
- [`docs/architecture/package-layout.md`](docs/architecture/package-layout.md) for package ownership and contributor layout rules
- [`docs/architecture/preprocessing-owned-paths.md`](docs/architecture/preprocessing-owned-paths.md) for preprocessing ownership boundaries and copy-churn benchmark guard rules
- [`docs/roadmap.md`](docs/roadmap.md) for likely next areas of work

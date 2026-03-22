# PhosPy

`PhosPy` is an unofficial Python port of selected PhosR workflow components for phosphoproteomics.

It brings a practical subset of the PhosR workflow into Python, including preprocessing utilities, downstream
kinase-analysis summaries, and a native end-to-end kinase scoring and prediction path. The aim is to make these
workflows easier to use in Python while staying clear about attribution, scope, and current parity limits.

## What PhosPy Can Do Today

PhosPy currently supports four main areas of work.

### Core Data Processing

- collapse duplicate total-protein genes by highest mean signal after gene normalisation
- replace sentinel missing values
- apply minimum-observation filtering
- correct phosphosites against total proteome
- add optional pairwise comparison columns from corrected phosphosites
- construct site matrices with duplicate site collapse

### Downstream Kinase Analysis

- calculate weighted kinase activity from `predMat`
- calculate KSEA-style scores and kinase target counts from `predMat`

### Native Kinase Workflow

- build kinase substrate profiles with `KinaseProfileBuilder`
- score phosphosites from substrate profiles with `KinaseScorer`
- score phosphosite motifs with `KinaseMotifScorer`
- combine profile and motif scores with `combine_profile_and_motif_scores()`
- select candidate substrates and run adaptive SVM prediction with `KinasePredictor`
- run the full native workflow through `KinaseWorkflow` or `run_kinase_workflow()`

### Project Tooling

- a minimal CLI
- a class-based public API
- parity tests against R-generated fixtures for selected workflow seams

## Still in Progress

PhosPy is not yet a full Python replacement for PhosR. Important areas still in progress include:

- broader native Python coverage across the PhosR package
- a native Python replacement for `Signalomes()`
- broader numerical parity claims for the newer native kinase workflow beyond the fixture-backed seams documented in [`docs/parity.md`](docs/parity.md)
- full feature, behaviour, and output parity with the original R package

## Install

Base install:

```bash
pip install -e .
```

For the native prediction workflow, install the machine-learning extra:

```bash
pip install -e ".[ml]"
```

For tests:

```bash
pip install -e ".[test]"
pytest
```

For development checks:

```bash
pip install -e ".[test,dev]"
pre-commit install
pre-commit run --all-files
```

## Quick Start

### Core Preprocessing

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.process_core()

site_matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

### Adding Pairwise Comparisons

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "total.tsv",
    "phospho.tsv",
    comparisons=[("group1", "group4"), ("group2", "group5")],
)
core = dataset.process_core()
```

If you do not pass `comparisons=...`, the core pipeline does not add pairwise comparison columns.

### Downstream Kinase Activity From `predMat`

```python
from phospy import KinaseActivityAnalyzer, PhosphoDataset

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.process_core()

analyzer = KinaseActivityAnalyzer.from_csv("predMat.csv")
kinase = analyzer.analyze(core.site_matrix.matrix)

weighted_activity = kinase.weighted_activity
ksea_scores = kinase.ksea_scores
target_counts = kinase.target_counts
```

### End-to-End Pipeline

```python
from phospy import PhosRPipeline

pipeline = PhosRPipeline.from_files(
    total_path="total.tsv",
    phospho_path="phospho.tsv",
    pred_mat_path="predMat.csv",
)
outputs = pipeline.run(outdir="output")

site_matrix = outputs.core.site_matrix.matrix
kinase_activity = outputs.kinase_activity
```

### Native End-to-End Kinase Workflow

This workflow requires the machine-learning extra:

```bash
pip install -e ".[ml]"
```

Then you can run the native workflow in one orchestration call:

```python
from phospy import run_kinase_workflow

result = run_kinase_workflow(
    phospho_matrix=phospho_matrix,
    substrate_map={
        "KINASE_A": ["SITE_1", "SITE_2", "SITE_3", "SITE_4"],
        "KINASE_B": ["SITE_5", "SITE_6", "SITE_7", "SITE_8"],
    },
    site_sequences={
        "SITE_1": "QQAAAAAYY",
        "SITE_2": "QQAAAAAYY",
        "SITE_3": "QQAAAAAYY",
        "SITE_4": "QQAAAAAYY",
        "SITE_5": "QQTTTTTYY",
        "SITE_6": "QQTTTTTYY",
        "SITE_7": "QQTTTTTYY",
        "SITE_8": "QQTTTTTYY",
    },
    motif_sequences={
        "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY", "QQAAAAAYY"],
        "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY", "QQTTTTTYY"],
    },
    min_substrates=2,
    min_motif_size=2,
    ensemble_size=3,
    top=4,
    score_threshold=0.75,
    inclusion=3,
    n_iterations=2,
    random_state=17,
    flank_size=2,
)

profile_matrix = result.profile_result.profile_matrix
motif_scores = result.motif_result.motif_scores
combined_scores = result.scoring_result.combined_scores
pred_matrix = result.prediction_result.pred_matrix
```

For a profile-only fallback path, omit `motif_sequences` and pass `allow_profile_only_fallback=True`.

If you want prediction settings that more closely track the R learner seam, you can also pass `svm_mode="r_parity"`.

### Lower-Level Native Building Blocks

If you need finer control, PhosPy also exposes lower-level building blocks for profile construction, motif scoring,
profile scoring, score combination, and prediction:

```python
from phospy import (
    KinaseMotifScorer,
    KinasePredictor,
    KinaseProfileBuilder,
    KinaseScorer,
    combine_profile_and_motif_scores,
)
```

## CLI Demo

A small synthetic dataset is included.

After installation, you can run:

```bash
phospy \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --pred-mat examples/data/predMat.csv \
  --outdir examples/output
```

This produces:

- `df_total_unique.csv`
- `df_total_filtered.csv`
- `df_phospho_filtered.csv`
- `df_phospho_corrected.csv`
- `phosr_input.csv`
- `mat_phospho_corrected.csv`
- `site_sequences.csv`
- `kinase_activity_matrix.csv`
- `ksea_scores.csv`
- `ksea_counts.csv`
- `kinase_target_counts.csv`
- `kinase_target_table.csv`

## Project Status

PhosPy is currently a structured Python package for PhosR-style preprocessing and downstream kinase-analysis summaries.
It also includes a live R-backed parity harness and a growing native kinase workflow for kinase scoring and prediction.

The newer native workflow should be described as an evolving port, not a parity-complete reimplementation. For the
current parity contract, fixture paths, and maintenance rules, see [`docs/parity.md`](docs/parity.md).

## Development

Code style and linting are enforced with Ruff through `pre-commit`. The local workflow is intentionally small:

- `ruff check --fix` for linting, import sorting, and safe fixes
- `ruff format` for formatting
- `pytest` for unit tests
- `pytest -m parity` for the R-backed parity layer

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for local setup and [`docs/parity.md`](docs/parity.md) for the parity contract.

## Attribution

All scientific credit for the original methods, package design, and biological workflow belongs to the PhosR authors and
maintainers.

Please cite and acknowledge the original PhosR work when using this repository:

- Kim, H. J., Kim, T., Hoffman, N. J., Xiao, D., James, D. E., Humphrey, S. J., & Yang, P. (2021). *PhosR enables
  processing and functional analysis of phosphoproteomic data*. Cell Reports, 34(8), 108771.
- Kim, H., Kim, T., Xiao, D., & Yang, P. (2021). *Protocol for the processing and downstream analysis of
  phosphoproteomic data with PhosR*. STAR Protocols, 2(2), 100585.
- Original R package: `PYangLab/PhosR`

PhosPy should be described as an unofficial port unless and until the original PhosR authors choose to endorse or
participate in it.

## License

This repository is distributed under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [`LICENSE`](LICENSE).

That choice is deliberate. PhosR is distributed under GPL-3, and the GNU GPL FAQ treats translation of a program into
another programming language as a kind of modification or translation under copyright law. This project therefore uses
GPL-3.0-only as the conservative licensing position for a Python port.
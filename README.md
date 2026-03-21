# PhosPy

`PhosPy` is an unofficial Python port of selected PhosR workflow components for phosphoproteomics.

It brings a practical subset of the PhosR workflow into Python, including preprocessing utilities, downstream
kinase-analysis summaries, and a native end-to-end kinase scoring and prediction path. The aim is to make these
workflows easier to use in Python while staying clear about attribution, scope, and current parity limits.

## What PhosPy Can Do Today

PhosPy currently includes:

- collapse duplicate total-protein genes by highest mean signal after gene normalisation
- replace sentinel missing values
- minimum-observation filtering
- phosphosite correction against total proteome
- optional pairwise comparison columns from corrected phosphosites
- site-matrix construction with duplicate site collapse
- weighted kinase activity from `predMat`
- KSEA-style scores and kinase target counts from `predMat`
- a minimal CLI
- a class-based public API
- parity tests against R-generated fixtures for selected workflow seams
- native kinase substrate-profile construction via `KinaseProfileBuilder`
- native profile-based kinase scoring via `KinaseScorer.score_phosphosite_profiles()`
- native motif-frequency scoring via `KinaseMotifScorer` and `score_phosphosite_motifs()`
- native rank-weighted motif/profile score combination via `combine_profile_and_motif_scores()`
- native candidate-substrate selection and adaptive SVM prediction via `KinasePredictor`
- a dedicated end-to-end native orchestration layer via `KinaseWorkflow` and `run_kinase_workflow()`

## Still in Progress

PhosPy is not yet a full Python replacement for PhosR. In particular, the following areas are still in progress:

- full native Python coverage of the broader PhosR package
- a native Python replacement for `Signalomes()`
- broader numerical parity claims for the newer native kinase workflow beyond the fixture-backed seams documented in [
  `docs/parity.md`](docs/parity.md)
- full feature, behaviour, and output parity with the original R package

## Install

Base install:

```bash
pip install -e .
```

For the native prediction path, install the machine-learning extra:

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
```

### Native End-to-End Kinase Workflow

```python
from phospy import run_kinase_workflow

result = run_kinase_workflow(
    phospho_matrix=phospho_matrix,
    substrate_map={
        "PRKACA": ["SITE_1", "SITE_2"],
        "BTK": ["SITE_3", "SITE_4"],
    },
    site_sequences={
        "SITE_1": "QQAAAAAYY",
        "SITE_2": "QQAAAAAYY",
        "SITE_3": "QQTTTTTYY",
        "SITE_4": "QQTTTTTYY",
    },
    motif_sequences={
        "PRKACA": ["QQAAAAAYY", "QQAAAAAYY"],
        "BTK": ["QQTTTTTYY", "QQTTTTTYY"],
    },
    flank_size=2,
    ensemble_size=10,
    random_state=42,
)

profile_matrix = result.profile_result.profile_matrix
motif_scores = result.motif_result.motif_scores
combined_scores = result.scoring_result.combined_scores
pred_matrix = result.prediction_result.pred_matrix
```

For a profile-only fallback path, omit `motif_sequences` and pass `allow_profile_only_fallback=True`.

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

A small synthetic dataset is included:

```bash
PYTHONPATH=src python -m phospy.cli \
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

Today, PhosPy is best described as a structured Python package for PhosR-style preprocessing and downstream
kinase-analysis summaries, with a live R-backed parity harness and a growing native kinase workflow covering profile
construction, motif scoring, profile scoring, rank-weighted score combination, and adaptive SVM prediction.

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
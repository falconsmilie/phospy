# PhosPy

`PhosPy` is an unofficial Python port of selected PhosR workflow components for phosphoproteomics.

It brings a practical subset of the PhosR workflow into Python, including core preprocessing, downstream kinase-analysis
summaries, and a growing native workflow for kinase scoring and prediction.

## What PhosPy Covers

PhosPy currently supports three main areas of work:

- **Core data processing** for PhosR-style preprocessing, phosphosite correction, optional pairwise comparisons, and
  site-matrix construction
- **Downstream kinase analysis** from `predMat`, including weighted kinase activity, KSEA-style scores, and target
  counts
- **Native kinase workflow** for substrate-profile construction, motif scoring, score combination, candidate-substrate
  selection, and adaptive SVM prediction

Alongside that, the repository includes a small CLI, a class-based public API, and fixture-backed parity tests for
selected workflow seams.

## Install

Base install:

```bash
pip install -e .
```

For the native prediction workflow:

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

If your analysis needs pairwise comparisons, pass them explicitly:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "total.tsv",
    "phospho.tsv",
    comparisons=[("group1", "group4"), ("group2", "group5")],
)
core = dataset.process_core()
```

### Downstream Kinase Analysis From `predMat`

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
    top=4,
    score_threshold=0.75,
    inclusion=3,
    random_state=17,
)

pred_matrix = result.prediction_result.pred_matrix
```

For a profile-only fallback path, omit `motif_sequences` and pass `allow_profile_only_fallback=True`.

If you need more control, the lower-level building blocks are also public, including `KinaseProfileBuilder`,
`KinaseMotifScorer`, `KinaseScorer`, and `KinasePredictor`.

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

The example output directory under `examples/output/` shows the generated CSV files.

## Testing and Reference Data

The supporting documentation is split by topic:

- [`docs/parity.md`](docs/parity.md) explains what parity means here, how to run the fixture-backed suite, and which
  options affect pytest output
- [`docs/fixtures.md`](docs/fixtures.md) explains the fixture and trace directories, generation commands, and which
  outputs are committed reference data
- [`CONTRIBUTING.md`](CONTRIBUTING.md) covers local setup, linting, tests, and CI expectations

## Project Status

PhosPy is currently a structured Python package for PhosR-style preprocessing and downstream kinase-analysis summaries.
It also includes a live R-backed parity harness and a growing native kinase workflow for kinase scoring and prediction.

It is not yet a full Python replacement for PhosR. In particular, parity claims for the newer native workflow should
stay limited to the fixture-backed seams documented in [`docs/parity.md`](docs/parity.md).

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

This repository is distributed under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [
`LICENSE`](LICENSE).

That choice is deliberate. PhosR is distributed under GPL-3, and the GNU GPL FAQ treats translation of a program into
another programming language as a kind of modification or translation under copyright law. This project therefore uses
GPL-3.0-only as the conservative licensing position for a Python port.
# PhosPy

`PhosPy` 1.0.0 is an unofficial Python implementation of selected PhosR-style workflows for phosphoproteomics.

It is designed for people who want a small, Python-native way to:

- preprocess phosphoproteomics tables
- analyse kinase activity from `predMat`
- run a native kinase workflow from scoring through prediction

PhosPy is deliberately narrow. It is **not** a full replacement for the R `PhosR` package.

## Install

Install the supported Python API and the `phospy` CLI:

```bash
pip install phospy
```

## What You Can Do With PhosPy

### Preprocess Phosphoproteomics Data

Start from total and phospho input tables and produce corrected phosphosite matrices for downstream use.

### Analyse Kinase Activity From `predMat`

Generate weighted activity scores, KSEA-style summaries, and target counts from predicted kinase–substrate
relationships.

### Run a Native Kinase Workflow

Construct substrate profiles, score motifs, combine evidence, select candidates, and perform adaptive SVM-based kinase
prediction.

## Supported Public API for 1.0.0

The stable root-level API for 1.0.0 is intentionally small:

- `PhosphoDataset`
- `PhosRPipeline`
- `KinaseActivityAnalyzer`
- `KinaseWorkflow`

Returned result dataclasses:

- `CoreProcessingResult`
- `SiteMatrixResult`
- `CoreOutputs`
- `KinaseActivityResult`
- `KinasePredictionResult`
- `KinaseWorkflowResult`

The examples below use only those imports.

## Quick Start

The quickest way to get started is to use the bundled example data in `examples/data/`.

### Core Preprocessing

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
)
core = dataset.process_core(max_unmatched_fraction=0.1)

site_matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

For the bundled example data, `site_matrix.index.tolist()` is `['BTK;Y551;']`.

If your analysis needs explicit pairwise comparisons, pass them when you build the dataset:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    comparisons=[("group1", "group4"), ("group2", "group5")],
)
core = dataset.process_core(max_unmatched_fraction=0.1)
```

### Downstream Kinase Analysis From `predMat`

```python
from phospy import KinaseActivityAnalyzer, PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
)
core = dataset.process_core(max_unmatched_fraction=0.1)

analyzer = KinaseActivityAnalyzer.from_csv("examples/data/predMat.csv")
kinase = analyzer.analyze(
    core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)

target_counts = kinase.target_counts
ksea_scores = kinase.ksea_scores
```

For the bundled example data, `target_counts.to_dict()` is `{'PRKACA': 3, 'BTK': 2}`.

### End-to-End Pipeline

```python
from phospy import PhosRPipeline

pipeline = PhosRPipeline.from_files(
    total_path="examples/data/total.tsv",
    phospho_path="examples/data/phospho.tsv",
    pred_mat_path="examples/data/predMat.csv",
    max_unmatched_fraction=0.1,
)
outputs = pipeline.run(outdir="examples/output")
```

This writes the core CSV outputs together with downstream kinase-analysis tables, including
`kinase_target_table.csv`.

### Native End-to-End Kinase Workflow

A complete runnable native-workflow example is included at
[`examples/native_workflow_demo.py`](examples/native_workflow_demo.py):

```bash
python examples/native_workflow_demo.py
```

That example uses only the supported 1.0.0 root API and prints a small prediction matrix for a synthetic two-kinase
setup.

## Command-Line Demo

After installation, you can run the bundled example from the command line:

```bash
phospy \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --pred-mat examples/data/predMat.csv \
  --max-unmatched-fraction 0.1 \
  --outdir examples/output
```

The example output directory in `examples/output/` shows the generated CSV files.

`--max-unmatched-fraction` defaults to `0.0`. That means protein correction fails if the inner join would silently drop
any phosphosite rows. Raise it only when you want to allow a small, bounded amount of row loss.

## Where to Go Next

If you want more detail, these are the most useful follow-on docs:

- [`docs/validation-and-parity.md`](docs/validation-and-parity.md) explains how validation is approached in PhosPy
- [`docs/parity.md`](docs/parity.md) explains what parity means here, especially for the native kinase workflow
- [`docs/fixtures.md`](docs/fixtures.md) maps the committed fixture and trace directories
- [`docs/roadmap.md`](docs/roadmap.md) outlines the most likely next steps after 1.0.0
- [`CHANGELOG.md`](CHANGELOG.md) contains the 1.0.0 release notes

If you want to contribute or work from a local checkout, see [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Known Limitations

A few boundaries are worth knowing up front:

- **Selective scope only.** PhosPy 1.0.0 covers the workflows documented above and nothing broader.
- **Parity is seam-level, not package-wide.** Validation claims are limited to the committed fixture-backed seams
  described in [`docs/validation-and-parity.md`](docs/validation-and-parity.md) and [`docs/parity.md`](docs/parity.md).
- **`KinaseWorkflow` is native first.** It includes an `svm_mode="r_parity"` option for narrower learner-seam
  comparison, but the default mode is the preferred Python-native path and is not claimed to numerically match every
  PhosR result.
- **The CLI is intentionally small.** It covers the core preprocessing and `predMat`-driven downstream path. The
  native kinase workflow is currently exposed through the Python API and example script.
- **R is only required for fixture regeneration.** You do not need R to install PhosPy or run the committed Python test
  suite.

## For Contributors

Most users can ignore this section.

To work from a local checkout:

```bash
pip install -e .
```

To run tests:

```bash
pip install -e ".[test]"
pytest -m "not parity"
pytest -m parity
```

To run the usual contributor checks:

```bash
pip install -e ".[dev]"
pre-commit install
pre-commit run --all-files
```

### R Requirements for Fixture Regeneration

The committed parity fixtures are already included in the repository. You only need R if you want to regenerate or
extend them.

Current R package requirements:

- `PhosR`
- `SummarizedExperiment`
- `e1071`
- `readr`
- `dplyr`
- `tidyr`
- `tibble`
- `janitor`

A practical greenfield setup is:

```r
install.packages(c("BiocManager", "devtools", "e1071", "readr", "dplyr", "tidyr", "tibble", "janitor"))
BiocManager::install("SummarizedExperiment")
devtools::install_github("PYangLab/PhosR")
```

To regenerate the committed R reference fixtures:

```bash
Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R
```

## Attribution

All scientific credit for the original methods, package design, and biological workflow belongs to the PhosR authors
and maintainers.

Please cite and acknowledge the original PhosR work when using this repository:

- Kim, H. J., Kim, T., Hoffman, N. J., Xiao, D., James, D. E., Humphrey, S. J., & Yang, P. (2021). *PhosR enables
  processing and functional analysis of phosphoproteomic data*. Cell Reports, 34(8), 108771.
- Kim, H., Kim, T., Xiao, D., & Yang, P. (2021). *Protocol for the processing and downstream analysis of
  phosphoproteomic data with PhosR*. STAR Protocols, 2(2), 100585.
- Original R package: `PYangLab/PhosR`

PhosPy should be described as an unofficial implementation unless and until the original PhosR authors choose to
endorse or participate in it.
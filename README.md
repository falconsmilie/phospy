# PhosPy

`PhosPy` v1 is an unofficial Python implementation of selected PhosR-style workflows for phosphoproteomics.

The package is deliberately narrow. v1 supports:

- **core preprocessing** from total and phospho input tables to corrected phosphosite matrices
- **downstream kinase analysis from `predMat`**, including weighted activity, KSEA-style summaries, and target counts
- **a native kinase workflow** for profile construction, motif scoring, score combination, candidate selection, and adaptive SVM prediction

PhosPy is useful as a Python-native, test-backed subset. It is **not** presented as a full replacement for the R/Bioconductor `PhosR` package.

## Known Limitations

Before you adopt PhosPy, the important boundaries are:

- **Selective scope only.** v1 covers the workflows documented below and nothing broader.
- **Parity is seam-level, not package-wide.** Validation claims are limited to the committed fixture-backed seams described in [`docs/validation-and-parity.md`](docs/validation-and-parity.md) and [`docs/parity.md`](docs/parity.md).
- **`KinaseWorkflow` is native first.** It includes an `svm_mode="r_parity"` option for narrower learner-seam comparisons, but the default mode is the preferred Python-native path and is not claimed to numerically match every PhosR result.
- **The CLI is intentionally small.** It covers the core preprocessing and `predMat`-driven downstream path. The native kinase workflow is currently exposed through the Python API and example script.
- **R is only required for fixture regeneration.** Running the committed Python test suite does not require R, but regenerating the R reference fixtures does.

## Install

### Python Requirements

Base install:

```bash
pip install -e .
```

For tests:

```bash
pip install -e ".[test]"
pytest -m "not parity"
pytest -m parity
```

For development checks:

```bash
pip install -e ".[dev]"
pre-commit install
pre-commit run --all-files
```

### R Requirements for Fixture Regeneration

The committed parity fixtures are already included in the repository. You only need R when you want to regenerate or extend
those fixtures.

The current scripts use these R packages:

- `PhosR`
- `SummarizedExperiment`
- `e1071`
- `readr`
- `dplyr`
- `tidyr`
- `tibble`
- `janitor`

A practical greenfield setup in R is:

```r
install.packages(c("BiocManager", "devtools", "e1071", "readr", "dplyr", "tidyr", "tibble", "janitor"))
BiocManager::install("SummarizedExperiment")
devtools::install_github("PYangLab/PhosR")
```

The bundled fixture scripts check for the packages they need and stop with a clear error if anything is missing.

To regenerate the committed R reference fixtures:

```bash
Rscript scripts/generate_r_fixtures.R
Rscript scripts/generate_r_l6_fixtures.R
```

## Supported Public API for v1

The supported root-level public API is intentionally small:

- `PhosphoDataset`
- `PhosRPipeline`
- `KinaseActivityAnalyzer`
- `KinaseWorkflow`
- result dataclasses returned by those classes:
  - `CoreProcessingResult`
  - `SiteMatrixResult`
  - `CoreOutputs`
  - `KinaseActivityResult`
  - `KinasePredictionResult`
  - `KinaseWorkflowResult`

Examples below use only those supported root imports. Lower-level submodule imports may still exist for internal use and
testing, but they are not part of the stable public API unless documented here.

## Quick Start

The quickest path is to use the bundled example data under `examples/data/`.

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

If your analysis needs pairwise comparisons, pass them explicitly:

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

This writes the core CSV outputs plus downstream kinase-analysis tables, including `kinase_target_table.csv`.

### Native End-to-End Kinase Workflow

A complete runnable native-workflow example lives at [`examples/native_workflow_demo.py`](examples/native_workflow_demo.py):

```bash
python examples/native_workflow_demo.py
```

That example uses only the supported v1 root API and prints a small prediction matrix for a synthetic two-kinase setup.

## CLI Demo

A small synthetic dataset is included.

After installation, you can run:

```bash
phospy \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --pred-mat examples/data/predMat.csv \
  --max-unmatched-fraction 0.1 \
  --outdir examples/output
```

The example output directory under `examples/output/` shows the generated CSV files.

`--max-unmatched-fraction` defaults to `0.0`, which means protein correction fails if the inner join would silently drop any phosphosite rows. Raise it only when you deliberately want to allow a bounded amount of row loss.

## Testing, Validation, and Release Gate

The v1 release gate is intentionally simple:

```bash
pre-commit run --all-files
pytest -m "not parity"
pytest -m parity
```

That gate covers:

- linting and formatting via `pre-commit`
- the regular Python test suite, including the documented example smoke test
- the parity suite against the committed R-backed fixtures

Supporting documentation:

- [`docs/validation-and-parity.md`](docs/validation-and-parity.md) explains the validation layers, release gate, and test commands
- [`docs/parity.md`](docs/parity.md) explains what parity means here, especially for the native kinase workflow
- [`docs/fixtures.md`](docs/fixtures.md) explains the fixture and trace directories, generation commands, and which outputs are committed reference data
- [`CONTRIBUTING.md`](CONTRIBUTING.md) covers local setup, linting, tests, and CI expectations
- [`CHANGELOG.md`](CHANGELOG.md) contains the draft v1 release notes

## Attribution

All scientific credit for the original methods, package design, and biological workflow belongs to the PhosR authors and
maintainers.

Please cite and acknowledge the original PhosR work when using this repository:

- Kim, H. J., Kim, T., Hoffman, N. J., Xiao, D., James, D. E., Humphrey, S. J., & Yang, P. (2021). *PhosR enables
  processing and functional analysis of phosphoproteomic data*. Cell Reports, 34(8), 108771.
- Kim, H., Kim, T., Xiao, D., & Yang, P. (2021). *Protocol for the processing and downstream analysis of
  phosphoproteomic data with PhosR*. STAR Protocols, 2(2), 100585.
- Original R package: `PYangLab/PhosR`

PhosPy should be described as an unofficial implementation unless and until the original PhosR authors choose to endorse
or participate in it.

## License

This repository is distributed under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [`LICENSE`](LICENSE).

That choice is deliberate. PhosR is distributed under GPL-3, and the GNU GPL FAQ treats translation of a program into
another programming language as a kind of modification or translation under copyright law. This project therefore uses
GPL-3.0-only as the conservative licensing position for a Python implementation of selected PhosR-style workflows.

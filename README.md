# PhosPy

`PhosPy` is an unofficial Python implementation of selected PhosR-style workflows for phosphoproteomics.

It is designed for people who want a small, Python-native way to:

- preprocess phosphoproteomics tables
- analyse kinase activity from `predMat`
- run a native kinase workflow from scoring through prediction

PhosPy is deliberately narrow. It is **not** a full replacement for the R `PhosR` package.

## Install

PhosPy supports Python 3.10 and newer.

Install the supported Python API and the `phospy` CLI:

```bash
pip install phospy
```

The file-path examples below use `examples/data/...`, so they assume you are working from a
repository checkout. If you installed from PyPI, use the same code with paths to your own input files.

## What You Can Do With PhosPy

### Preprocess Phosphoproteomics Data

Start from total and phospho input tables and produce corrected phosphosite matrices for downstream use.

### Analyse Kinase Activity From `predMat`

Generate weighted activity scores, KSEA-style summaries, and target counts from predicted kinase–substrate
relationships.

### Run a Native Kinase Workflow

Construct substrate profiles, score motifs, combine evidence, select candidates, and perform adaptive SVM-based kinase
prediction.

## Supported Public API

The stable API is intentionally small:

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

For a compact guide to the supported classes, methods, and result objects, see [`docs/api.md`](docs/api.md).

## Input Tables at a Glance

PhosPy expects a small, fixed set of input shapes.

### Total-proteome table

Required columns:

- `genes`
- `group1` to `group6`

### Phosphoproteome table

Required columns:

- `uid`
- `gene_names`
- `gene_p_site`
- `localization_prob`
- `centralized_sequence`
- `p_group1` to `p_group6`

`gene_p_site` must look like `GENE_SITE`, for example `PRKACA_S339`.

### `predMat`

`predMat` must be a numeric matrix with:

- phosphosite IDs as the index, for example `BTK;Y551;`
- kinase names as columns
- scores in the range `[0, 1]`

When you load tables from files, PhosPy normalises input headers to lowercase snake case before validation. For example,
`Gene Names` and `gene-names` both become `gene_names`. That makes file input a little more forgiving, but it also
means loading fails if two raw headers collapse to the same cleaned name.

If you build `PhosphoDataset` from in-memory pandas data frames instead, those column names are validated as provided.

## Quick Start

The quickest way to get started from a source checkout is to use the bundled example data in `examples/data/`.

### Core Preprocessing

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

site_matrix = core.site_matrix.matrix
corrected = core.phospho_corrected
```

For the bundled example data, `site_matrix.index.tolist()` returns `['BTK;Y551;']`.

`dataset.preprocessing` is the bound preprocessing facade for the dataset. `process_core(...)` remains available, but `dataset.preprocessing.run(...)` is the preferred public entrypoint.

`dataset.preprocessing.run()` returns a `CoreProcessingResult` with:

- `total_unique`
- `total_filtered`
- `phospho_filtered`
- `phospho_corrected`
- `site_matrix`

If your analysis needs explicit pairwise comparisons, pass them when you build the dataset:

```python
from phospy import PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
    comparisons=[("group1", "group4"), ("group2", "group5")],
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)
```

If you do not pass `comparisons`, preprocessing still runs normally and no extra pairwise columns are added.

If you only want the phosphosite localisation filter as a standalone preprocessing step, use the public helper in
`phospy.preprocessing`:

```python
from phospy.preprocessing import filter_localized_sites

filtered = filter_localized_sites(phospho_df, threshold=0.75)
summary_result = filter_localized_sites(
    phospho_df,
    threshold=0.75,
    return_summary=True,
)
```

`summary_result.filtered` contains the retained rows and `summary_result.summary` reports how many rows were kept or
removed.

If you want to filter by observed data coverage before the broader workflow, use the standalone coverage helper:

```python
from phospy.preprocessing import filter_sites_by_coverage

coverage_result = filter_sites_by_coverage(
    phospho_df,
    columns=["p_group1", "p_group2", "p_group3", "p_group4", "p_group5", "p_group6"],
    min_coverage=0.5,
    return_summary=True,
)
```

`filter_localized_sites(...)` removes sites with weak localisation evidence, while
`filter_sites_by_coverage(...)` removes sites with too many missing sample values. Coverage filtering currently
operates across the sample columns you provide rather than from a separate group metadata model.

### Downstream Kinase Analysis From `predMat`

`KinaseActivityAnalyzer` is the public orchestration layer for standalone downstream kinase analysis. Use it when you
already have a phosphosite matrix and a `predMat` and want the downstream kinase summary tables without going through
`PhosRPipeline`.

```python
from phospy import KinaseActivityAnalyzer, PhosphoDataset

dataset = PhosphoDataset.from_files(
    "examples/data/total.tsv",
    "examples/data/phospho.tsv",
    phospho_encoding="utf-16le",
)
core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

analyzer = KinaseActivityAnalyzer()
kinase = analyzer.load_and_analyze(
    pred_mat_path="examples/data/predMat.csv",
    phospho_matrix=core.site_matrix.matrix,
    threshold=0.6,
    min_substrates=1,
    top_n_substrates=1,
)
analyzer.write_outputs(kinase, outdir="examples/output")

target_counts = kinase.target_counts
ksea_scores = kinase.ksea_scores
```

The bundled example uses `min_substrates=1` and `top_n_substrates=1` because the example matrix is intentionally tiny.
For larger real datasets, the defaults (`min_substrates=3`, `top_n_substrates=20`) are usually the better starting
point.

For the bundled example data, `target_counts.to_dict()` is `{'PRKACA': 3, 'BTK': 2}`.

`KinaseActivityAnalyzer.load_and_analyze(...)` returns a `KinaseActivityResult` with:

- `weighted_activity`
- `ksea_scores`
- `ksea_counts`
- `target_counts`
- `target_table`

### End-to-End Pipeline

```python
from phospy import PhosRPipeline

pipeline = PhosRPipeline.from_files(
    total_path="examples/data/total.tsv",
    phospho_path="examples/data/phospho.tsv",
    pred_mat_path="examples/data/predMat.csv",
    phospho_encoding="utf-16le",
    max_unmatched_fraction=0.1,
)
outputs = pipeline.run(outdir="examples/output")
```

`outputs` is a `CoreOutputs` object with:

- `outputs.core`
- `outputs.kinase_activity`

This writes the core CSV outputs together with downstream kinase-analysis tables, including:

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

If you omit `pred_mat_path`, the pipeline still runs the core preprocessing path and simply skips the downstream
kinase-analysis outputs.

### Native End-to-End Kinase Workflow

A complete runnable native-workflow example is included at
[`examples/native_workflow_demo.py`](examples/native_workflow_demo.py).

If PhosPy is installed in the environment, for example with `pip install phospy` or `pip install -e .` from a local
checkout, you can run it directly:

```bash
python examples/native_workflow_demo.py
```

From a local checkout, there is also a Make target that runs the example with the repository `src/` path configured for
that shell session:

```bash
make native-workflow-demo
```

That example uses only the supported API and prints a small prediction matrix for a synthetic two-kinase setup.

The native workflow expects:

- a phosphosite matrix
- a `substrate_map`
- `site_sequences` keyed by phosphosite ID when motif scoring is used
- `motif_sequences` for end-to-end motif-aware prediction

`site_sequences` can be passed as either a mapping keyed by phosphosite ID or a pandas Series with a phosphosite index.
If you want profile-only prediction, pass `allow_profile_only_fallback=True` and omit `motif_sequences`.

## Command-Line Demo

After installation, you can run the CLI on your own files. The example below uses the bundled tables from a source
checkout:

```bash
phospy \
  --total examples/data/total.tsv \
  --phospho examples/data/phospho.tsv \
  --pred-mat examples/data/predMat.csv \
  --phospho-encoding utf-16le \
  --max-unmatched-fraction 0.1 \
  --outdir examples/output
```

The example output directory in `examples/output/` shows the generated CSV files.

The CLI currently supports these options:

- `--total` and `--phospho` are required input files
- `--phospho-encoding` optionally overrides the default `utf-8` reader encoding
- `--outdir` is the required output directory
- `--pred-mat` is optional
- `--localization-threshold` defaults to `0.75`
- `--min-observed` defaults to `4`
- `--total-sentinel` defaults to `10.0`
- `--phospho-sentinel` defaults to `12.0`
- `--max-unmatched-fraction` defaults to `0.0`

`--max-unmatched-fraction=0.0` means protein correction fails if the inner join would silently drop any phosphosite
rows. Raise it only when you want to allow a small, bounded amount of row loss.

The CLI is intentionally small. It does not currently expose pairwise comparison generation or the native
`KinaseWorkflow` path.

## Validation Rules Worth Knowing

A few checks are especially useful to know up front:

- `localization_prob` must stay within `[0, 1]`.
- `predMat` values must stay within `[0, 1]`.
- file-loaded total and phospho headers are cleaned to lowercase snake case before validation, so duplicate cleaned
  names are rejected.
- `predMat` and the phosphosite matrix must overlap by at least one phosphosite row, and that overlap must cover at
  least 10% of the phosphosite matrix.
- Protein correction normalises gene identifiers before matching and, by default, refuses to drop unmatched phosphosite
  rows.
- Site-matrix construction drops rows with missing sequences or incomplete corrected values, then deduplicates repeated
  phosphosites by keeping the row with the highest mean corrected signal.
- In the native workflow, `motif_sequences` require matching `site_sequences`. If you omit motif data entirely, set
  `allow_profile_only_fallback=True`.

## Where to Go Next

If you want more detail, these are the most useful follow-on docs:

- [`docs/api.md`](docs/api.md) maps the supported public API
- [`docs/validation-and-parity.md`](docs/validation-and-parity.md) explains how validation is approached in PhosPy
- [`docs/parity.md`](docs/parity.md) explains what parity means here, especially for the native kinase workflow
- [`docs/fixtures.md`](docs/fixtures.md) maps the committed fixture and trace directories
- [`docs/roadmap.md`](docs/roadmap.md) outlines the most likely next steps
- [`CHANGELOG.md`](CHANGELOG.md) contains the release notes

If you want to contribute or work from a local checkout, see [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Known Limitations

A few boundaries are worth knowing up front:

- **Selective scope only.** PhosPy covers the workflows documented above and nothing broader.
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

If you want the parity suite to print its optional comparison metrics while you debug a seam, these environment
variables are available:

- `PHOSPY_SHOW_PARITY`: master switch for parity metrics output
- `PHOSPY_SHOW_PROFILE_CONSTRUCTION`: also print the optional profile-construction metrics
- `PHOSPY_SHOW_PREDICTION_MODE_COMPARISON`: also print default-versus-`r_parity` prediction comparison metrics
- `PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON`: also print replayed prediction comparison metrics

The three more specific flags only do anything when `PHOSPY_SHOW_PARITY` is enabled first. Truthy values are
case-insensitive and include `1`, `true`, `yes`, and `on`.

To see the printed summaries in the terminal, run pytest with `-s` (or `--capture=no`). If you enable all four flags
and run the full parity suite, PhosPy prints every available metrics block reached by those tests.

Linux or macOS quick example:

```bash
PHOSPY_SHOW_PARITY=1 PHOSPY_SHOW_PROFILE_CONSTRUCTION=1 PHOSPY_SHOW_PREDICTION_MODE_COMPARISON=1 PHOSPY_SHOW_REPLAYED_PREDICTION_MODE_COMPARISON=1 pytest -m parity -s
```

For Linux/macOS, PowerShell, and Command Prompt examples together with a sample of the bundled parity output, see
[`docs/parity.md`](docs/parity.md).

To run the usual contributor checks:

```bash
pip install -e ".[dev]"
pre-commit install
pre-commit run --all-files
```

### R Requirements for Fixture Regeneration

The committed parity fixtures are already included in the repository. You only need R if you want to regenerate or
extend them.

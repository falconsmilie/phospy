# phosrpy

`phosrpy` is an **unofficial Python port** of selected PhosR workflow components for phosphoproteomics.

The original **PhosR** project is an R/Bioconductor package from the PhosR team / Yang Lab for phosphoproteomic data processing and downstream kinase and signalling analysis. This repository aims to translate parts of that workflow into Python while keeping attribution clear and keeping scope claims honest.

## Attribution

All scientific credit for the original methods, package design, and biological workflow belongs to the **PhosR authors and maintainers**.

Please cite and acknowledge the original PhosR work when using this repository:

- Kim, H. J., Kim, T., Hoffman, N. J., Xiao, D., James, D. E., Humphrey, S. J., & Yang, P. (2021). *PhosR enables processing and functional analysis of phosphoproteomic data*. Cell Reports, 34(8), 108771.
- Kim, H., Kim, T., Xiao, D., & Yang, P. (2021). *Protocol for the processing and downstream analysis of phosphoproteomic data with PhosR*. STAR Protocols, 2(2), 100585.
- Original R package: `PYangLab/PhosR`

This Python repository should be described as an **unofficial port** unless and until the original PhosR authors choose to endorse or participate in it.

## License

This repository is distributed under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [`LICENSE`](LICENSE).

That choice is deliberate. PhosR is distributed under GPL-3, and the GNU GPL FAQ treats translation of a program into another programming language as a kind of modification or translation under copyright law. This project therefore uses GPL-3.0-only as the conservative licensing position for a Python port.

## Current scope

Implemented now:

- collapse duplicate total-protein genes by highest mean signal
- replace sentinel missing values
- minimum-observation filtering
- phosphosite correction against total proteome
- pairwise comparison columns from corrected phosphosites
- site-matrix construction with duplicate site collapse
- weighted kinase activity from `predMat`
- KSEA-style scores and kinase target counts from `predMat`
- a minimal CLI
- class-based public API
- parity-test harness for comparing Python outputs against R-generated fixtures

Not implemented yet:

- native Python replacement for `kinaseSubstrateScore()`
- native Python replacement for `kinaseSubstratePred()`
- native Python replacement for `Signalomes()`
- full numerical parity claims against PhosR for the unported methods

## Design notes

The package core keeps only structural defaults such as `group1` to `group6` and generic comparison pairs. Study-specific labels or experimental naming should live in caller code, example configuration, or fixture-generation scripts rather than inside `src/phosrpy`.

## Install

```bash
pip install -e .
```

For tests:

```bash
pip install -e ".[test]"
pytest
```

## Public API

```python
from phosrpy import PhosphoDataset, KinaseActivityAnalyzer, PhosRPipeline

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.process_core()

analyzer = KinaseActivityAnalyzer.from_csv("predMat.csv")
kinase = analyzer.analyze(core.site_matrix.matrix)

pipeline = PhosRPipeline.from_files(
    total_path="total.tsv",
    phospho_path="phospho.tsv",
    pred_mat_path="predMat.csv",
)
outputs = pipeline.run(outdir="output")
```

You can also override the default structural comparison pairs explicitly:

```python
dataset = PhosphoDataset(
    total_df=total_df,
    phospho_df=phospho_df,
    comparisons=[("group1", "group2"), ("group1", "group3")],
)
```

Legacy four-field comparison tuples are still accepted for compatibility, but only the first two group names are used by the Python core.

## Minimal demo

A small synthetic dataset is included:

```bash
PYTHONPATH=src python -m phosrpy.cli \
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

## R reference fixtures and parity tests

This repository has two fixture paths.

### 1. Small synthetic fixture path

Use this for deterministic preprocessing and core matrix-building parity:

```bash
Rscript scripts/generate_r_fixtures.R
```

This writes CSV fixtures into `tests/fixtures/r_reference/` for:

- corrected phosphosite values
- PhosR input rows and site matrix
- `predMat`
- weighted kinase activity
- KSEA scores
- substrate counts
- `sessionInfo()` for provenance

This path is useful for logic-level parity, but it is not strong evidence for downstream scoring equivalence beyond the implemented wrapper flow.

### 2. Richer bundled PhosR L6 fixture path

Use this for a more realistic downstream kinase-analysis parity path based on PhosR’s bundled rat L6 myotube example dataset, which is used throughout the original package examples and vignette.

Generate those fixtures with:

```bash
Rscript scripts/generate_r_l6_fixtures.R
```

This writes CSV fixtures into `tests/fixtures/r_reference_l6/` for:

- the filtered standardised L6 phosphosite matrix used for kinase analysis
- `predMat`
- weighted kinase activity
- KSEA scores
- kinase target counts
- `sessionInfo()` for provenance

This path is the better current evidence for parity of the implemented downstream kinase-analysis methods.

If the fixtures are present, parity tests can be run with:

```bash
pytest -m parity
```

## Honest project status

Today this repository is best described as a structured Python package for PhosR-style preprocessing and downstream kinase-analysis summaries, with a live R-backed parity harness and a clear path toward porting the core PhosR scoring methods.

It is not yet a full Python replacement for PhosR.

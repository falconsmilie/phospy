# phosrpy

A small Python package that reproduces a practical subset of the outputs used in a PhosR-shaped phosphoproteomics workflow.

The package now exposes a **class-based API** around the core domain concepts instead of centring everything on one pipeline file:

- `PhosphoDataset` owns input tables and preprocessing
- `KinaseActivityAnalyzer` owns downstream summaries derived from `predMat`
- `PhosRPipeline` orchestrates both when you want the full run

This is **not** yet a full Python replacement for PhosR. It is the start of a package that can be open-sourced and grown in a disciplined way.

## Current scope

Implemented now:

- collapse duplicate total-protein genes by highest mean signal
- replace sentinel missing values
- minimum-observation filtering
- phosphosite correction against total proteome
- pairwise comparison columns from corrected phosphosites
- site-matrix construction with duplicate site collapse
- weighted kinase activity from `predMat`
- KSEA scores and kinase target counts from `predMat`
- a minimal CLI
- class-based public API

Not implemented yet:

- native Python replacement for `kinaseSubstrateScore()`
- native Python replacement for `kinaseSubstratePred()`
- native Python replacement for `Signalomes()`

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

## Minimal demo

A small synthetic dataset is included:

```bash
PYTHONPATH=src python -m phosrpy.cli   --total examples/data/total.tsv   --phospho examples/data/phospho.tsv   --pred-mat examples/data/predMat.csv   --outdir examples/output
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

## Generating R reference fixtures

The next milestone is parity against fixed R outputs. A small R harness is included to generate those fixtures from the synthetic example data using real PhosR.

From the repository root:

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

Those outputs are meant to be checked into version control and used by later parity tests.

## Notes

The code is intentionally split into domain classes with small helper modules underneath so parity against R can be tested method by method.

The next sensible milestone is to add **golden parity tests** against fixed R outputs for:

- corrected phosphosite values
- site matrix rows / IDs
- kinase activity scores
- KSEA scores

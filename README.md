# phosrpy

`phosrpy` is an **unofficial Python port** of selected PhosR workflow components for phosphoproteomics.

The original **PhosR** project is an R/Bioconductor package from the PhosR team / Yang Lab for phosphoproteomic data processing and downstream kinase and signalling analysis. This repository is intended to translate parts of that workflow into Python while preserving clear attribution to the original project and staying honest about which pieces are already ported and which are not.

## Attribution

All scientific credit for the original methods, package design, and biological workflow belongs with the **PhosR authors and maintainers**.

Please cite and acknowledge the original PhosR work when using this repository:

- Kim, H. J., Kim, T., Hoffman, N. J., Xiao, D., James, D. E., Humphrey, S. J., & Yang, P. (2021). *PhosR enables processing and functional analysis of phosphoproteomic data*. Cell Reports, 34(8), 108771.
- Kim, H., Kim, T., Xiao, D., & Yang, P. (2021). *Protocol for the processing and downstream analysis of phosphoproteomic data with PhosR*. STAR Protocols, 2(2), 100585.
- Original R package: `PYangLab/PhosR`

This Python repository should be described as an **unofficial port** unless and until the original PhosR authors choose to endorse or participate in it.

## License

This repository is distributed under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [`LICENSE`](LICENSE).

That choice is deliberate. PhosR is distributed under GPL-3, and the GNU GPL FAQ states that translating a program into another programming language is treated as a kind of modification/translation under copyright law. For that reason, this project uses GPL-3.0-only as the conservative licensing position for a Python port.

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

## R reference fixtures and parity tests

A small R harness is included to generate reference fixtures from the synthetic example data using PhosR:

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

If the fixtures are present, parity tests can be run with:

```bash
pytest -m parity
```

Until the core scoring and prediction methods are ported, this repository should not claim to be a full PhosR replacement.

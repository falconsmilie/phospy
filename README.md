# PhosPy

`PhosPy` is an **unofficial Python port** of selected PhosR workflow components for phosphoproteomics.

The original **PhosR** project is an R/Bioconductor package from the PhosR team / Yang Lab for phosphoproteomic data 
processing and downstream kinase and signalling analysis. This repository translates parts of that workflow into Python, 
including a native end-to-end kinase scoring and prediction path, while keeping attribution clear and parity claims 
limited to fixture-backed seams.

## Attribution

All scientific credit for the original methods, package design, and biological workflow belongs to the **PhosR authors 
and maintainers**.

Please cite and acknowledge the original PhosR work when using this repository:

- Kim, H. J., Kim, T., Hoffman, N. J., Xiao, D., James, D. E., Humphrey, S. J., & Yang, P. (2021). *PhosR enables processing and functional analysis of phosphoproteomic data*. Cell Reports, 34(8), 108771.
- Kim, H., Kim, T., Xiao, D., & Yang, P. (2021). *Protocol for the processing and downstream analysis of phosphoproteomic data with PhosR*. STAR Protocols, 2(2), 100585.
- Original R package: `PYangLab/PhosR`

This Python repository should be described as an **unofficial port** unless and until the original PhosR authors choose 
to endorse or participate in it.

## License

This repository is distributed under the **GNU General Public License v3.0 only (GPL-3.0-only)**. See [`LICENSE`](LICENSE).

That choice is deliberate. PhosR is distributed under GPL-3, and the GNU GPL FAQ treats translation of a program into 
another programming language as a kind of modification or translation under copyright law. This project therefore uses 
GPL-3.0-only as the conservative licensing position for a Python port.

## Current scope

Implemented now:

- collapse duplicate total-protein genes by highest mean signal after gene normalisation
- replace sentinel missing values
- minimum-observation filtering
- phosphosite correction against total proteome
- optional pairwise comparison columns from corrected phosphosites
- site-matrix construction with duplicate site collapse
- weighted kinase activity from `predMat`
- KSEA-style scores and kinase target counts from `predMat`
- a minimal CLI
- class-based public API
- parity-test harness for comparing Python outputs against R-generated fixtures
- native kinase substrate-profile construction via `KinaseProfileBuilder`
- native profile-based kinase scoring via `KinaseScorer.score_phosphosite_profiles()`
- native motif-frequency scoring via `KinaseMotifScorer` and `score_phosphosite_motifs()`
- native rank-weighted motif/profile score combination via `combine_profile_and_motif_scores()`
- native candidate-substrate selection and adaptive SVM prediction via `KinasePredictor`
- a dedicated end-to-end native orchestration layer via `KinaseWorkflow` / `run_kinase_workflow()`

Not implemented yet:

- full native Python replacement for every PhosR method
- native Python replacement for `Signalomes()`
- broad numerical parity claims for the newer native kinase workflow beyond the fixture-backed seams documented in `docs/parity.md`
- full feature, behaviour, and output parity with the original R package

## Design notes

The package core keeps only structural defaults such as `group1` to `group6`. Comparison choices are analysis intent, 
so they should be supplied by caller code or live in example configuration and fixture-generation scripts rather than 
inside `src/phospy`.

The native kinase scoring path is intentionally modular:

- substrate-profile construction
- motif scoring
- profile scoring
- weighted score combination
- adaptive SVM prediction

Those pieces can be used independently or through the higher-level workflow API.

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

## Public API

### Core preprocessing and downstream summaries

```python
from phospy import (
    KinaseActivityAnalyzer,
    PhosphoDataset,
    PhosRPipeline,
)

dataset = PhosphoDataset.from_files("total.tsv", "phospho.tsv")
core = dataset.process_core()

# Add pairwise comparisons only when your analysis wants them.
dataset_with_comparisons = PhosphoDataset.from_files(
    "total.tsv",
    "phospho.tsv",
    comparisons=[("group1", "group4"), ("group2", "group5")],
)
core_with_comparisons = dataset_with_comparisons.process_core()

analyzer = KinaseActivityAnalyzer.from_csv("predMat.csv")
kinase = analyzer.analyze(core.site_matrix.matrix)

pipeline = PhosRPipeline.from_files(
    total_path="total.tsv",
    phospho_path="phospho.tsv",
    pred_mat_path="predMat.csv",
)
outputs = pipeline.run(outdir="output")
```

You can also request pairwise comparisons explicitly:

```python
dataset = PhosphoDataset(
    total_df=total_df,
    phospho_df=phospho_df,
    comparisons=[("group1", "group2"), ("group1", "group3")],
)
```

Without `comparisons=...`, the core pipeline does not add any pairwise comparison columns. Comparison definitions are 
plain two-tuples of group names.

### Native kinase scoring building blocks

```python
from phospy import (
    KinaseMotifScorer,
    KinasePredictor,
    KinaseProfileBuilder,
    KinaseScorer,
    combine_profile_and_motif_scores,
)

profile_builder = KinaseProfileBuilder()
profile_result = profile_builder.build(
    substrate_map={
        "PRKACA": ["SITE_1", "SITE_2"],
        "BTK": ["SITE_3", "SITE_4"],
    },
    phospho_matrix=phospho_matrix,
    min_substrates=1,
)

motif_scorer = KinaseMotifScorer.from_substrate_sequences(
    {
        "PRKACA": ["QQAAAAAYY"],
        "BTK": ["QQTTTTTYY"],
    },
    flank_size=2,
)
motif_result = motif_scorer.score_sequences(site_sequences)

scorer = KinaseScorer(profile_result.profile_matrix)
scoring_result = scorer.score(
    phospho_matrix=phospho_matrix,
    motif_scores=motif_result.motif_scores,
    motif_sizes=motif_result.motif_sizes,
    profile_sizes=profile_result.substrate_counts.astype(float),
)

predictor = KinasePredictor()
prediction_result = predictor.predict_from_scoring_result(
    scoring_result=scoring_result,
    ensemble_size=10,
    random_state=42,
)
```

### Native end-to-end kinase workflow

If you want one orchestration call instead of stitching the pieces together manually:

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

This workflow is implemented natively in Python and is intended to provide a coherent PhosR-style kinase scoring path. 
It should not be read as a blanket claim of numerical equivalence to the R package beyond the specific fixture-backed 
seams described in [`docs/parity.md`](docs/parity.md).

## Minimal demo

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

## R reference fixtures and parity tests

This repository has two fixture paths. The detailed parity model, limits, and maintenance rule live in [`docs/parity.md`](docs/parity.md).

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

This path is useful for logic-level parity, but it is not strong evidence for downstream scoring equivalence beyond the 
implemented wrapper flow.

### 2. Richer bundled PhosR L6 fixture path

Use this for a more realistic downstream kinase-analysis parity path based on PhosR’s bundled rat L6 myotube example 
dataset, which is used throughout the original package examples and vignette.

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

## Development checks

Code style and linting are enforced with Ruff through `pre-commit`. The repository policy is intentionally small:

- `ruff check --fix` for linting, import sorting, and safe fixes
- `ruff format` for formatting
- `pytest` for unit tests
- `pytest -m parity` for the R-backed parity layer

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for local setup and [`docs/parity.md`](docs/parity.md) for the parity contract.

## Honest project status

Today this repository is best described as a structured Python package for PhosR-style preprocessing and downstream 
kinase-analysis summaries, with a live R-backed parity harness and a growing native kinase workflow covering profile 
construction, motif scoring, profile scoring, rank-weighted score combination, and adaptive SVM prediction.

It is still not a full Python replacement for PhosR, and the newer native workflow should be described as an evolving 
port rather than a parity-complete reimplementation.

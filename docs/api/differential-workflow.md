# Differential Workflow

This page explains the public differential workflow API.

## Purpose

`DifferentialAnalysisWorkflow` is a first-class PhosPy workflow entrypoint
exposed from top-level `phospy`. It runs moderated differential analysis from an
`AnalysisReadyPhosphoDataset` plus explicit experimental design and contrast
definitions.

The supported public import route is:

```python
from phospy import DifferentialAnalysisWorkflow
```

`from phospy.api import DifferentialAnalysis` is not a supported public route.

Differential analysis is separate from kinase and signalome workflows.

```python
differential_result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
)
```

## Imports

```python
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    MultipleTestingConfig,
    SampleDesignRecord,
)
```

## Input Contract

`DifferentialAnalysisWorkflow.run(...)` accepts `DifferentialAnalysisRequest`.

Required inputs:

- `dataset`: `AnalysisReadyPhosphoDataset`
- `design`: `ExperimentalDesign`
- `contrasts`: `tuple[Contrast, ...]`

Optional inputs:

- `allow_design_subset` (`False` by default)
- `minimum_condition_replicates` (`2` by default)
- `empirical_bayes` (`EmpiricalBayesConfig()` by default)
- `multiple_testing` (`MultipleTestingConfig()` by default)

### Matrix Shape Expectations

- `dataset.phospho` must be a numeric phosphosite-by-sample matrix.
- Rows are phosphosite features (site IDs); columns are sample IDs.
- Values are expected to be finite and analysis-ready (no hidden preprocessing
  is run inside the differential workflow).

### Design Matrix Expectations

- `ExperimentalDesign.samples[].sample_id` must align to dataset sample IDs.
- `condition` labels are required and cannot be empty.
- Duplicate `sample_id` values are rejected.
- By default, all dataset samples must appear in the design.
- With `allow_design_subset=True`, design may define a strict sample subset.

### Contrast Expectations

- Contrasts are explicit condition-vs-condition definitions (`numerator` and
  `denominator`).
- Contrast conditions must exist in the design.
- Minimum replicate requirements are enforced per condition.
- Invalid contrast definitions fail before statistical execution.

### Sample Alignment Rules

- Alignment is label-based on sample IDs, not positional.
- Reordered sample columns are allowed when labels remain consistent.
- Design/sample mismatch fails with explicit validation errors.

### Intensity Scale and Missing Values

- Differential analysis assumes valid upstream preprocessing.
- It does not build datasets, perform sequence resolution, localisation
  filtering, imputation, normalisation, or batch correction.
- Missing-value policy is inherited from the analysis-ready dataset boundary.

## Design and Contrast Requirements

- Conditions are never inferred from sample names.
- `ExperimentalDesign.samples[].sample_id` must align to dataset sample IDs.
- By default, every dataset sample must appear in design.
- Duplicate design sample IDs are rejected.
- Empty condition labels are rejected.
- Contrast conditions must exist in the design.
- Each contrast must satisfy minimum replicate counts for numerator and
  denominator conditions.
- Batch/block metadata is modeled in the contract but not yet executable in the
  differential engine; such requests fail with explicit unsupported-feature
  errors.

## Empirical-Bayes Settings

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `method` | `str` | `"standard"` | `"standard"`, `"robust"` | `standard` applies limma-style moderation. `robust` applies winsorized robust hyperparameter estimation. |
| `trend` | `bool` | `False` | `True`, `False` | Enables mean-intensity trend fitting (`limma-trend` style). |
| `winsor_tail_p` | `tuple[float, float]` | `(0.05, 0.1)` | Two probabilities in `[0, 1)` summing to `< 1.0` | Tail clipping parameters used by robust moderation. |

## Public Workflow Shape

`DifferentialAnalysisWorkflow` is a thin public shell. Internal execution
follows the
same stage pattern as other PhosPy workflows:

1. `validated = validator.run(request)`
2. `interpreted = interpreter.run(validated)`
3. `result = executor.run(interpreted)`

The design contract is validated before statistical execution. The interpreter
receives resolved matrix-ready design/contrast structures.

## Output Model

`DifferentialAnalysisWorkflow.run(...)` returns `DifferentialAnalysisResult`.

Common outputs include:

- `result.table_for("B_vs_A")` with columns:
  - `logFC`
  - `t`
  - `P.Value`
  - `adj.P.Val`
- `result.prior_diagnostics`
- `result.mean_variance_trend_diagnostics` (when trend is enabled)

Each contrast result table is row-aligned to the input site IDs.

### Statistical Method and Multiple Testing

- Effect estimates are OLS contrast estimates (`logFC`).
- Test statistics are moderated t-statistics (`t`) using empirical-Bayes
  variance moderation.
- Raw p-values are reported in `P.Value`.
- Multiple-testing-adjusted p-values are reported in `adj.P.Val`.
- Current multiple-testing policy supports Benjamini-Hochberg adjustment.

### Result Structure and Metadata

- Per-contrast result tables (`result.table_for(contrast_name)`)
- Prior diagnostics (`result.prior_diagnostics`)
- Method metadata:
  - `result.empirical_bayes_method`
  - `result.empirical_bayes_robust`
  - `result.empirical_bayes_trend`
- Trend diagnostics when enabled (`result.mean_variance_trend_diagnostics`)

## Limitations and Non-goals

- Assumes valid upstream preprocessing and quantitative inputs.
- Assumes valid design matrix, contrast definitions, and replicate structure.
- Does not resolve peptide/site ambiguity.
- Does not perform localisation filtering unless explicitly implemented upstream.
- Does not perform missing-value imputation unless explicitly implemented
  upstream.
- Does not perform batch correction unless explicitly implemented upstream.
- Statistical interpretation depends on design, contrast specification, and
  replicate structure.

## Minimal Example

```python
from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord(
            sample_id="A_1",
            condition="A",
            biological_replicate_id="A_r1",
        ),
        SampleDesignRecord(
            sample_id="A_2",
            condition="A",
            biological_replicate_id="A_r2",
        ),
        SampleDesignRecord(
            sample_id="B_1",
            condition="B",
            biological_replicate_id="B_r1",
        ),
        SampleDesignRecord(
            sample_id="B_2",
            condition="B",
            biological_replicate_id="B_r2",
        ),
    )
)
contrasts = (
    Contrast(
        name="B_vs_A",
        numerator_condition="B",
        denominator_condition="A",
    ),
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
)

print(result.table_for("B_vs_A").head())
```

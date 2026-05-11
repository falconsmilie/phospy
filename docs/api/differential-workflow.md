# Differential Workflow

This page explains the differential workflow API. Differential analysis is a
first-class downstream workflow over `AnalysisReadyPhosphoDataset`.

## Purpose

`DifferentialAnalysisWorkflow` runs moderated differential analysis from an
analysis-ready phosphosite matrix plus an explicit, typed experimental design
contract and typed contrast definitions.

It is separate from kinase and signalome workflows. Differential results can be
consumed by downstream workflows where explicitly supported.

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
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    MultipleTestingConfig,
    SampleDesignRecord,
)
```

## Request Parameters

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | None | Yes | Dataset returned by `AnalysisReadyDatasetBuilder.run(...)`. Differential analysis consumes `dataset.phospho` as the quantitative matrix. |
| `design` | `ExperimentalDesign` | None | Yes | Typed sample-level design records (`sample_id`, `condition`, and optional replicate/batch/block fields). |
| `contrasts` | `tuple[Contrast, ...]` | None | Yes | Typed condition-vs-condition contrasts. |
| `allow_design_subset` | `bool` | `False` | No | If `True`, design sample IDs may be a strict subset of dataset samples. |
| `minimum_condition_replicates` | `int` | `2` | No | Minimum required replicate count per contrast condition. |
| `empirical_bayes` | `EmpiricalBayesConfig` | `EmpiricalBayesConfig()` | No | Moderation policy (`standard` or `robust`) and optional trend settings. |
| `multiple_testing` | `MultipleTestingConfig` | `MultipleTestingConfig(method="benjamini_hochberg")` | No | Multiple-testing adjustment policy. Current release supports Benjamini-Hochberg. |

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

`DifferentialAnalysisWorkflow` follows the same public stage pattern as other
workflows:

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

## Minimal Example

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
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

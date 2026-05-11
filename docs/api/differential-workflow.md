# Differential Workflow

This page explains the differential workflow API. Differential analysis is a
first-class downstream workflow over `AnalysisReadyPhosphoDataset`.

## Purpose

`DifferentialAnalysisWorkflow` runs moderated differential analysis from an
analysis-ready phosphosite matrix plus an explicit design and contrast
definition.

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
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    EmpiricalBayesConfig,
    MultipleTestingConfig,
)
```

## Request Parameters

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | None | Yes | Dataset returned by `AnalysisReadyDatasetBuilder.run(...)`. Differential analysis consumes `dataset.phospho` as the quantitative matrix. |
| `design` | `DesignMatrix` or `pandas.DataFrame` | None | Yes | Sample-by-term design matrix. Rows are samples and must align to `dataset.phospho.columns` by sample labels. |
| `contrasts` | `ContrastMatrix` or `pandas.DataFrame` | None | Yes | Design-term-by-contrast matrix. Rows must align to `design.columns`. |
| `empirical_bayes` | `EmpiricalBayesConfig` | `EmpiricalBayesConfig()` | No | Moderation policy (`standard` or `robust`) and optional trend settings. |
| `multiple_testing` | `MultipleTestingConfig` | `MultipleTestingConfig(method="benjamini_hochberg")` | No | Multiple-testing adjustment policy. Current release supports Benjamini-Hochberg. |

## Design and Contrast Requirements

- `design.index` must represent the same sample set as `dataset.phospho.columns`.
- Sample order may differ; labels must still match as a set.
- `design` must be full column rank.
- Residual degrees of freedom must be positive.
- `contrasts.index` must match `design.columns`.
- Each contrast vector must be estimable under the resolved design.

## Empirical-Bayes Settings

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `method` | `str` | `"standard"` | `"standard"`, `"robust"` | `standard` applies limma-style moderation. `robust` applies winsorized robust hyperparameter estimation. |
| `trend` | `bool` | `False` | `True`, `False` | Enables mean-intensity trend fitting (`limma-trend` style). |
| `winsor_tail_p` | `tuple[float, float]` | `(0.05, 0.1)` | Two probabilities in `[0, 1)` summing to `< 1.0` | Tail clipping parameters used by robust moderation. |

## Multiple-Testing Settings

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `method` | `str` | `"benjamini_hochberg"` | `"benjamini_hochberg"` | Adjusts per-contrast p-values as FDR q-values. |

## Public Workflow Shape

`DifferentialAnalysisWorkflow` follows the same public stage pattern as other
workflows:

1. `validated = validator.run(request)`
2. `interpreted = interpreter.run(validated)`
3. `result = executor.run(interpreted)`

The interpreter stage prepares the statistical execution plan (alignment,
resolved design/contrast, and executability checks) before computation.

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

## Validation Behaviour

- Invalid request types fail at validator stage.
- Unknown contrast terms fail at validator stage.
- Sample/design misalignment fails at interpreter stage.
- Non-estimable contrast definitions fail at interpreter stage.
- Execution receives interpreted inputs only.

## Scientific Assumptions and Interpretation Notes

- `logFC` values are OLS contrast estimates.
- Moderation changes variance, t-statistics, and p-values, not fold-change
  estimates.
- `adj.P.Val` reports Benjamini-Hochberg adjusted p-values per contrast.
- Differential scores are statistical associations under the supplied design and
  contrast definitions; they are not direct evidence of causal regulation.

## Minimal Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    DifferentialAnalysisWorkflow,
    Organism,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)

design = pd.DataFrame(
    {"A": [1.0, 1.0, 0.0, 0.0], "B": [0.0, 0.0, 1.0, 1.0]},
    index=["A_1", "A_2", "B_1", "B_2"],
)
contrasts = pd.DataFrame({"B_vs_A": [-1.0, 1.0]}, index=["A", "B"])

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
)

print(result.table_for("B_vs_A").head())
```

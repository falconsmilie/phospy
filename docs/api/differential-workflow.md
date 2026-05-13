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

`from phospy import DifferentialAnalysis` and
`from phospy.api import DifferentialAnalysis` are not supported public routes.

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
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    SampleDesignRecord,
    TechnicalReplicatePolicy,
)
```

## Input Contract

`DifferentialAnalysisWorkflow.run(...)` accepts `DifferentialAnalysisRequest`.

Required inputs:

- `dataset`: `AnalysisReadyPhosphoDataset`
- `design`: `ExperimentalDesign`
- `contrasts`: `tuple[Contrast, ...]`

Optional inputs:

- `config` (`DifferentialAnalysisConfig()` by default), including:
  - `technical_replicate_policy` (`TechnicalReplicatePolicy.REJECT`)
  - `allow_design_subset` (`False`)
  - `minimum_condition_replicates` (`2`)
  - `empirical_bayes` (`EmpiricalBayesConfig()`)

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
- With `config.allow_design_subset=True`, design may define a strict sample subset.

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

- Differential analysis requires an established log2 phospho intensity scale
  at the `AnalysisReadyPhosphoDataset` boundary before it can emit `logFC`.
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
- Use at least two biological replicates per condition for meaningful
  differential examples and interpretation.
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

## Worked Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    IntensityScaleKind,
    Organism,
    SampleDesignRecord,
)

phospho = pd.DataFrame(
    {
        "control_rep1": [8200.0, 9100.0],
        "control_rep2": [8000.0, 9000.0],
        "treatment_rep1": [16200.0, 9150.0],
        "treatment_rep2": [15800.0, 9050.0],
    },
    index=["MAPK14;Y182;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3B"],
        "site": ["Y182", "S9"],
        "site_sequence": [
            "MPRKSLVGTPYWMNQYAVNQKQTLRDLKQEN",
            "ATMSGRPRTTSFAESSKPVQQPSAFGQAAAL",
        ],
        "protein_id": ["MAPK14", "GSK3B"],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(
                policy="log2",
                pseudocount=1.0,
            )
        ),
    )
)

assert dataset.intensity_scale_state.kind is IntensityScaleKind.LOG2
assert dataset.intensity_scale_state.is_established

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord(
            sample_id="control_rep1",
            condition="control",
            biological_replicate_id="control_r1",
        ),
        SampleDesignRecord(
            sample_id="control_rep2",
            condition="control",
            biological_replicate_id="control_r2",
        ),
        SampleDesignRecord(
            sample_id="treatment_rep1",
            condition="treatment",
            biological_replicate_id="treatment_r1",
        ),
        SampleDesignRecord(
            sample_id="treatment_rep2",
            condition="treatment",
            biological_replicate_id="treatment_r2",
        ),
    )
)
contrasts = (
    Contrast(
        name="treatment_vs_control",
        numerator_condition="treatment",
        denominator_condition="control",
    ),
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
)

print(result.table_for("treatment_vs_control").loc[:, ["logFC", "adj.P.Val"]])
```

Interpretation notes for this tiny synthetic matrix:

- `MAPK14;Y182;` has higher treatment intensity than control, so `logFC` should
  be positive for `treatment_vs_control`.
- `GSK3B;S9;` is approximately unchanged across conditions, so its `logFC`
  should be near zero.
- The example demonstrates workflow contracts and mechanics, not biological
  discovery or study-level statistical power.

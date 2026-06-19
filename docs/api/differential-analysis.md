# Differential Analysis Workflow

`DifferentialAnalysisWorkflow` helps you test explicit condition contrasts on an
`AnalysisReadyPhosphoDataset`. It keeps the scientific boundary clear: build and
preprocess the dataset first, then pass an explicit experimental design and
contrast list to the workflow.

## When to Use This Workflow

Use this workflow when you have a strict analysis-ready phosphosite matrix and
want moderated differential phosphorylation statistics for named condition
contrasts.

Good fits:

- two-condition unpaired contrasts
- explicit fixed-effect covariates such as batch, categorical covariates, or
  continuous covariates
- complete fixed-block designs when `paired_design_policy="fixed_block"` is set
  and every block covers the requested contrast

This workflow is not a dataset builder, peptide-to-site resolver, imputation
engine, broad batch-correction method, `duplicateCorrelation` implementation, or
mixed-effects model.

## Inputs

`DifferentialAnalysisRequest.dataset` must be an
`AnalysisReadyPhosphoDataset`. The dataset must already have:

- numeric phosphosite-by-sample values in `dataset.phospho`
- rows keyed by `site_key`
- complete analysis-ready values at the dataset boundary
- an established log2 intensity scale when interpreting `logFC`
- site metadata carrying `site_key`, `display_id`, `organism`,
  `protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`

`dataset.sample_metadata` may exist, but differential design semantics come from
`ExperimentalDesign`, not passive sample metadata.

Differential results report `logFC`, so build or declare a log2 dataset before
running this workflow. For example, use
`DatasetIntensityTransformConfig(policy="log2")` during dataset building when
your input values are linear intensities.

For site-level differential analysis, configure localisation at dataset build
time so low-confidence phosphosite assignments fail fast before statistics are
reported:

```python
from phospy.api import DatasetLocalisationConfig, DatasetPreprocessingConfig

preprocessing = DatasetPreprocessingConfig(
    localisation=DatasetLocalisationConfig(
        mode="require_threshold",
        confidence_column="localisation_confidence",
        min_confidence=0.75,
    )
)
```

## Request Object

Use `DifferentialAnalysisRequest`.

Important fields:

| Field | Meaning |
| --- | --- |
| `dataset` | The `AnalysisReadyPhosphoDataset` to test. |
| `design` | An `ExperimentalDesign` with sample IDs, conditions, and optional replicate, batch, block, or covariate metadata. |
| `contrasts` | A tuple of `Contrast` objects defining condition-vs-condition tests. |
| `config` | A `DifferentialAnalysisConfig` controlling replicate, paired/block, imputation, empirical-Bayes, and multiple-testing policy. |

Constructing the request records intent only. `DifferentialAnalysisWorkflow.run`
validates dataset eligibility, sample/design alignment, contrast validity,
replicate requirements, and config coherence before execution.

## Request Configuration

Use `DifferentialAnalysisConfig`.

Important fields:

| Field | Default | Notes |
| --- | --- | --- |
| `technical_replicate_policy` | `TechnicalReplicatePolicy.REJECT` | Controls explicit technical-replicate aggregation when repeated biological replicate groups are present. |
| `paired_design_policy` | `"reject"` | Use `"fixed_block"` only for complete fixed-block designs with explicit `SampleDesignRecord.block_id`. |
| `imputed_value_policy` | `"reject"` | The default rejects upstream-imputed datasets. `"withhold_imputed_features"` requires dataset-owned imputation observation metadata. |
| `imputed_value_max_fraction` | `0.0` | Threshold used by the withhold policy. |
| `allow_design_subset` | `False` | Allows an explicit subset of dataset samples when set to `True`. |
| `minimum_condition_replicates` | `2` | Minimum replicate count per contrast side. |
| `empirical_bayes` | `EmpiricalBayesConfig()` | Controls `method`, `trend`, and robust winsor tails. |
| `multiple_testing` | `MultipleTestingConfig()` | Current public method is Benjamini-Hochberg. |

Related request classes:

- `ExperimentalDesign`
- `SampleDesignRecord`
- `Contrast`
- `BatchCovariate`
- `CategoricalCovariate`
- `ContinuousCovariate`
- `EmpiricalBayesConfig`
- `MultipleTestingConfig`

## Running the Workflow

Import the workflow from top-level `phospy` and request/config classes from
`phospy.api`.

```python
from phospy import DifferentialAnalysisWorkflow
from phospy.api import DifferentialAnalysisRequest

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=contrasts,
    )
)
```

`from phospy import DifferentialAnalysis` and
`from phospy.api import DifferentialAnalysis` are not supported public routes.

## Result Object

`DifferentialAnalysisWorkflow.run(...)` returns `DifferentialAnalysisResult`.

Important fields and helpers:

| Field or helper | Meaning |
| --- | --- |
| `table_for(contrast_name)` | Defensive snapshot of one contrast result table. |
| `contrast_tables()` | Defensive snapshots for all contrast tables. |
| `residual_variance_series()` | Feature-level residual variance snapshot. |
| `posterior_residual_variance_series()` | Moderated residual variance snapshot. |
| `prior_diagnostics` | Empirical-Bayes prior diagnostics. |
| `mean_variance_trend_diagnostics` | Trend diagnostics when trend moderation is enabled. |
| `policy_provenance` | Structured design, contrast, replicate, imputation, and testing provenance. |
| `workflow_provenance` | Workflow-level execution metadata. |
| `input_dataset_preprocessing_report` | Preprocessing report carried from the input dataset when available. |

Each contrast result table is indexed by the input `site_key`. The minimum public
identity columns are `site_key`, `display_id`, `organism`,
`protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`. Typical
statistics columns include `logFC`, `t`, `P.Value`, and `adj.P.Val`.

stat-only computation payloads are internal. They are not a public scientific
result object and are not valid `DifferentialAnalysisResult` tables.

## Interpreting the Result

`logFC` is the fitted contrast estimate on the established log2 scale. `t` is a
moderated t-statistic. `P.Value` is the raw p-value and `adj.P.Val` is the
multiple-testing adjusted value for the implemented correction policy.

Repeated `display_id` values can appear when different `site_key` rows preserve
different protein context. Interpret rows by `site_key`; `display_id` is for
readability.

## Provenance and Reproducibility

Result provenance records the resolved design, contrast vectors, fixed-effect
covariates, replicate policy, empirical-Bayes settings, multiple-testing method,
imputation policy, and unsupported-design rejection policy. Table exports return
defensive in-memory snapshots; mutating them does not mutate the result object.

## Limitations

- Conditions, batches, blocks, and covariates are not inferred from sample names.
- The workflow does not run preprocessing, imputation, localisation filtering,
  sequence resolution, or batch correction.
- Fixed-effect batch terms are model covariates, not ComBat, RUV,
  `removeBatchEffect`, `duplicateCorrelation`, or mixed-effects modelling.
- `paired_design_policy="fixed_block"` is fixed-effect block modelling only. It
  is not random subject modelling.
- Upstream-imputed datasets are rejected by default. The supported withhold
  policy does not implement observed-only fitting.

## Minimal Example

```python
from phospy import DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    SampleDesignRecord,
)

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

table = result.table_for("treatment_vs_control")
print(table.loc[:, ["site_key", "display_id", "logFC", "adj.P.Val"]])
```

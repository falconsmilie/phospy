# Differential Analysis

`DifferentialAnalysisWorkflow` tests named condition contrasts for each
phosphosite in an `AnalysisReadyPhosphoDataset`.

Use this workflow when you want effect estimates, moderated statistics, and
multiple-testing-adjusted *p* values for explicit comparisons such as treated
versus control.

!!! info "At a Glance"
    **Input:** An analysis-ready log2 dataset, an `ExperimentalDesign`, and one
    or more `Contrast` objects  
    **Run:** `DifferentialAnalysisWorkflow().run(request)`  
    **Returns:** A `DifferentialAnalysisResult` with one table per contrast,
    diagnostics, provenance, and caveats

The current implementation is limited to tested design and contrast envelopes;
it is not full limma or PhosR parity. Familiar column names such as `logFC`,
`P.Value`, and `adj.P.Val` describe the output, not broad limma compatibility.

## Before You Begin

Build the dataset with
[`AnalysisReadyDatasetBuilder`](dataset-build-workflow.md). Differential
analysis requires:

- unique `site_key` rows and complete analysis-ready values;
- required site metadata, including `site_sequence`;
- sample columns that match the explicit design;
- an established log2 intensity scale;
- at least two biological replicates per contrasted condition in production
  mode; and
- localization evidence that meets the dataset policy.

Use `DatasetIntensityTransformConfig(policy="log2")` when the source matrix is
linear. A declared log2 scale with suspicious diagnostics fails by default
unless you deliberately enable the recorded override.

Configure localization while building the dataset. A low-confidence
phosphosite fails before differential fitting:

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

PhosPy does not infer conditions, replicates, batches, blocks, or covariates from
sample names or passive `dataset.sample_metadata`. Put these values in the
`ExperimentalDesign`.

The withdrawn post-hoc peptide-to-site estimate-combination route fails closed.
Future public support requires executable peptide-to-site mapping semantics, a
coherent combined effect and inferential estimand, same-experiment dependence
handling, multiple-testing semantics, and provenance. Resolve peptide evidence
to site-level sample intensities during dataset preparation instead.

## Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

phospho = pd.DataFrame(
    {
        "control_rep1": [1000.0, 900.0, 800.0],
        "control_rep2": [1050.0, 880.0, 820.0],
        "treatment_rep1": [1800.0, 930.0, 760.0],
        "treatment_rep2": [1750.0, 920.0, 740.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A", "TSC2"],
        "site": ["Y182", "S21", "S939"],
        "site_sequence": [
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
        ],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
    },
    index=phospho.index,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            localisation=DatasetLocalisationConfig(
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            ),
        ),
    )
)

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord("control_rep1", "control", "control_r1"),
        SampleDesignRecord("control_rep2", "control", "control_r2"),
        SampleDesignRecord("treatment_rep1", "treatment", "treatment_r1"),
        SampleDesignRecord("treatment_rep2", "treatment", "treatment_r2"),
    )
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=(
            Contrast(
                name="treatment_vs_control",
                numerator_condition="treatment",
                denominator_condition="control",
            ),
        ),
    )
)

print(
    result.table_for("treatment_vs_control").loc[
        :, ["display_id", "logFC", "P.Value", "adj.P.Val"]
    ]
)
```

## Request

Create a `DifferentialAnalysisRequest`.

| Parameter | Type | Required or Default | Description | Main Constraint |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Required | Dataset to test. | Must be complete, `site_key` indexed, and established as log2. |
| `design` | `ExperimentalDesign` | Required | Explicit sample and model design. | Sample IDs must match dataset columns unless subsetting is enabled. |
| `contrasts` | `tuple[Contrast, ...]` | Required | Comparisons to estimate. | Names must be unique; numerator and denominator conditions must exist and differ. |
| `config` | `DifferentialAnalysisConfig` | `DifferentialAnalysisConfig()` | Reliability, replicate, imputation, moderation, and multiple-testing policy. | Unsupported combinations fail before fitting. |

<details markdown="1">
<summary><strong>Experimental Design Parameters</strong></summary>

### `ExperimentalDesign`

| Parameter | Type | Required or Default | Description |
| --- | --- | --- | --- |
| `samples` | `tuple[SampleDesignRecord, ...]` | Required | One record per analyzed sample. Sample IDs must be unique. |
| `fixed_effects` | `tuple[FixedEffectCovariate, ...]` | `()` | Optional categorical, continuous, or batch fixed effects. |

### `SampleDesignRecord`

| Parameter | Type | Required or Default | Description |
| --- | --- | --- | --- |
| `sample_id` | `str` | Required | Dataset sample-column name. |
| `condition` | `str` | Required | Condition label used by contrasts. |
| `biological_replicate_id` | `str` or `None` | `None` | Biological replicate identity. Technical replicates never replace this requirement. |
| `technical_replicate_id` | `str` or `None` | `None` | Technical replicate identity. Requires biological replicate IDs when used. |
| `batch` | `str` or `None` | `None` | Batch value available to `BatchCovariate`. |
| `block_id` | `str` or `None` | `None` | Block identity for explicit paired-design policies: fixed block effects or duplicate correlation. |
| `covariates` | `Mapping[str, str or int or float]` | `{}` | Values for named fixed-effect covariates. |

### `Contrast`

| Parameter | Type | Required or Default | Description |
| --- | --- | --- | --- |
| `name` | `str` | Required | Stable result-table name. |
| `numerator_condition` | `str` | Required | Positive side of the contrast. |
| `denominator_condition` | `str` | Required | Negative side of the contrast. |

Contrast direction is `numerator_condition - denominator_condition`.

### Fixed-Effect Covariates

`FixedEffectCovariate` accepts `name`, `kind`, `required=True`, and
`include_in_model=True`. Prefer `CategoricalCovariate`, `ContinuousCovariate`,
or `BatchCovariate` for common cases. Batch covariates are model terms, not
batch correction.

PhosPy supports two explicit paired-design policies. `fixed_block` represents
block identity with ordinary fixed nuisance coefficients and requires every
block to have complete within-block contrast coverage. Incomplete or partially
covered blocks are rejected; PhosPy does not silently drop those blocks or
samples.

`duplicate_correlation` estimates one consensus within-block correlation by
feature-wise REML, then fits a compound-symmetry GLS model with condition terms
and supported fixed covariates only. Block IDs are retained as covariance-group
metadata and are not added as fixed block coefficients. This is not a general
mixed-effects framework, feature-specific random-effects fitting, random
slopes, or automatic policy selection. The non-block fixed-effects design must
leave more than two residual degrees of freedom for REML correlation estimation.
Simple unpaired workflows remain the default.

</details>

<details markdown="1">
<summary><strong>Differential Configuration Parameters</strong></summary>

### `DifferentialAnalysisConfig`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `reliability_profile` | `"production"` or `"exploratory_single_replicate"` | `"production"` | Selects the production lane or the explicit exploratory single-replicate lane. |
| `technical_replicate_policy` | `"reject"`, `"mean"`, or `"median"` | `"reject"` | Rejects or combines technical replicates. Combining them does not create biological replication. |
| `paired_design_policy` | `"reject"`, `"fixed_block"`, or `"duplicate_correlation"` | `"reject"` | Selects explicit paired-design handling. `fixed_block` adds block nuisance coefficients; `duplicate_correlation` estimates one consensus compound-symmetry correlation and fits GLS. |
| `imputed_value_policy` | `"reject"` or `"withhold_imputed_features"` | `"reject"` | Rejects upstream imputation or withholds affected features using dataset-owned observation metadata. |
| `imputed_value_max_fraction` | `float` | `0.0` | Maximum imputed-cell fraction for a tested feature under the withhold policy. |
| `allow_design_subset` | `bool` | `False` | Allows the design to use an intentional subset of dataset samples. |
| `allow_suspicious_declared_input_scale` | `bool` | `False` | Allows a suspicious declared log2 scale and records the override. |
| `minimum_condition_replicates` | `int` | `2` | Minimum biological replicates per contrasted condition after policy resolution. |
| `empirical_bayes` | `EmpiricalBayesConfig` | `EmpiricalBayesConfig()` | Moderation settings. |
| `multiple_testing` | `MultipleTestingConfig` | `MultipleTestingConfig()` | Per-contrast *p*-value adjustment. |

### `EmpiricalBayesConfig`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `method` | `"standard"` or `"robust"` | `"standard"` | Moderation method. |
| `trend` | `bool` | `False` | Estimates a mean-variance trend when enabled. |
| `winsor_tail_p` | `tuple[float, float]` | `(0.05, 0.1)` | Tail proportions used by robust moderation. |

### `MultipleTestingConfig`

| Parameter | Default | Supported Values |
| --- | --- | --- |
| `method` | `"benjamini_hochberg"` | `"none"`, `"benjamini_hochberg"`, `"bonferroni"`, `"holm"`, `"benjamini_yekutieli"` |

Adjustment is performed separately for each contrast.

</details>

## Run the Workflow

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

`from phospy import DifferentialAnalysis` is not supported. `from phospy.api
import DifferentialAnalysis` is also not supported.

The workflow is deterministic for the same inputs, configuration, and package
version. It raises `WorkflowValidationError` before fitting when the dataset,
design, contrasts, replication, scale, imputation state, or covariates do not
meet the selected policy.

## Response

`DifferentialAnalysisWorkflow.run(...)` returns a
`DifferentialAnalysisResult`.

| Attribute or Helper | Format | Meaning |
| --- | --- | --- |
| `table_for(name)` | `pandas.DataFrame` | Independent snapshot of one contrast table. |
| `contrast_tables` | `dict[str, pandas.DataFrame]` | Independent snapshots of all contrast tables. |
| `feature_eligibility` | `pandas.DataFrame` or `None` | Row-level tested or withheld status. |
| `residual_variance`, `residual_variance_series()` | `pandas.Series` | Feature-level residual variance; use the helper for an independent snapshot. |
| `posterior_residual_variance`, `posterior_residual_variance_series()` | `pandas.Series` | Moderated residual variance; use the helper for an independent snapshot. |
| `prior_residual_variance`, `prior_residual_variance_series()` | `pandas.Series` | Prior residual variance by feature. |
| `prior_degrees_of_freedom_series_value`, `prior_degrees_of_freedom_series()` | `pandas.Series` | Prior degrees of freedom by feature. |
| `prior_variance`, `prior_degrees_of_freedom`, `residual_degrees_of_freedom` | Numeric | Model-wide moderation summaries. |
| `empirical_bayes_method`, `empirical_bayes_robust`, `empirical_bayes_trend` | String and booleans | Resolved empirical Bayes settings used by the fit. |
| `prior_diagnostics` | `EmpiricalBayesPriorDiagnostics` | Feature-level prior estimates and fitting diagnostics. |
| `mean_variance_trend_diagnostics` | `MeanVarianceTrendDiagnostics` or `None` | Trend diagnostics when trend moderation is enabled. |
| `diagnostics` | `DifferentialModelDiagnostics` | Design, contrast, scale, and model diagnostics. |
| `policy_provenance`, `workflow_provenance` | Typed or mapping-like provenance | Resolved scientific policy and execution metadata. |
| `caveats` | `tuple[ResultCaveat, ...]` | Structured interpretation limits. |
| `input_dataset_preprocessing_report` | Report or `None` | Dataset preprocessing report carried into the result. |
| `to_payload()` | JSON-compatible mapping | Serializable tables, diagnostics, provenance, and caveats. |
| `scientifically_equals(...)` | `bool` | Scientific-result comparison helper. |

Each contrast result table is indexed by the input `site_key`. The minimum
public identity columns are `site_key`, `display_id`, `organism`,
`protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`.

Stat-only computation payloads are internal. A stat-only computation payload is
not a public scientific result and is not a valid differential result table.

### Contrast Table Format

| Column or Index | Meaning | Always Present? |
| --- | --- | --- |
| index / `site_key` | Unique phosphosite row identity. | Yes |
| `display_id` | Readable site label; may repeat. | Yes |
| `organism`, `protein_namespace`, `protein_identifier` | Protein-scoped identity context. | Yes |
| `gene_symbol`, `site` | Readable gene and residue-position labels. | Yes |
| `logFC` | Fitted numerator-minus-denominator contrast on the log2 scale. | Yes; may be missing for withheld rows |
| `t` | Moderated *t* statistic. | Yes; may be missing for withheld rows |
| `P.Value` | Raw *p* value. | Yes; may be missing for withheld rows |
| `adj.P.Val` | Per-contrast adjusted *p* value. | Yes; may be missing for withheld rows |
| `result_status`, `result_status_reason` | Tested or withheld status and explanation. | Present when eligibility metadata is available |

<details markdown="1">
<summary><strong>Optional Imputation and Eligibility Columns</strong></summary>

When the withhold policy is active, a contrast table can also include
`imputed_cell_count`, `observed_cell_count`, `imputed_fraction`,
`imputation_policy`, `imputation_fraction_threshold`,
`contains_imputed_cells`, `observed_only_fit`,
`residual_df_adjusted_for_imputation`, and `inferential_status`.

`feature_eligibility`, when present, is indexed by `site_key` and contains
`site_key`, `result_status`, and `result_status_reason`. Status values distinguish
tested rows from rows withheld for constant values, invalid numeric values,
high imputation, insufficient observed values, or another recorded reason.

</details>

Stat-only computation payloads are internal. They are not a public scientific
result object and are not valid `DifferentialAnalysisResult` tables.

Returned DataFrames are independent snapshots. Editing one does not change the
result object.

## Interpret the Result

`logFC` is the fitted log2 difference between numerator and denominator.
Positive values indicate higher fitted phosphorylation in the numerator
condition; negative values indicate lower fitted phosphorylation.

`P.Value` is the raw *p* value for the contrast. `adj.P.Val` is adjusted within
that contrast. Smaller values indicate stronger evidence under the fitted model,
but neither value is an effect size.

Repeated `display_id` values may represent distinct protein-scoped rows. Use
`site_key` for analysis and joins.

The workflow does not perform preprocessing, localization filtering, sequence
resolution, batch correction, random-effects modelling, or post-hoc
peptide-to-site differential aggregation. A fixed-effect batch term is not
ComBat, not RUV, not limma `removeBatchEffect` parity, and not mixed-effects
modelling.

## Common Issues

| Issue | What to Check |
| --- | --- |
| A sample is missing from the dataset. | Match every `sample_id` to a dataset column, or deliberately enable `allow_design_subset`. |
| Replication fails. | Production contrasts need at least two biological replicates per condition. Technical replicates do not count. |
| The scale is unsupported. | Build or declare a log2 dataset before running the workflow. |
| Localization fails. | Add valid localization metadata and configure `DatasetLocalisationConfig` during dataset preparation. |
| A fixed-block design fails. | Confirm that every block covers both sides of every requested contrast. |
| Imputation is rejected. | Keep the default rejection policy, or use withholding only with builder-owned observation metadata. |
| Results contain withheld rows. | Review `feature_eligibility`, caveats, and the preprocessing report before interpreting absence biologically. |

## Related Guides

- [Prepare a Dataset](dataset-build-workflow.md)
- [Scientific Interpretation and Limitations](../scientific-interpretation.md)
- [Scientific Coverage](../scientific-coverage.md)
- [Kinase Analysis](kinase.md)

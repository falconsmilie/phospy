# Differential analysis workflow

## Plain-language introduction

`DifferentialAnalysisWorkflow` tests explicit condition contrasts on an
`AnalysisReadyPhosphoDataset`.

Use it when you have phosphosite intensities for named samples and want fitted
condition differences with moderated statistics. You provide the dataset, an
explicit sample design, and one or more contrasts. The workflow returns a
`DifferentialAnalysisResult` with one result table per contrast, diagnostics,
provenance, and caveats.

Current differential analysis is limited to tested design and contrast envelopes;
it is not full limma or PhosR parity. Familiar column names such as `logFC`,
`P.Value`, and `adj.P.Val` are reporting names, not a broad limma compatibility
claim.

## Input and dataset requirements

The request expects an `AnalysisReadyPhosphoDataset`. Build it with
[`AnalysisReadyDatasetBuilder`](dataset-build-workflow.md) unless you already
own trusted analysis-ready tables.

For differential analysis, the dataset must provide:

- phosphosite rows keyed by `site_key`;
- numeric phosphosite values with sample columns matching the design;
- complete analysis-ready values at the dataset boundary;
- required site metadata, including `site_key`, `display_id`, `organism`,
  `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`, and
  required `site_sequence`;
- an established log2 intensity scale for `logFC` interpretation.

Use `DatasetIntensityTransformConfig(policy="log2")` during dataset building
when your input values are linear and you want differential `logFC` output.
Declared log2 input that carries suspicious declaration diagnostics is rejected
by default unless `allow_suspicious_declared_input_scale=True` is set
deliberately.

For site-level analysis, configure localisation before the workflow. This fails
when a low-confidence phosphosite assignment would otherwise reach the
statistics.

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

Conditions, batches, blocks, and covariates are not inferred from sample names
or passive `dataset.sample_metadata`. They come from `ExperimentalDesign`.

The withdrawn post-hoc peptide-to-site differential estimate-combination route
fails closed. Future public support requires executable peptide-to-site mapping
semantics, a coherent combined effect/inference estimand, same-experiment
peptide-estimate dependence handling, multiple-testing semantics, and
provenance semantics. Until then, resolve peptide evidence into site-level
sample intensities during dataset building.

## Minimal end-to-end example

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
        "display_id": ["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;"],
        "organism": ["rat", "rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2"],
        "localisation_confidence": [0.95, 0.94, 0.96],
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
            ),
            localisation=DatasetLocalisationConfig(
                mode="require_threshold",
                confidence_column="localisation_confidence",
                min_confidence=0.75,
            ),
        ),
    )
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
print(table.loc[:, ["site_key", "display_id", "logFC", "P.Value", "adj.P.Val"]])
```

## Request model

Use `DifferentialAnalysisRequest`.

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Required | Analysis-ready phosphosite dataset to test. | Must be `site_key` indexed, complete, and established as log2 for differential output. |
| `design` | `ExperimentalDesign` | Required | Explicit sample design used to build the model. | Design sample IDs must match dataset sample columns unless `allow_design_subset=True`. |
| `contrasts` | `tuple[Contrast, ...]` | Required | Condition comparisons to estimate. | Contrast names should be unique. Each numerator and denominator condition must exist in the design and must differ. |
| `config` | `DifferentialAnalysisConfig` | Default: `DifferentialAnalysisConfig()` | Replicate, paired/block, imputation, empirical-Bayes, and multiple-testing policy. | Invalid policy combinations fail before model execution. |

`ExperimentalDesign`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `samples` | `tuple[SampleDesignRecord, ...]` | Required | One record per sample used by the workflow. | Sample IDs must be unique. |
| `fixed_effects` | `tuple[FixedEffectCovariate, ...]` | Default: `()` | Optional fixed-effect covariates included in the model. | Covariate names must resolve in the sample records. |

`SampleDesignRecord`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `sample_id` | `str` | Required | Dataset sample column name. | Must match a dataset column, unless design subsetting is explicitly allowed. |
| `condition` | `str` | Required | Experimental condition label. | Contrasts use these labels. |
| `biological_replicate_id` | `str | None` | Default: `None` | Biological replicate identity. | Production contrasts require at least two biological replicates per contrasted condition. |
| `technical_replicate_id` | `str | None` | Default: `None` | Technical replicate identity. | Technical replicates do not count as biological replicates. Declaring them requires biological replicate IDs. |
| `batch` | `str | None` | Default: `None` | Batch label available to a `BatchCovariate`. | Batch must not make the design invalid or non-estimable. |
| `block_id` | `str | None` | Default: `None` | Fixed block or paired-design label. | Used only when `paired_design_policy="fixed_block"` and blocks cover each requested contrast completely. |
| `covariates` | `Mapping[str, str | int | float]` | Default: `{}` | Extra named covariate values. | Required covariates must be present for every analysed sample. |

`Contrast`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `name` | `str` | Required | Result-table name for this contrast. | Use stable, unique names for downstream reporting. |
| `numerator_condition` | `str` | Required | Condition on the positive side of the contrast. | Must exist in `design.samples`. |
| `denominator_condition` | `str` | Required | Condition on the negative side of the contrast. | Must exist in `design.samples` and differ from the numerator. |

Contrast direction is `numerator_condition - denominator_condition`.

`FixedEffectCovariate` and helpers:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `name` | `str` | Required | Covariate name. | `BatchCovariate` uses the name `"batch"`. |
| `kind` | `"categorical" | "continuous" | "batch"` | Required | How the covariate is encoded. | Use the helper classes `CategoricalCovariate`, `ContinuousCovariate`, or `BatchCovariate` when possible. |
| `required` | `bool` | Default: `True` | Whether every analysed sample must provide the covariate. | Missing required values fail before fitting. |
| `include_in_model` | `bool` | Default: `True` | Whether to include the covariate as a model term. | Excluded covariates may still be recorded as design metadata. |

`DifferentialAnalysisConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `reliability_profile` | `"production" | "exploratory_single_replicate"` | Default: `"production"` | Selects the supported production lane or the explicit exploratory single-replicate lane. | Production requires at least two biological replicates per contrasted condition. |
| `technical_replicate_policy` | `TechnicalReplicatePolicy` or string | Default: `"reject"` | Controls technical-replicate handling. | Accepted values are `"reject"`, `"mean"`, and `"median"`. Aggregation never turns technical replicates into biological replicates. |
| `paired_design_policy` | `"reject" | "fixed_block"` | Default: `"reject"` | Enables complete fixed-block paired designs. | Incomplete or partially covered blocks are rejected. Fixed blocks are not random effects. |
| `imputed_value_policy` | `"reject" | "withhold_imputed_features"` | Default: `"reject"` | Controls upstream-imputed datasets. | Withholding requires dataset-owned imputation observation metadata. |
| `imputed_value_max_fraction` | `float` | Default: `0.0` | Maximum imputed-cell fraction allowed for a tested feature under the withhold policy. | Must be between `0.0` and `1.0`. |
| `allow_design_subset` | `bool` | Default: `False` | Allows the design to name a subset of dataset sample columns. | Keep `False` unless the analysed sample subset is intentional. |
| `allow_suspicious_declared_input_scale` | `bool` | Default: `False` | Allows declared log2 input despite suspicious scale diagnostics. | Successful overrides are recorded in policy provenance. |
| `minimum_condition_replicates` | `int` | Default: `2` | Minimum biological replicates per contrasted condition after policy resolution. | Lower values are rejected in production mode. |
| `empirical_bayes` | `EmpiricalBayesConfig` | Default: `EmpiricalBayesConfig()` | Moderation settings. | See table below. |
| `multiple_testing` | `MultipleTestingConfig` | Default: `MultipleTestingConfig()` | P-value adjustment settings. | Applied separately to each contrast table. |

`EmpiricalBayesConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `method` | `"standard" | "robust"` | Default: `"standard"` | Empirical-Bayes moderation method. | `"robust"` uses winsor-tail settings. |
| `trend` | `bool` | Default: `False` | Whether to estimate a mean-variance trend. | Trend diagnostics are reported when enabled. |
| `winsor_tail_p` | `tuple[float, float]` | Default: `(0.05, 0.1)` | Lower and upper tail proportions for robust moderation. | Each value must be in `[0, 1)` and the sum must be less than `1`. |

`MultipleTestingConfig`:

| Parameter | Type | Required or default | Description | Constraints |
| --- | --- | --- | --- | --- |
| `method` | `"none" | "benjamini_hochberg" | "bonferroni" | "holm" | "benjamini_yekutieli"` | Default: `"benjamini_hochberg"` | Converts raw `P.Value` into `adj.P.Val`. | Correction is per contrast, not pooled across contrasts. |

## Running the workflow

Call `DifferentialAnalysisWorkflow().run(request)`.

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

`from phospy import DifferentialAnalysis` is not a supported public route.
`from phospy.api import DifferentialAnalysis` is not a supported public route.

The workflow is deterministic for the same inputs, configuration, and package
version. Configuration objects may reject invalid local values at construction.
The workflow can raise `WorkflowValidationError` before fitting when the dataset,
design, contrasts, replicate policy, intensity scale, imputation state, or
covariate structure is not supported.

## Response model and output formats

`DifferentialAnalysisWorkflow.run(...)` returns `DifferentialAnalysisResult`.

| Attribute or helper | Python type | Always present? | Meaning |
| --- | --- | --- | --- |
| `table_for(contrast_name)` | `pandas.DataFrame` | Yes, for a known contrast | Defensive snapshot of one contrast result table. |
| `contrast_tables()` | `dict[str, pandas.DataFrame]` | Yes | Defensive snapshots for all contrast tables. |
| `feature_eligibility` | `pandas.DataFrame | None` | Optional | Row-level testing/withholding status when available. |
| `residual_variance` | DataFrame-like table | Yes | Feature-level residual variance payload. Prefer `residual_variance_series()` for user code. |
| `posterior_residual_variance` | DataFrame-like table | Yes | Moderated residual variance payload. Prefer `posterior_residual_variance_series()`. |
| `prior_residual_variance` | DataFrame-like table | Yes | Prior residual variance payload. Prefer `prior_residual_variance_series()`. |
| `prior_degrees_of_freedom_series_value` | DataFrame-like table | Yes | Prior degrees-of-freedom payload. Prefer `prior_degrees_of_freedom_series()`. |
| `prior_variance` | `float` or numeric value | Yes | Empirical-Bayes prior variance summary. |
| `prior_degrees_of_freedom` | `float` or numeric value | Yes | Empirical-Bayes prior degrees of freedom. |
| `residual_degrees_of_freedom` | `int` or numeric value | Yes | Residual degrees of freedom used by the fitted model. |
| `empirical_bayes_method` | `str` | Yes | Resolved empirical-Bayes method. |
| `empirical_bayes_robust` | `bool` | Yes | Whether robust moderation was used. |
| `empirical_bayes_trend` | `bool` | Yes | Whether trend moderation was used. |
| `prior_diagnostics` | object or mapping | Yes | Empirical-Bayes prior diagnostics. |
| `mean_variance_trend_diagnostics` | object or mapping | Optional by method | Mean-variance trend diagnostics when available. |
| `diagnostics` | object | Yes | Model scope, design, contrast, scale, warning, and unsupported-assumption diagnostics. |
| `policy_provenance` | object | Yes | Resolved design, contrast, replicate, imputation, input-scale, and testing policy. |
| `workflow_provenance` | object | Yes | Workflow-level execution metadata. |
| `caveats` | `tuple` | Yes | Structured caveats, such as exploratory single-replicate execution or input-scale overrides. |
| `input_dataset_preprocessing_report` | object or `None` | Optional | Preprocessing report carried from the input dataset. |
| `to_payload()` | mapping | Yes | Serializable payload containing tables, diagnostics, provenance, and caveats. |
| `scientifically_equals(...)` | `bool` | Yes | Comparison helper for scientific/result equivalence checks. |

Each contrast result table is indexed by the input `site_key`. The minimum public
identity columns are `site_key`, `display_id`, `organism`,
`protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`.

stat-only computation payloads are internal. They are not a public scientific
result object and are not valid `DifferentialAnalysisResult` tables.

Contrast table schema:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Phosphosite row identity. | `site_key` string labels | Yes |
| `site_key` | Same row identity as the index. | string | Yes |
| `display_id` | Human-readable site label. May repeat. | string | Yes |
| `organism` | Dataset organism label. | string | Yes |
| `protein_namespace` | Protein identifier namespace. | string | Yes |
| `protein_identifier` | Protein identifier in the declared namespace. | string | Yes |
| `gene_symbol` | Display gene symbol. | string | Yes |
| `site` | Residue-position site label. | string | Yes |
| `logFC` | Fitted numerator-minus-denominator contrast on the log2 scale. | numeric; may be missing for withheld rows | Yes |
| `t` | Moderated t-statistic. | numeric; may be missing for withheld rows | Yes |
| `P.Value` | Raw p-value for the fitted contrast. | numeric in `[0, 1]`; may be missing for withheld rows | Yes |
| `adj.P.Val` | Multiple-testing adjusted p-value for this contrast. | numeric in `[0, 1]`; may be missing for withheld rows | Yes |
| `imputed_cell_count` | Number of imputed cells in the analysed sample subset. | integer-like | Only with imputation withholding metadata |
| `observed_cell_count` | Number of observed cells in the analysed sample subset. | integer-like | Only with imputation withholding metadata |
| `imputed_fraction` | Imputed-cell fraction used by the withhold policy. | numeric in `[0, 1]` | Only with imputation withholding metadata |
| `imputation_policy` | Resolved imputation policy label. | string | Only with imputation withholding metadata |
| `imputation_fraction_threshold` | Configured withhold threshold. | numeric in `[0, 1]` | Only with imputation withholding metadata |
| `contains_imputed_cells` | Whether a tested retained row contains imputed cells. | boolean-like | Only with imputation withholding metadata |
| `observed_only_fit` | Whether fitting used observed-only values. | boolean-like; current withhold policy reports `False` | Only with imputation withholding metadata |
| `residual_df_adjusted_for_imputation` | Whether residual degrees of freedom were adjusted for imputation. | boolean-like; current withhold policy reports `False` | Only with imputation withholding metadata |
| `inferential_status` | Inference status for tested/withheld rows. | string | Only with imputation withholding metadata |
| `result_status` | Row status such as `tested` or withheld status. | string | Only with feature eligibility metadata |
| `result_status_reason` | Human-readable row status reason. | string | Only with feature eligibility metadata |

`feature_eligibility`, when present:

| Column or index | Meaning | Type or format | Always present |
| --- | --- | --- | --- |
| index | Phosphosite row identity. | `site_key` string labels | Yes |
| `site_key` | Same row identity as the index. | string | Yes |
| `result_status` | Testing status. | `tested`, `withheld_all_constant`, `withheld_invalid_numeric_values`, `withheld_other`, `withheld_high_imputation`, or `withheld_insufficient_observed_values` | Yes |
| `result_status_reason` | User-facing explanation for the status. | string | Yes |

PhosPy returns defensive snapshots from table helpers. Mutating a returned
DataFrame does not mutate the result object.

## Interpreting the result

`logFC` is the fitted difference between the numerator and denominator
conditions on the established log2 scale. Positive `logFC` means the numerator
condition is higher than the denominator for that phosphosite; negative means it
is lower.

`P.Value` is the raw p-value for the tested contrast. `adj.P.Val` is adjusted
within that contrast only. Smaller adjusted values indicate stronger statistical
evidence after multiple testing, but they are not effect sizes.

Repeated `display_id` values can appear when different `site_key` rows preserve
different protein context. Interpret rows by `site_key`.

The workflow does not run preprocessing, imputation, localisation filtering,
sequence resolution, batch correction, random-effects modelling, or
peptide-to-site differential aggregation. Fixed-effect batch terms are model
covariates, not ComBat, not RUV, not limma `removeBatchEffect` parity, and not
mixed-effects modelling.

## Common problems

| Problem | What to check |
| --- | --- |
| Missing sample error | Confirm every `SampleDesignRecord.sample_id` matches a dataset column, or deliberately set `allow_design_subset=True`. |
| Replicate error | Production contrasts need at least two biological replicates per contrasted condition. Technical replicates do not count. |
| Unsupported intensity scale | Build or declare a log2 dataset before differential analysis. |
| Low-confidence phosphosite failure | Add localisation metadata and use `DatasetLocalisationConfig` during dataset building. |
| Block or paired-design error | Use `paired_design_policy="fixed_block"` only when every block covers both sides of every requested contrast. |
| Imputation error | The default rejects upstream-imputed datasets. Use the withhold policy only with dataset-owned imputation observation metadata. |
| Empty or mostly missing result rows | Inspect `feature_eligibility`, caveats, and preprocessing reports before treating absence from a table as biological absence. |

## Related documentation

- [Preparing a dataset](dataset-build-workflow.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)
- [Scientific coverage](../scientific-coverage.md)
- [Kinase analysis](kinase.md)

# Differential Analysis Workflow

`DifferentialAnalysisWorkflow` helps you test explicit condition contrasts on an
`AnalysisReadyPhosphoDataset`. It keeps the scientific boundary clear: build and
preprocess the dataset first, then pass an explicit experimental design and
contrast list to the workflow.

## When to Use This Workflow

Use this workflow when you have a strict analysis-ready phosphosite matrix and
want moderated differential phosphorylation statistics for named condition
contrasts.

Current differential analysis is not full PhosR or limma parity. Supported
designs are limited to tested design and contrast envelopes documented in
[Scientific Coverage](../scientific-coverage.md) and protected by parity tests.
Column names such as `logFC`, `P.Value`, and `adj.P.Val` are familiar reporting
names, not a broad limma compatibility claim.

Good fits:

- two-condition unpaired contrasts
- explicit fixed-effect covariates such as batch, categorical covariates, or
  continuous covariates
- complete fixed-block designs when `paired_design_policy="fixed_block"` is set
  and every block covers the requested contrast

This workflow is not a dataset builder, peptide-to-site resolver, imputation
engine, broad batch-correction method, `duplicateCorrelation` implementation, or
mixed-effects model.
For PhosPy-origin peptide evidence, the supported peptide-to-site route is still
to resolve peptide intensities into analysis-ready site rows during dataset
building, then run this workflow. The former advanced post-hoc peptide
differential estimate-combination route has been withdrawn from public support
and its compatibility shell fails closed. Future public support requires
executable peptide-to-site mapping semantics and a coherent combined estimand and
inferential result.

## Inputs

`DifferentialAnalysisRequest.dataset` must be an
`AnalysisReadyPhosphoDataset`. The dataset must already have:

- numeric phosphosite-by-sample values in `dataset.phospho`
- rows keyed by `site_key`
- complete analysis-ready values at the dataset boundary
- an established log2 intensity scale for differential `logFC` reporting
- site metadata carrying `site_key`, `display_id`, `organism`,
  `protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`

`dataset.sample_metadata` may exist, but differential design semantics come from
`ExperimentalDesign`, not passive sample metadata.

Differential results report `logFC`, so build or declare a log2 dataset before
running this workflow. Linear or unestablished intensity-scale state is rejected
before model execution. For example, use
`DatasetIntensityTransformConfig(policy="log2")` during dataset building when
your input values are linear intensities.
Declared log2 state that carries suspicious declaration diagnostics is also
rejected by default. Use
`DifferentialAnalysisConfig(allow_suspicious_declared_input_scale=True)` only
when the declaration is scientifically trusted; this override is recorded in
differential policy provenance.

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
| `allow_suspicious_declared_input_scale` | `False` | Explicit override for declared log2 input-scale provenance that recorded suspicious diagnostics. |
| `minimum_condition_replicates` | `2` | Minimum replicate count per contrast side. |
| `empirical_bayes` | `EmpiricalBayesConfig()` | Controls `method`, `trend`, and robust winsor tails. |
| `multiple_testing` | `MultipleTestingConfig()` | Controls how raw p-values become `adj.P.Val`. Default method is `"benjamini_hochberg"`. |

Related request classes:

- `ExperimentalDesign`
- `SampleDesignRecord`
- `Contrast`
- `BatchCovariate`
- `CategoricalCovariate`
- `ContinuousCovariate`
- `EmpiricalBayesConfig`
- `MultipleTestingConfig`

## Multiple-Testing Correction

`MultipleTestingConfig(method=...)` accepts:

- `"benjamini_hochberg"` (default)
- `"bonferroni"`
- `"holm"`
- `"benjamini_yekutieli"`
- `"none"`

Correction is applied separately for each contrast table. PhosPy does not pool
p-values across different contrasts. If an imputation withhold policy removes
features from testing, correction uses only the tested features for that
contrast.

`adj.P.Val` is the p-value after the selected correction method. Smaller values
indicate stronger evidence after accounting for the number of tested features in
that contrast, but the adjusted value is not an effect size. Report the method
and thresholds used in your analysis.

## Peptide-to-Site Differential Evidence

Preferred PhosPy-origin lane:

1. Build the dataset with `site_resolution_mode="peptide_evidence"`,
   `peptide_evidence_sample_intensity_columns=...`, and an explicit
   `multi_site_policy`.
2. Let dataset building resolve peptide evidence into site-level sample
   intensities with provenance.
3. Run `DifferentialAnalysisWorkflow` on the resulting
   `AnalysisReadyPhosphoDataset`.

The post-hoc peptide differential estimate-combination lane is withdrawn from
public support. The retained compatibility shell is internal/experimental and
fails closed with an error explaining that coherent combined effect/inference
semantics and executable peptide-to-site mapping semantics are not implemented.
It must not silently execute mapping policies such as equal splitting or
statistical-model exclusion as ordinary evidence.

Future public support requires a new scientific design decision and executable
implementation that defines:

- peptide-to-site mapping semantics, including ambiguous and excluded evidence;
- the combined estimand for site-level effects;
- the inferential result and uncertainty model;
- dependence handling for same-experiment peptide estimates;
- multiple-testing and provenance semantics.

Until then, resolve peptide evidence at sample-intensity level during dataset
building and then run the core `DifferentialAnalysisWorkflow`.

## Contrast Helpers

Manual `Contrast` objects are still the most explicit option. For common
condition comparisons, PhosPy also provides helpers that return ordinary
`Contrast` tuples. They do not run the workflow or add statistical behaviour.

Contrast direction is `numerator_condition - denominator_condition`. A contrast
named `treatment_vs_control` means treatment is the numerator and control is the
denominator.

All pairwise contrasts use the condition order from the design. For each pair,
the later condition is compared with the earlier condition:

```python
from phospy.api import all_pairwise_contrasts

contrasts = all_pairwise_contrasts(design)

assert tuple(contrast.name for contrast in contrasts) == (
    "B_vs_A",
    "C_vs_A",
    "C_vs_B",
)
```

Treatment-versus-control contrasts use each non-control condition as the
numerator and the named control as the denominator:

```python
from phospy.api import contrasts_vs_control

contrasts = contrasts_vs_control(
    design,
    control_condition="control",
)

assert tuple(contrast.name for contrast in contrasts) == (
    "treatment_a_vs_control",
    "treatment_b_vs_control",
)
```

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
| `diagnostics` | Explicit model scope, design, contrast, sample/site count, method, imputation, scale, normalisation, unsupported-assumption, and warning diagnostics. |
| `policy_provenance` | Structured design, contrast, replicate, imputation, input-scale, and testing provenance. |
| `workflow_provenance` | Workflow-level execution metadata. |
| `input_dataset_preprocessing_report` | Preprocessing report carried from the input dataset when available. |

Each contrast result table is indexed by the input `site_key`. The minimum public
identity columns are `site_key`, `display_id`, `organism`,
`protein_namespace`, `protein_identifier`, `gene_symbol`, and `site`. Typical
statistics columns include `logFC`, `t`, `P.Value`, and `adj.P.Val`.

stat-only computation payloads are internal. They are not a public scientific
result object and are not valid `DifferentialAnalysisResult` tables.

## Filtering and Ranking Result Tables

Use table helpers when you want a smaller reporting table after the workflow has
finished. These helpers only filter or sort existing columns. They do not refit
the model, recompute p-values, or change the `DifferentialAnalysisResult`.

Filter by adjusted p-value and absolute `logFC`:

```python
from phospy.api import filter_differential_results

table = result.table_for("treatment_vs_control")

reported = filter_differential_results(
    table,
    adjusted_p_value_max=0.05,
    min_abs_effect_size=1.0,
)
```

This example reports rows with `adj.P.Val <= 0.05` under the configured
multiple-testing method and `abs(logFC) >= 1.0`. Choose and state thresholds as
part of your reporting; they are not model settings.

Rank by raw p-value, or by absolute `logFC` for the largest fitted effects:

```python
from phospy.api import rank_differential_results

ranked_by_p = rank_differential_results(
    table,
    by="P.Value",
)

ranked_by_effect = rank_differential_results(
    table,
    by="logFC",
    ascending=False,
    absolute=True,
)
```

Missing requested columns raise a clear input error.

## Interpreting the Result

`logFC` is the fitted condition contrast on the established log2 phosphosite
intensity scale. `t` is a moderated t-statistic. `P.Value` is the raw p-value,
and `adj.P.Val` is the multiple-testing adjusted value for the configured
method.

PhosPy uses its own moderated OLS-style implementation. The result table follows
familiar limma/topTable-style column names (`logFC`, `P.Value`, `adj.P.Val`), but
that wording does not mean exact limma or `topTable` numerical parity.
`result.policy_provenance.statistical_testing.input_intensity_scale` records the
validated input scale, and `logfc_interpretation` records how `logFC` should be
read.

Repeated `display_id` values can appear when different `site_key` rows preserve
different protein context. Interpret rows by `site_key`; `display_id` is for
readability.

## Provenance and Reproducibility

Result diagnostics and provenance record the resolved design, contrast vectors,
fixed-effect covariates, replicate policy, empirical-Bayes settings,
multiple-testing method, imputation policy, normalisation state, unsupported
assumptions, and unsupported-design rejection policy. Warnings are exposed on
`result.diagnostics.warnings`; they are not only logged. Table exports return
defensive in-memory snapshots; mutating them does not mutate the result object.

## Limitations

- Conditions, batches, blocks, and covariates are not inferred from sample names.
- The workflow does not run preprocessing, imputation, localisation filtering,
  sequence resolution, or batch correction.
- Fixed-effect batch terms are model covariates, not ComBat, RUV,
  `removeBatchEffect`, `duplicateCorrelation`, or mixed-effects modelling.
- Fixed-effect covariates are not full batch correction; they are ordinary
  design terms in the fitted model.
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

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

## Empirical Bayes Variance Moderation

Differential analysis estimates one residual variance per feature, then uses
empirical Bayes moderation to borrow variance information across features
before computing moderated *t* statistics and *p* values. This stabilizes
inference when each site has limited replicate information.

PhosPy supports three empirical-Bayes prior modes:

- `EmpiricalBayesConfig()` uses the existing global prior. All fitted features
  contribute to a shared prior-variance distribution, and no trend diagnostics
  are emitted.
- `EmpiricalBayesConfig(trend=True)` keeps the existing mean-intensity trend
  mode. The prior variance is fitted as a population-level function of the
  feature mean intensity, and `mean_variance_trend_diagnostics` is populated.
- `EmpiricalBayesConfig(trend=True, trend_covariate="quantification_depth",
  quantification_depth_kind=...)` uses the existing trend machinery with
  quantification depth as the variance-trend covariate, and
  `quantification_depth_trend_diagnostics` is populated.

Quantification depth can matter in proteomics and phosphoproteomics because
measurement precision is often related to how much evidence supports a
feature. A phosphosite quantified from many peptide-spectrum matches or unique
peptides may belong to a different variance population from a site supported
by sparse evidence. The fitted depth trend estimates that population-level
relationship; it is not a guarantee that every site with a higher PSM or
peptide count has lower variance than every site with lower depth.

Depth-aware moderation is opt-in and requires two pieces of information:

- `dataset.site_metadata["quantification_depth"]`: one feature-aligned numeric
  count per site that will be tested. Counts must be finite, integer-valued,
  and at least one.
- `EmpiricalBayesConfig.quantification_depth_kind`: the declared meaning of
  the counts. Use `"psm_count"` for peptide-spectrum matches supporting the
  feature, and `"peptide_count"` for unique peptides or peptide forms
  supporting the feature.

PhosPy validates and stores the raw count, then fits the variance trend against
`log2(quantification_depth)`. The diagnostics therefore expose both
`quantification_depth` and the transformed `log2_quantification_depth` trend
covariate.

### Manual Quantification-Depth Metadata Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DifferentialAnalysisConfig,
    EmpiricalBayesConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

phospho = pd.DataFrame(
    {
        "control_1": [1000.0, 900.0, 700.0],
        "control_2": [1040.0, 880.0, 720.0],
        "treated_1": [1800.0, 930.0, 760.0],
        "treated_2": [1760.0, 920.0, 750.0],
    },
    index=["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
        "site": ["Y182", "T308", "S9"],
        "protein_id": ["P53778", "P31749", "P49841"],
        "site_sequence": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        ],
        "localisation_confidence": [0.95, 0.96, 0.97],
        "quantification_depth": [12, 4, 7],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.HUMAN,
        input_intensity_scale="linear",
        preprocessing_config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2")
        ),
    )
)

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord("control_1", "control", "control_r1"),
        SampleDesignRecord("control_2", "control", "control_r2"),
        SampleDesignRecord("treated_1", "treated", "treated_r1"),
        SampleDesignRecord("treated_2", "treated", "treated_r2"),
    )
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=(
            Contrast(
                name="treated_vs_control",
                numerator_condition="treated",
                denominator_condition="control",
            ),
        ),
        config=DifferentialAnalysisConfig(
            empirical_bayes=EmpiricalBayesConfig(
                trend=True,
                trend_covariate="quantification_depth",
                quantification_depth_kind="psm_count",
            )
        ),
    )
)

depth_trend = result.quantification_depth_trend_diagnostics
assert depth_trend is not None

print(result.policy_provenance.empirical_bayes.trend_covariate)
print(result.policy_provenance.empirical_bayes.quantification_depth_kind)
print(depth_trend.quantification_depth_series())
print(depth_trend.trend_covariate_series())
print(result.diagnostics.moderation_method)
```

MaxQuant and FragPipe/PTMProphet importers can populate
`site_metadata["quantification_depth"]` only when the caller maps an explicit
source column and declares its `quantification_depth_kind`. They do not infer
PSM or peptide depth from arbitrary numeric columns, and depth is left
unavailable when multi-site peptide evidence is split in a way that would
require inventing per-site counts.

Depth-aware moderation is validated against a checked-in R/DEqMS
`spectraCounteBayes` reference fixture within documented tolerances. Because
PhosPy preserves its existing deterministic trend smoother while DEqMS uses
its own loess implementation, describe the feature as
quantification-depth-aware empirical Bayes moderation inspired by DEqMS, not
exact DEqMS-compatible numerical equivalence.

## Complete Public-API Example: Paired Duplicate-Correlation Design

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DifferentialAnalysisConfig,
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
)
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
        "control_rep1": [1000.0, 900.0, 800.0, 700.0],
        "control_rep2": [1050.0, 880.0, 820.0, 710.0],
        "control_rep3": [980.0, 910.0, 790.0, 690.0],
        "treatment_rep1": [1800.0, 930.0, 760.0, 740.0],
        "treatment_rep2": [1750.0, 920.0, 740.0, 730.0],
        "treatment_rep3": [1825.0, 940.0, 770.0, 760.0],
    },
    index=["MAPK14;Y182;", "GSK3A;S21;", "TSC2;S939;", "GSK3B;S9;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "GSK3A", "TSC2", "GSK3B"],
        "site": ["Y182", "S21", "S939", "S9"],
        "site_sequence": [
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
        ],
        "protein_namespace": ["protein_id", "protein_id", "protein_id", "protein_id"],
        "protein_identifier": ["MAPK14", "GSK3A", "TSC2", "GSK3B"],
        "localisation_confidence": [0.95, 0.94, 0.96, 0.93],
    },
    index=phospho.index,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale="linear",
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
        SampleDesignRecord(
            "control_rep1", "control", "control_r1", block_id="subject_1"
        ),
        SampleDesignRecord(
            "treatment_rep1", "treatment", "treatment_r1", block_id="subject_1"
        ),
        SampleDesignRecord(
            "control_rep2", "control", "control_r2", block_id="subject_2"
        ),
        SampleDesignRecord(
            "treatment_rep2", "treatment", "treatment_r2", block_id="subject_2"
        ),
        SampleDesignRecord(
            "control_rep3", "control", "control_r3", block_id="subject_3"
        ),
        SampleDesignRecord(
            "treatment_rep3", "treatment", "treatment_r3", block_id="subject_3"
        ),
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
        config=DifferentialAnalysisConfig(
            paired_design_policy=PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
        ),
    )
)

duplicate = result.policy_provenance.duplicate_correlation
assert duplicate is not None

print(
    result.table_for("treatment_vs_control").loc[
        :, ["display_id", "logFC", "P.Value", "adj.P.Val"]
    ]
)
print("consensus correlation:", duplicate.consensus.consensus_correlation)
print("estimator features:", duplicate.consensus.estimated_feature_count)
```

## Experimental Protein-Aware Differential Analysis

`protein_covariate_adjusted_moderated_linear_model_v1` is an experimental,
opt-in differential estimator. It reports the requested phosphosite condition
contrast conditional on measured matched total-protein abundance under
`y_s = X beta_s + z_p(s) gamma_s + error`. It requires established log2
phosphosite and total-protein inputs; centres, but does not standardize,
normalize, subtract, or impute, the total-protein covariate; withholds
inadmissible sites rather than falling back to the ordinary estimator; and is
not stoichiometry, occupancy, causal decomposition, MSstatsPTM parity, or
MSstatsPTM-style joint PTM/protein inference.

Use this lane when the dataset already carries dataset-owned protein-aware
preparation from
`DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")` and the
scientific question calls for a condition effect adjusted for the measured
matched total-protein abundance. The workflow request still has only
`dataset`, `design`, `contrasts`, and `config`; do not pass a protein matrix,
mapping table, or preparation result to `DifferentialAnalysisRequest`.
Differential analysis does not rerun protein mapping or automatically
preprocess total-protein values.

The ordinary differential lane remains the default. A dataset may carry
`dataset.protein_aware_preparation`, but `DifferentialAnalysisWorkflow` ignores
that sidecar unless `DifferentialAnalysisConfig.protein_aware_model` is set.

### Complete Protein-Aware Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.advanced import (
    DatasetProteinAwarePreparationConfig,
    DifferentialAnalysisConfig,
    DifferentialProteinAwareModelConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)

phospho = pd.DataFrame(
    {
        "A_1": [1.00, 2.05, 1.48],
        "A_2": [1.15, 2.10, 1.50],
        "A_3": [0.95, 1.92, 1.46],
        "B_1": [1.75, 2.48, 1.55],
        "B_2": [1.83, 2.57, 1.58],
        "B_3": [1.69, 2.41, 1.62],
    },
    index=["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
        "site": ["Y182", "T308", "S9"],
        "protein_id": ["P53778", "P31749", "P49841"],
        "site_sequence": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        ],
        "localisation_confidence": [0.95, 0.96, 0.97],
    },
    index=phospho.index.copy(),
)

total = pd.DataFrame(
    {
        "A_1": [10.0, 8.0, 12.5],
        "A_2": [10.4, 8.5, 12.5],
        "A_3": [9.7, 7.8, 12.5],
        "B_1": [11.8, 8.9, 12.5],
        "B_2": [12.2, 9.7, 12.5],
        "B_3": [11.4, 8.8, 12.5],
    },
    index=pd.Index(["P53778", "P31749", "P49841"], name="protein_id"),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        organism=Organism.RAT,
        input_intensity_scale="log2",
        preprocessing_config=DatasetPreprocessingConfig(
            protein_aware_preparation=DatasetProteinAwarePreparationConfig(
                policy="prepare_model_inputs"
            )
        ),
    )
)

design = ExperimentalDesign(
    samples=(
        SampleDesignRecord(
            sample_id="A_1",
            condition="A",
            biological_replicate_id="A_1_bio",
        ),
        SampleDesignRecord(
            sample_id="A_2",
            condition="A",
            biological_replicate_id="A_2_bio",
        ),
        SampleDesignRecord(
            sample_id="A_3",
            condition="A",
            biological_replicate_id="A_3_bio",
        ),
        SampleDesignRecord(
            sample_id="B_1",
            condition="B",
            biological_replicate_id="B_1_bio",
        ),
        SampleDesignRecord(
            sample_id="B_2",
            condition="B",
            biological_replicate_id="B_2_bio",
        ),
        SampleDesignRecord(
            sample_id="B_3",
            condition="B",
            biological_replicate_id="B_3_bio",
        ),
    )
)

result = DifferentialAnalysisWorkflow().run(
    DifferentialAnalysisRequest(
        dataset=dataset,
        design=design,
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
        config=DifferentialAnalysisConfig(
            protein_aware_model=DifferentialProteinAwareModelConfig(
                method="protein_covariate_adjusted_moderated_linear_model_v1"
            )
        ),
    )
)

table = result.table_for("B_vs_A")
tested = table.loc[table["result_status"] == "tested"]
withheld = table.loc[table["result_status"] != "tested"]
diagnostics = result.protein_aware_diagnostics
assert diagnostics is not None

print(tested.loc[:, ["display_id", "logFC", "P.Value", "adj.P.Val"]])
print(withheld.loc[:, ["display_id", "result_status", "result_status_reason"]])
print("method:", diagnostics.method_id)
print("claim:", diagnostics.claim_status)
print("tested sites:", diagnostics.tested_site_count)
print("withheld sites:", diagnostics.withheld_site_count)
print("fallback policy:", diagnostics.fallback_policy)
print(
    diagnostics.per_site_diagnostics_dataframe().loc[
        :, ["total_protein_row_key", "protein_covariate_coefficient"]
    ]
)
print([caveat.code for caveat in result.caveats])
```

### Protein-Aware Compatibility

| Capability | Protein-Aware Version 1 |
| --- | --- |
| Fixed condition design | Supported when each augmented design is admissible. |
| Declared fixed covariates | Supported when each augmented design is admissible. |
| `fixed_block` | Supported when each augmented design is admissible. |
| `duplicate_correlation` | Rejected before fitting. |
| No-op technical-replicate policy | Supported. |
| Actual technical-replicate aggregation | Rejected before fitting. |
| `allow_design_subset=True` | Supported by exact deterministic protein-covariate subsetting and reordering. |
| Existing phosphosite imputation policy | Preserved for phosphosite eligibility. |
| Protein-covariate imputation | Not performed; non-finite covariates are withheld. |
| Prior `subtract_log_total` | Rejected before fitting. |
| Protein-aware fallback to ordinary fitting | Not allowed. Withheld rows stay visible in the full result index. |
| Mixed effects | Not supported in version 1. |

### Protein-Aware Results

Each contrast table keeps the full `site_key` index. Tested rows contain the
usual `logFC`, `t`, `P.Value`, and `adj.P.Val` values. Withheld rows keep their
identity columns and receive missing statistics with typed `result_status` and
`result_status_reason` values such as:

- `withheld_protein_preparation_ineligible`
- `withheld_protein_covariate_invalid`
- `withheld_protein_augmented_design_invalid`
- `withheld_protein_contrast_non_estimable`

The protein covariate coefficient is reported only as a nuisance coefficient in
diagnostics. It is not a separate inferential result and has no public
protein-covariate *p* value or adjusted *p* value claim.

`result.protein_aware_diagnostics` records the selected method and
`experimental` claim status, preparation and mapping policies, execution sample
order, centring/no-standardization policy, no-imputation and no-fallback
policies, tested/withheld counts, reason counts, base-design diagnostics,
augmented-design condition-number summaries, and full per-site diagnostics.
Use `result.protein_aware_site_diagnostics_dataframe()` or
`result.protein_aware_diagnostics.per_site_diagnostics_dataframe()` for an
independent diagnostics table.

`result.policy_provenance.protein_aware` records the method version, model
formula, conditional `logFC` interpretation, transformation and prior
total-protein-correction state, design-subset behavior, duplicate-correlation
and technical-aggregation rejection policies, order-sensitive input
fingerprints for the phosphosite matrix, matched pairs, protein-covariate
matrix, eligibility table, design matrix, and contrast matrix, plus limitations
and unsupported claims. Scalar diagnostics on `result.diagnostics` describe the
base design; augmented-design summaries are aggregates across successfully
fitted total-protein-row groups.

## Request

Create a `DifferentialAnalysisRequest`.

| Parameter | Type | Required or Default | Description | Main Constraint |
| --- | --- | --- | --- | --- |
| `dataset` | `AnalysisReadyPhosphoDataset` | Required | Dataset to test. | Must be complete, `site_key` indexed, and established as log2. |
| `design` | `ExperimentalDesign` | Required | Explicit sample and model design. | Sample IDs must match dataset columns unless subsetting is enabled. |
| `contrasts` | `tuple[Contrast, ...]` | Required | Comparisons to estimate. | Names must be unique; numerator and denominator conditions must exist and differ. |
| `config` | `DifferentialAnalysisConfig` | `DifferentialAnalysisConfig()` | Reliability, replicate, imputation, moderation, multiple-testing, and optional protein-aware model policy. | Unsupported combinations fail before fitting. |

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

Paired and repeated-sample handling is controlled by
`DifferentialAnalysisConfig.paired_design_policy`; the supported values,
constructed designs, limitations, and provenance fields are documented below.

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
| `protein_aware_model` | `DifferentialProteinAwareModelConfig` or `None` | `None` | Selects the experimental protein-covariate-adjusted lane. When absent, ordinary differential analysis ignores any dataset-owned protein-aware preparation sidecar. |

### `DifferentialProteinAwareModelConfig`

| Parameter | Default | Supported Values |
| --- | --- | --- |
| `method` | `"protein_covariate_adjusted_moderated_linear_model_v1"` | `"protein_covariate_adjusted_moderated_linear_model_v1"` |

### `paired_design_policy`

Set this on the request config:
`DifferentialAnalysisRequest(config=DifferentialAnalysisConfig(paired_design_policy=...))`.
The default remains `"reject"`. PhosPy does not inspect the design and choose a
blocked model automatically.

| Value | Use When | Design Constructed | Main Failure Behaviour |
| --- | --- | --- | --- |
| `"reject"` | Use for ordinary unblocked designs, or when paired/repeated `block_id` metadata must fail unless the caller deliberately selects a supported blocked model. | Condition terms and supported fixed covariates only. No block fixed effects and no block correlation structure are created. | If analysed samples carry `block_id` metadata, validation fails before fitting. |
| `"fixed_block"` | Use for complete paired or blocked designs where block identity should be adjusted as an ordinary fixed nuisance effect. | Condition terms, supported fixed covariates, and block dummy columns are included in the fixed-effects design. One block level is the reference. | Every block must cover both sides of every requested contrast; incomplete contrast coverage, rank deficiency, or non-estimable contrasts fail before fitting. |
| `"duplicate_correlation"` | Use for repeated, multi-condition, or incomplete blocked designs where one shared within-block correlation is appropriate and fixed block coefficients are not desired. | Condition terms and supported fixed covariates are included in the fixed-effects design. `block_id` is supplied separately as a correlation group and is not added as fixed block coefficients. | Missing block IDs, no repeated blocks, rank-deficient non-block design, insufficient residual degrees of freedom, unsupported weights, failed consensus estimation, or failed final GLS fitting fail without fallback. |

`fixed_block` remains a valid supported policy. Block identity is included as a
fixed nuisance effect, and the condition estimate is based on within-block
comparisons. For a complete two-condition paired experiment, it is the
linear-model analogue of a paired test. It does not estimate within-block
correlation.

`duplicate_correlation` supplies block identity as a correlation group rather
than a fixed coefficient. PhosPy estimates feature-wise correlations by REML,
combines eligible estimates into one consensus, and refits all features by GLS
using that consensus. The model assumes one shared compound-symmetry
correlation. It can be useful for repeated, multi-condition, or incomplete
blocked designs, but it is not a general mixed-effects model and is not
universally superior to `fixed_block`.

The implemented estimator is the residual-space variance-component REML
formulation recorded in ADR-0048. For each feature, finite observations are
selected, the fixed-effects design is removed through a residual-space basis,
and residual/block variance components are fitted before mapping the block
component ratio to a raw correlation. Feature estimates are then clamped before
Fisher aggregation: for maximum observed repeated-block size `m >= 2`, the
feature-level lower clamp is `-1 / (m - 1) + 0.01` and the upper clamp is
`0.99`. Boundary-clamped finite estimates can contribute to the consensus. The
workflow-level consensus is checked again against the full block structure
before GLS.

!!! warning "First duplicate-correlation implementation limits"
    The first implementation does not support multiple random effects, random
    slopes, feature-specific final covariance models, arbitrary longitudinal
    covariance, time-dependent correlation, nested or crossed random-effects
    syntax, simultaneous fixed block coefficients and duplicate correlation,
    automatic policy selection, user-supplied consensus correlation,
    user-configurable trimming, or unsupported precision weights.

    Do not add block dummy variables and a block-correlation structure to the
    same model. Select `fixed_block` or `duplicate_correlation` explicitly.

### Missing, Imputed, and Authoritative Matrix Policy

Differential fitting consumes the workflow-approved analysis matrix. Actual
missing values (`NaN`) are rejected before model fitting. By default,
upstream-imputed datasets are also rejected.

If `imputed_value_policy="withhold_imputed_features"` is selected, PhosPy uses
dataset-owned observation metadata to decide which features remain `tested`.
Withheld rows do not contribute to duplicate-correlation consensus estimation,
GLS fitting, or multiple-testing adjustment. Tested rows are fitted on the
workflow-approved matrix. If tested rows contain retained imputed cells, those
cells participate in REML and GLS, and the result records this through
imputation provenance and caveats. This is not observed-only fitting, does not
use feature-specific residual degrees of freedom, and does not change the
duplicate-correlation covariance model.

For `duplicate_correlation`, provenance records the
`analysis_matrix_fingerprint` and `authoritative_matrix_fingerprint`. They
match by design: the approved fitting matrix is the matrix authority for both
REML consensus estimation and final GLS.

### `EmpiricalBayesConfig`

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `method` | `"standard"` or `"robust"` | `"standard"` | Moderation method. |
| `trend` | `bool` | `False` | Enables feature-specific prior variances through a variance trend. |
| `trend_covariate` | `"mean_intensity"` or `"quantification_depth"` | `"mean_intensity"` | Selects the variance-trend covariate when `trend=True`. The default preserves the existing mean-intensity trend mode. |
| `quantification_depth_kind` | `"psm_count"`, `"peptide_count"`, or `None` | `None` | Required only when `trend_covariate="quantification_depth"`; declares what the depth counts mean. |
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

### Failure Behaviour

PhosPy fails closed for unsupported blocked designs.

- `paired_design_policy="reject"` rejects paired or repeated designs when the
  caller does not want PhosPy to select a supported blocked model.
- `paired_design_policy="fixed_block"` rejects missing block IDs, incomplete
  within-block contrast coverage, rank-deficient designs, non-estimable
  contrasts, and insufficient residual degrees of freedom before fitting.
- `paired_design_policy="duplicate_correlation"` rejects missing block IDs, no
  repeated blocks, fixed block coefficients in the duplicate-correlation
  design, unsupported precision weights, rank-deficient non-block designs,
  non-estimable contrasts, and designs with two or fewer residual degrees of
  freedom for REML. If REML consensus estimation or final GLS fitting fails,
  the workflow reports a fitting failure.

There is no silent fallback from `duplicate_correlation` to `fixed_block`,
ordinary least squares, or correlation zero.

The protein-aware lane also fails closed before fitting when the selected
method lacks a current dataset-owned preparation sidecar, the sidecar is
unsupported or cannot be proven to belong to the current dataset, the dataset
has no total-protein matrix, sample identities are inconsistent, phosphosite or
total-protein values are not established log2, prior `subtract_log_total` was
applied, actual technical-replicate aggregation is required,
`paired_design_policy="duplicate_correlation"` is selected, or every site is
withheld. Row-level protein mapping, protein-covariate numeric, augmented-rank,
conditioning, residual-DF, and contrast-estimability failures are reported as
typed withheld rows when at least one site remains testable.

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
| `mean_variance_trend_diagnostics` | `MeanVarianceTrendDiagnostics` or `None` | Legacy mean-intensity trend diagnostics when mean-intensity trend moderation is enabled. |
| `quantification_depth_trend_diagnostics` | `QuantificationDepthTrendDiagnostics` or `None` | Depth-aware trend diagnostics when quantification-depth trend moderation is enabled; includes original depth values, log2 trend covariate, depth kind, residual log variance, and fitted log prior variance. |
| `diagnostics` | `DifferentialModelDiagnostics` | Design, contrast, scale, and model diagnostics. |
| `protein_aware_diagnostics` | `ProteinAwareDifferentialDiagnostics` or `None` | Experimental protein-aware method, counts, group summaries, and per-site diagnostics when the opt-in lane is selected. |
| `protein_aware_site_diagnostics_dataframe()` | `pandas.DataFrame` or `None` | Independent snapshot of per-site protein-aware diagnostics. |
| `policy_provenance`, `workflow_provenance` | Typed or mapping-like provenance | Resolved scientific policy and execution metadata. |
| `caveats` | `tuple[ResultCaveat, ...]` | Structured interpretation limits. |
| `input_dataset_preprocessing_report` | Report or `None` | Dataset preprocessing report carried into the result. |
| `to_payload()` | JSON-compatible mapping | Serializable tables, diagnostics, provenance, and caveats. |
| `scientifically_equals(...)` | `bool` | Scientific-result comparison helper. |

### Policy and Duplicate-Correlation Provenance

`result.policy_provenance` records the resolved design and fitting policy. The
design summary includes `paired_design_policy`, `block_id_field_name`,
`block_count`, `block_levels`, `block_column_names`,
`condition_coverage_rule`, validation statuses, and limitations. For
`fixed_block`, `block_column_names` names the fixed nuisance columns used in the
linear model. For `duplicate_correlation`, `block_column_names` is empty because
the block IDs are supplied as correlation groups instead.

`result.policy_provenance.empirical_bayes` records the empirical-Bayes method,
robust status, winsor settings, whether variance trending occurred, the trend
covariate, and the covariate transformation. Global-prior fits record no trend
covariate. Mean-intensity trend fits record `trend_covariate="mean_intensity"`
with `trend_covariate_transformation="identity"`. Depth-aware trend fits record
`trend_covariate="quantification_depth"`,
`trend_covariate_transformation="log2"`, and the selected
`quantification_depth_kind` (`psm_count` or `peptide_count`).

Depth-aware variance moderation is opt-in. The workflow reads
`site_metadata["quantification_depth"]` as a feature-aligned count and requires
finite integer-valued values greater than or equal to one for every feature
that enters differential testing and empirical-Bayes moderation. Features
withheld before testing do not require valid depth and are not included when
fitting the depth-aware trend. They may still appear in re-expanded contrast
tables with missing or NaN quantification-depth trend diagnostics. Use
`quantification_depth_kind="psm_count"` when the depth is the number of
peptide-spectrum matches supporting the feature, and
`quantification_depth_kind="peptide_count"` when it is the number of unique
peptides or peptide forms supporting the feature. PhosPy does not infer these
counts from arbitrary metadata columns; callers must provide the exact
feature-level evidence count they want treated as quantification depth.

The depth-aware path fits the existing empirical-Bayes trend estimator with
`log2(quantification_depth)` as the trend covariate. It is a
quantification-depth-aware moderation strategy, but this release does not claim
DEqMS numerical compatibility. The fitted trend describes a population-level
mean relationship between depth and residual variance; it should not be read as
a guarantee that every higher-depth feature has lower variance than every
lower-depth feature.

`result.prior_residual_variance`, `result.prior_diagnostics.prior_variance`,
and `result.to_payload()["empirical_bayes"]["prior_residual_variance_by_feature"]`
expose the feature-level prior variances used by moderation. Trend diagnostics
are split by covariate: mean-intensity analyses populate
`mean_variance_trend_diagnostics`, depth-aware analyses populate
`quantification_depth_trend_diagnostics`, and non-trend analyses populate
neither trend diagnostics payload.

When `paired_design_policy="duplicate_correlation"`,
`result.policy_provenance.duplicate_correlation` records:

- requested and normalised paired-design policy;
- `block_treatment="consensus_correlation"` and
  `covariance_structure="compound_symmetry"`;
- estimator name, estimator policy version, and fixed trim fraction;
- matrix, design, and block-assignment fingerprints;
- `block_structure`, including sample count, block count, repeated-block count,
  singleton-block count, correlated-pair count, block levels, and minimum and
  maximum block size;
- `consensus.consensus_correlation`, eligible feature count, estimated feature
  count, failed feature count, non-finite feature count, and the consensus
  failure-reason field, which is normally empty in a returned successful result;
- attempted, trimmed, and retained estimator-feature counts;
- `failure_reason_counts`, `convergence_summary`, and `boundary_summary`;
- `gls_fit_status`; and
- `imputed_values_participated`, `imputed_feature_count`, and
  `imputed_cell_count`.

`result.workflow_provenance` also carries a serializable
`duplicate_correlation` entry when the policy is selected. `result.caveats`
contains an informational `differential_duplicate_correlation_consensus` caveat
stating that one consensus compound-symmetry correlation was used and that
feature-specific random effects were not fitted.

When `protein_aware_model` is selected,
`result.policy_provenance.protein_aware` records the experimental estimator
policy and exact input fingerprints. `result.workflow_provenance` also carries
protein-aware row-attrition metrics, and `result.caveats` records the
experimental interpretation limits and unsupported claims.

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

For `protein_covariate_adjusted_moderated_linear_model_v1`, `logFC` is the
requested phosphosite condition contrast conditional on the matched measured
total-protein covariate. It is not a phosphosite/protein subtraction,
stoichiometry, occupancy, or a causal decomposition of protein abundance and
phosphorylation regulation. Adjustment can change the estimand when protein
abundance is part of the biological pathway, and it does not establish causal
independence.

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

For `fixed_block`, interpret the condition contrast as an estimate adjusted for
ordinary block fixed effects. The estimate is driven by within-block
comparisons where the requested contrast is covered by every block. This policy
is valid and supported; it is not an unpaired analysis and is not deprecated by
`duplicate_correlation`.

For `duplicate_correlation`, interpret
`policy_provenance.duplicate_correlation.consensus.consensus_correlation` as
the one within-block compound-symmetry correlation used for all final GLS
feature fits. Positive values increase the modelled covariance among samples in
the same block; negative values decrease it within the valid
compound-symmetry bounds for the observed block sizes. The value is a modelling
parameter for this workflow, not a per-feature biological correlation.

A feature can fail to contribute to the REML consensus and still receive final
GLS statistics if the workflow obtains a valid consensus from other eligible
features. Those estimator failures are summarized by
`consensus.failed_feature_count`, `consensus.non_finite_feature_count`, and
`failure_reason_counts`. A final GLS or moderated-statistic failure is
different: it means the tested feature could not be fitted under the selected
final model. Under `duplicate_correlation`, such a final fitting failure stops
the workflow rather than being treated as only a consensus-contribution
failure.

Review `convergence_summary` and `boundary_summary` when interpreting a
duplicate-correlation result. They report numerical optimisation failures,
non-finite objectives or estimates, and estimates that converged at the valid
compound-symmetry boundary. These summaries are caveats for the consensus
estimate, not additional public estimator APIs.

If `imputed_values_participated` is true, at least one tested row contained
retained imputed cells under `imputed_value_policy="withhold_imputed_features"`.
Those values participated in REML consensus estimation and GLS fitting because
the workflow-approved matrix is authoritative. Review the imputation caveat and
feature-eligibility table before interpreting those rows.

`duplicate_correlation` is limited to one shared compound-symmetry correlation.
It is not multiple random effects, random slopes, feature-specific final
covariance, arbitrary longitudinal covariance, time-dependent correlation,
nested or crossed random-effects syntax, simultaneous fixed block coefficients
and duplicate correlation, automatic policy selection, user-supplied consensus
correlation, user-configurable trimming, or unsupported precision weights.

The experimental protein-aware lane is different from MSstatsPTM. MSstatsPTM
uses separate PTM-site and global-protein model outputs for adjusted PTM
inference. PhosPy version 1 instead fits one phosphosite-row model with a
matched total-protein abundance vector as a row-specific nuisance covariate.
Do not interpret the PhosPy result as MSstatsPTM parity, an MSstatsPTM-style
joint model, or a substitute for that workflow.

The committed differential limma parity fixtures demonstrate implementation
agreement only for the fixture-scoped model envelopes they cover. They are not
an independent scientific validation and do not make PhosPy generally identical
to limma.

For the duplicate-correlation fixture lane, fixtures A-C cover the complete
supported public path through final moderated statistics. Fixture D is narrower:
it checks internal feature-level REML and compound-symmetry GLS behavior for
controlled missingness/failure cases. It is not a public empirical-Bayes parity
fixture because public analysis-ready inputs with actual missing values fail
closed before fitting, and PhosPy does not return partial moderated output after
final duplicate-correlation GLS failures.

## Common Issues

| Issue | What to Check |
| --- | --- |
| A sample is missing from the dataset. | Match every `sample_id` to a dataset column, or deliberately enable `allow_design_subset`. |
| Replication fails. | Production contrasts need at least two biological replicates per condition. Technical replicates do not count. |
| The scale is unsupported. | Build or declare a log2 dataset before running the workflow. |
| Localization fails. | Add valid localization metadata and configure `DatasetLocalisationConfig` during dataset preparation. |
| A fixed-block design fails. | Confirm that every block covers both sides of every requested contrast. |
| A duplicate-correlation design fails. | Confirm every analysed sample has `block_id`, at least one block is repeated, the non-block fixed-effects design is full rank with more than two residual degrees of freedom, and no fixed block columns or unsupported weights are present. |
| Imputation is rejected. | Keep the default rejection policy, or use withholding only with builder-owned observation metadata. |
| Results contain withheld rows. | Review `feature_eligibility`, caveats, and the preprocessing report before interpreting absence biologically. |
| A protein-aware request reports a missing sidecar. | Build the dataset with `DatasetProteinAwarePreparationConfig(policy="prepare_model_inputs")`; do not pass the preparation result separately to the workflow request. |
| A protein-aware request reports incompatible scales. | Rebuild or transform both phosphosite and total-protein inputs so their log2 scale is established before differential analysis. |
| A protein-aware request reports mapping or preparation ineligibility. | Inspect `dataset.protein_aware_preparation`, the preparation report, and the row-level `result_status_reason` values. |
| A protein-aware request reports augmented-design rank or contrast failures. | Check fixed covariates, fixed blocks, sample subset order, and protein covariate variation for the affected total-protein-row groups. |
| A protein-aware request uses duplicate correlation or technical aggregation. | Use the ordinary lane for those policies, or use the protein-aware lane only with no actual technical aggregation and without `duplicate_correlation`. |

## Reference

Smyth, G. K., Michaud, J., & Scott, H. S. (2005). Use of within-array replicate
spots for assessing differential expression in microarray experiments.
*Bioinformatics, 21*(9), 2067–2075. https://doi.org/10.1093/bioinformatics/bti270

Kohler, D., Tsai, T.-H., Verschueren, E., Huang, T., Hinkle, T., Phu, L.,
Choi, M., & Vitek, O. (2023). MSstatsPTM: Statistical relative quantification
of posttranslational modifications in bottom-up mass spectrometry-based
proteomics. *Molecular & Cellular Proteomics, 22*(1), 100477.
https://doi.org/10.1016/j.mcpro.2022.100477

## Related Guides

- [Prepare a Dataset](dataset-build-workflow.md)
- [Scientific Interpretation and Limitations](../scientific-interpretation.md)
- [Scientific Coverage](../scientific-coverage.md)
- [Kinase Analysis](kinase.md)

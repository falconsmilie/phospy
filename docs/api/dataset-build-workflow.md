# Prepare a Dataset

Use `AnalysisReadyDatasetBuilder` to turn phosphosite intensities and metadata
into an `AnalysisReadyPhosphoDataset`. The builder checks row identity, sample
alignment, localisation, intensity scale, missing values, and any preprocessing
policy you select.

!!! info "At a Glance"
    **Input:** Phosphosite intensities, site metadata, and optional sample or
    total-protein tables  
    **Request:** `DatasetBuildRequest`  
    **Run:** `AnalysisReadyDatasetBuilder().run(request)`  
    **Returns:** An `AnalysisReadyPhosphoDataset` ready for supported workflows

For routine use, always go through the builder. Direct
`AnalysisReadyPhosphoDataset` construction raises immediately. The separate
advanced/trusted route,
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)`, is for callers that
already own fully validated, analysis-ready tables and complete typed evidence.
It is not a shortcut around the builder.

## Before You Begin

PhosPy expects a site-by-sample intensity matrix. Each row must resolve to one
phosphosite in a protein context. Human-readable labels such as `TSC2;S939;`
are useful for display, but they are not sufficient row identity when the same
label can refer to more than one protein or isoform.

The built dataset therefore uses:

- `site_key` as the unique, protein-scoped row index;
- `display_id` as the readable `GENE;SITE;` label; and
- `site_sequence` as required sequence context for analysis-ready rows.

Duplicate `display_id` values are allowed when their `site_key` values differ.
Duplicate `site_key` values require an explicit resolution policy because
merging them changes the represented evidence.

Site-level workflows also need localisation evidence. The default policy
requires a confidence value of at least 0.75 in the
`localisation_confidence` column.

## Example

This example starts with linear intensities, applies a log2 transform, and
requires site-level localisation evidence.

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.advanced import DatasetIntensityTransformConfig
from phospy.api import (
    DatasetBuildRequest,
    DatasetLocalisationConfig,
    DatasetPreprocessingConfig,
    Organism,
)

phospho = pd.DataFrame(
    {
        "control_1": [1000.0, 800.0],
        "control_2": [1050.0, 820.0],
        "treated_1": [1800.0, 760.0],
        "treated_2": [1750.0, 740.0],
    },
    index=["MAPK14;Y182;", "TSC2;S939;"],
)

site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "TSC2"],
        "site": ["Y182", "S939"],
        "site_sequence": [
            ("A" * 15) + "Y" + ("A" * 15),
            ("A" * 15) + "S" + ("A" * 15),
        ],
        "protein_identifier": ["MAPK14", "TSC2"],
        "localisation_confidence": [0.95, 0.96],
    },
    index=phospho.index,
)

request = DatasetBuildRequest(
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

dataset = AnalysisReadyDatasetBuilder().run(request)

print(dataset.phospho.shape)
print(dataset.intensity_scale_state.label)
print(dataset.site_metadata.loc[:, ["display_id", "site_sequence"]])
```

The builder may derive missing identity fields from safe input metadata, but the
finished dataset always has strict residue/position identity and the required
`AnalysisReadyPhosphoDataset` metadata.

## Input Tables

`DatasetBuildRequest` accepts pandas data frames or paths to `.csv`, `.tsv`,
`.txt`, and `.parquet` files. CSV, TSV, and TXT files use the first column as
the row index.

| Table | Shape | Purpose | Main Requirements |
| --- | --- | --- | --- |
| `phospho` | Sites × samples | Quantitative phosphosite matrix. | Numeric values; unique sample columns; row labels aligned with `site_metadata`. |
| `site_metadata` | Sites × fields | Site, sequence, protein, and localisation context. | Must align with phospho rows and provide enough information to derive strict `site_key` values. |
| `sample_metadata` | Samples × fields | Descriptive labels and preprocessing metadata. | Index must align with phospho columns; column names must be unique. |
| `total` | Proteins × samples | Optional total-protein matrix. | Sample columns must align with phospho columns. |
| `peptide_evidence` | Peptides × fields | Alternative peptide-evidence input lane. | Requires explicit intensity columns and multi-site handling. |

Text metadata is preserved as text: values such as `"NA"` and identifiers with
leading zeros are not silently converted. Quantitative tables must contain
numeric cells; invalid cells fail with row and column context.

`sample_metadata` does not define a differential design automatically. The
[differential workflow](differential-analysis.md) uses an explicit
`ExperimentalDesign` and `Contrast` objects.

## Request

Create a `DatasetBuildRequest`.

| Parameter | Type | Required or Default | Description | Main Constraint |
| --- | --- | --- | --- | --- |
| `phospho` | `pandas.DataFrame`, `str`, `Path`, or `None` | Required for `site_level_resolved` | Site-by-sample intensity matrix. | Must align with `site_metadata`. |
| `site_metadata` | `pandas.DataFrame`, `str`, `Path`, or `None` | Required for `site_level_resolved` | Site and protein metadata. | Must provide or safely support derivation of analysis-ready identity and `site_sequence`. |
| `sample_metadata` | `pandas.DataFrame`, `str`, `Path`, or `None` | `None` | Sample metadata used by selected preprocessing policies. | Index must align with phospho columns. |
| `total` | `pandas.DataFrame`, `str`, `Path`, or `None` | `None` | Total-protein matrix. | Used only by total-protein correction or protein-aware preparation. |
| `site_resolution_mode` | `str` | `"site_level_resolved"` | Selects site-level or peptide-evidence input. | Supports `"site_level_resolved"` and `"peptide_evidence"`. |
| `peptide_evidence` | `pandas.DataFrame`, `str`, `Path`, or `None` | Required for `peptide_evidence` | Peptide-level evidence table. | Must contain the metadata needed for peptide-to-site interpretation. |
| `peptide_evidence_sample_intensity_columns` | `tuple[str, ...]` or `None` | Required for `peptide_evidence` | Names the quantitative sample columns in `peptide_evidence`. | Values must be unique and present. |
| `peptide_site_mapping` | `pandas.DataFrame`, `str`, `Path`, or `None` | `None` | Optional explicit peptide-to-site mapping. | Must agree with the selected resolution policy. |
| `multi_site_policy` | `str` or `None` | Required when ambiguous peptide evidence is possible | Controls ambiguous multi-site peptides. | Use one of the supported policies listed below. |
| `allow_opaque_site_values` | `bool` | `False` | Allows site values outside strict serine, threonine, or tyrosine residue/position notation. | Enable only for a documented nonstandard use case. |
| `organism` | `Organism` or `None` | `None` | Species represented by the dataset. | Must agree with supplied identity metadata. |
| `preprocessing_config` | `DatasetPreprocessingConfig` | `DatasetPreprocessingConfig()` | Groups dataset preprocessing policies. | Incompatible stages fail before construction. |
| `input_intensity_scale` | `IntensityScaleKind`, `str`, or `None` | `None` | Declares an already established `"linear"` or `"log2"` scale when no transform establishes it. | Suspicious declarations fail by default. |
| `allow_suspicious_declared_input_intensity_scale` | `bool` | `False` | Allows a suspicious declared log2 scale with recorded warnings. | Use only when independent evidence supports the declaration. |
| `quantitative_meaning` | `QuantitativeMeaning`, `str`, or `None` | `None` | Declares what the supplied matrix represents. | Must agree with the established scale and cannot claim an operation PhosPy did not perform. |
| `corrected_preprocessing_output` | `CorrectedPreprocessingOutput` or `None` | `None` | Supplies externally corrected output at the supported boundary. | Cannot be combined with downstream matrix-changing preprocessing stages. |

Constructing the request stores these choices. Validation happens when
`AnalysisReadyDatasetBuilder().run(request)` is called.

## Preprocessing

`DatasetPreprocessingConfig` keeps related policies together. Defaults are
conservative: no transform, no normalisation, no correction, missing values
forbidden, and localisation required at 0.75.

| Configuration Group | Default | Purpose |
| --- | --- | --- |
| `intensity_transform` | `policy="identity"` | Keeps values unchanged or applies a log2 transform. |
| `normalisation` | `policy="none"` | Keeps distributions unchanged, median-centres samples, or applies quantile normalisation. |
| `missing_data` | `policy="forbid"` | Rejects missing values or applies a selected imputation method. |
| `group_coverage_filter` | `enabled=False` | Filters rows by finite-value coverage within sample groups. |
| `localisation` | `mode="require_threshold"` | Enforces localisation evidence at dataset build time. |
| `batch_correction` | `method="none"` | Applies fixed-effect residualisation or explicit native SPS/RUV-style correction. |
| `total_protein_correction` | `policy="none"` | Optionally subtracts matched log2 total-protein values. |
| `protein_aware_preparation` | `policy="disabled"` | Prepares aligned phosphosite/protein model inputs without changing intensities. |
| `site_matrix` | `policy="as_input"` | Controls duplicate-site handling and optional metadata-based construction. |
| `site_sequence_resolution` | no FASTA path | Optionally validates or fills sequences from a local FASTA file. |
| `comparisons` | `policy="none"` | Optionally builds descriptive pairwise comparison columns. |
| `ruv_readiness` | `enabled=False` | Records report-only readiness metadata; it does not modify the matrix. |

Common presets are available as
`DatasetPreprocessingConfig.strict()` and
`DatasetPreprocessingConfig.from_raw_phosphosite_table()`.

<details markdown="1">
<summary><strong>Common Preprocessing Parameters</strong></summary>

### Intensity Transform

`DatasetIntensityTransformConfig` accepts `policy="identity"` or
`policy="log2"`. For log2 transformation, `pseudocount` defaults to 1.0 and
must be finite and nonnegative.

When you keep the identity policy, declare the supplied scale explicitly:

```python
from phospy.api import IntensityScaleKind

request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    input_intensity_scale=IntensityScaleKind.LINEAR,
)
```

A suspicious declared log2 scale fails by default rather than being silently
accepted or transformed again. Any deliberate override is recorded in the
preprocessing report and provenance.

### Normalisation

`DatasetNormalisationConfig.policy` supports:

| Value | Effect |
| --- | --- |
| `"none"` | Keeps sample distributions unchanged. |
| `"median_center"` | Subtracts each sample's median. |
| `"quantile"` | Makes sample columns share one empirical distribution. |

Use quantile normalisation only when that distributional assumption is
scientifically appropriate.

### Missing Data

`DatasetMissingDataConfig.policy` supports:

| Value | Main Parameters | Notes |
| --- | --- | --- |
| `"forbid"` | None | Fails when missing values remain. |
| `"impute_row_median"` | `min_observed_values`, `input_scale` | Replaces missing cells with the row median after the row passes the observation threshold. |
| `"impute_minprob"` | `q`, `width`, `seed`, `max_missing_fraction_per_row`, `input_scale` | Stochastic lower-tail imputation; set a seed for reproducibility. |
| `"impute_knn"` | `k`, `distance`, `max_missing_fraction_per_row`, `input_scale`, `no_overlap_policy` | KNN imputation with explicit missingness limits. |

The dataset retains observation metadata so downstream workflows can distinguish
originally observed values from imputed replacements. Imputation can affect
scientific inference; inspect the workflow-specific policy before analysis.

### Group Coverage Filter

`DatasetGroupCoverageFilterConfig` uses `sample_metadata[group_column]` to keep
rows with sufficient finite observations in a requested number of groups.
Choose either `min_finite_observations_per_group` or
`min_finite_fraction_per_group`, then set `min_groups_passing_threshold`.
PhosPy does not infer groups from sample names.

### Localisation

`DatasetLocalisationConfig` accepts:

| Mode | Meaning |
| --- | --- |
| `"require_threshold"` | Requires every retained row to contain a valid confidence value at or above `min_confidence`. |
| `"allow_missing_with_waiver"` | Allows missing evidence only with a nonempty `waiver_reason`. |
| `"ignore"` | Deliberately ignores localisation evidence. |

The strict mode fails when the confidence column is absent, values are invalid
or missing, or any retained value is below `min_confidence`.

</details>

<details markdown="1">
<summary><strong>Site Construction, Sequences, and Comparisons</strong></summary>

### Site Matrix

`DatasetSiteMatrixConfig` controls how rows become a site matrix.

| Parameter | Default | Supported Values |
| --- | --- | --- |
| `policy` | `"as_input"` | `"as_input"`, `"build_from_metadata"` |
| `duplicate_site_policy` | `"error"` | `"error"`, `"max_mean_signal"`, `"first"`, `"aggregate_mean"`, `"aggregate_median"` |
| `missing_data_policy` | `"drop_any_missing"` | `"drop_any_missing"` |
| `minimum_observed_values` | `None` | Must remain `None` at the strict public boundary. |

The default rejects duplicate `site_key` rows. When you choose a non-error
policy, inspect `dataset.preprocessing_report.duplicate_site_resolution` and
`dataset.preprocessing_report.metadata_conflicts`.

### Site-Sequence Resolution

`DatasetSiteSequenceResolutionConfig` can validate, fill, or replace
`site_sequence` values from a local FASTA file.

| Parameter | Default | Supported Values or Meaning |
| --- | --- | --- |
| `fasta_path` | `None` | Local FASTA path; no online lookup occurs. |
| `mode` | `"validate_existing_and_fill_missing"` | Also supports `"fill_missing_only"`, `"validate_existing_only"`, and `"replace_existing"`. |
| `conflict_policy` | `None` | `"error"`, `"preserve_existing"`, or `"replace_existing"`. |
| `flank_size` | `7` | Residues requested on each side of the modified residue. |
| `accession_column` | `"protein_accession"` | Site-metadata accession column. |
| `site_column` | `"site"` | Site label column, such as `S939`. |

### Comparison Building

`DatasetComparisonBuildingConfig(policy="sample_metadata_pairs")` can build
pairwise dataset columns from `sample_metadata`. Use `sample_group_column` to
name the group field and `pairs` for explicit group pairs. These descriptive
comparisons do not replace the differential workflow's design and contrast
contracts.

### RUV Readiness

`ruv_readiness` diagnostics are report-only. `DatasetRuvReadinessConfig`
records whether named control-feature, replicate-group, and optional batch
columns are available. It does not select controls, estimate factors, or change
the matrix.

</details>

## Batch Correction

PhosPy keeps two correction methods separate because they answer different
technical questions.

### Limited Fixed-Effect Residualisation

`linear_residualize_batch`, a limited fixed-effect residualisation, removes
configured batch terms while including condition terms in the design. It is
not ComBat, not RUV, not limma `removeBatchEffect` parity, and not mixed-effects
modelling.

```python
from phospy.advanced import DatasetBatchCorrectionConfig

batch_correction = DatasetBatchCorrectionConfig(
    method="linear_residualize_batch",
    batch_column="batch",
    condition_column="condition",
    preserve_condition_effects=True,
)
```

This method requires aligned `sample_metadata`. A confounded batch and condition
design fails before correction because the condition effect cannot be
preserved. Inspect `dataset.preprocessing_report.batch_correction` for its
status, observed levels, confounding-check status, matrix shapes, warnings, and
limitations.

### Native SPS/RUV-Style Correction

Native PhosPy SPS/RUV-style preprocessing correction estimates unwanted factors
from eligible control-site residuals after protected-design handling. The
caller supplies the controls, protected condition columns, missingness policy,
and number of factors. PhosPy does not fetch or silently choose controls.

Batch terms are resolved for validation and diagnostics; they are not directly
residualized as fixed effects by the native correction. This implementation is
not PhosR-equivalent SPS/RUV-III parity. Replicate-aware RUV-III semantics are
not implemented.

The optional `replicate_column` is checked and recorded for diagnostics and
provenance. The `replicate_column` metadata is not used for numerical
unwanted-factor estimation and does not enable RUV-III behavior.

<details markdown="1">
<summary><strong>Native Correction Configuration and Safety Boundaries</strong></summary>

Control metadata is explicit. Caller-supplied controls must provide auditable
control-source metadata or field-level `metadata_missing_reason` rationale. The
audited fields include `organism`, `identifier_namespace`, source identity,
`source_version`, `license`, and `redistribution`. PhosPy rejects incompatible
organism or namespace metadata when present. Formal or external source names
require source version. Packaged controls, if added, require complete organism,
namespace, source, version, license, and redistribution metadata. PhosPy does
not infer metadata from `site_key` strings and does not fetch metadata online.

Native correction requires:

- a caller-supplied `ControlSiteSet` with auditable source metadata;
- aligned `sample_metadata` with `batch_column` and all protected
  `condition_columns`;
- an explicit `CorrectionMissingnessPolicy`;
- `n_unwanted_factors >= 1`, supported by the eligible controls and residual
  rank; and
- diagnostics plus `provenance_enabled=True`.

Multiple protected condition columns are represented as observed joint condition
strata, such as `condition=treated|timepoint=early`. This is not additive
protected-condition modelling, and the method does not fit additive
protected-condition terms. Correction runs after missing-data handling
and before downstream matrix consumers. Other stage placements are rejected.

#### Minimal Valid Native-Correction Example

Minimal valid native-correction example:

Replace the example keys with `site_key` values present in your dataset input.

```python
from phospy.advanced import (
    ControlSiteSet,
    ControlSiteSourceMetadata,
    CorrectionMissingnessPolicy,
    SpsRuvBatchCorrectionConfig,
)
from phospy.api import DatasetPreprocessingConfig

control_site_keys = (
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=MAPK14|residue=Y|position=182",
    "phospy:v1|organism=rat|protein_namespace=protein_id|"
    "protein_identifier=SRC|residue=Y|position=416",
)

control_source = ControlSiteSourceMetadata(
    organism="rat",
    identifier_namespace="site_key",
    source_name="manual-curated-controls",
    source_version="manual-v1",
    license="caller local use",
    redistribution="not redistributed",
)

native_correction = SpsRuvBatchCorrectionConfig(
    control_site_set=ControlSiteSet.from_site_keys(
        control_site_keys,
        source_metadata=control_source,
    ),
    batch_column="batch",
    condition_columns=("condition",),
    replicate_column="replicate",
    missingness_policy=CorrectionMissingnessPolicy(),
    n_unwanted_factors=1,
    diagnostics_enabled=True,
    provenance_enabled=True,
)

preprocessing = DatasetPreprocessingConfig(
    batch_correction=native_correction,
)
```

The public native SPS/RUV-style workflow requires a complete correction-stage
matrix and rejects actual missing values (NaNs) before executor invocation.
Upstream-imputed cells may remain identified with observation-mask provenance
through an `ObservationMask`; they are not treated as observed evidence.

```python
from phospy.advanced import (
    ObservationMask,
    OriginallyMissingCellTracking,
    TemporaryImputationMethod,
    TemporaryImputationPolicy,
)

mask = ObservationMask(
    feature_ids=("site_a", "site_b", "site_c"),
    sample_ids=("sample_1", "sample_2", "sample_3", "sample_4"),
    originally_missing_cells=(("site_b", "sample_2"),),
)

missingness_policy = CorrectionMissingnessPolicy(
    temporary_imputation=TemporaryImputationPolicy(
        allowed=True,
        method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
        method_parameters={"min_observed_values": 2},
    ),
    originally_missing_cells_tracked_by=(
        OriginallyMissingCellTracking.OBSERVATION_MASK
    ),
    observation_mask=mask,
)
```

`TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY` is a recorded correction
mechanic. It does not let actual NaNs pass through the public native workflow,
and it is not a way to pass actual NaNs through the public native workflow or
to reclassify imputed cells as observed evidence.

#### Rejected Unsafe Example

Rejected unsafe example:

This configuration requests one factor from only one control site. Validation
rejects it instead of selecting fallback controls.

```python
unsafe_correction = SpsRuvBatchCorrectionConfig(
    control_site_set=ControlSiteSet.from_site_keys(
        (control_site_keys[0],),
        source_metadata=control_source,
    ),
    batch_column="batch",
    condition_columns=("condition",),
    replicate_column="replicate",
    missingness_policy=CorrectionMissingnessPolicy(),
    n_unwanted_factors=1,
    provenance_enabled=True,
)
```

Successful correction records typed provenance, selected and rejected controls,
factor diagnostics, singular values, variance summaries, observation-mask
fingerprints, warnings, and input/output matrix fingerprints.

</details>

Externally supplied `CorrectedPreprocessingOutput` is accepted only at its safe
boundary. Provide it as the only matrix-changing preprocessing input; do not
combine it with total-protein correction, site-matrix construction,
normalisation, or comparison building.

## Total-Protein Options

`DatasetTotalProteinCorrectionConfig(policy="subtract_log_total")` calculates
log2 phosphosite minus matched log2 total protein. It requires a log2 transform,
an aligned `total` table, and explicit identity matching.

```python
from phospy.advanced import (
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
)

correction = DatasetTotalProteinCorrectionConfig(
    policy="subtract_log_total",
    identity=DatasetTotalProteinCorrectionIdentityConfig(
        mode="direct",
        phosphosite_key="protein_id",
        total_protein_key="__index__",
    ),
)
```

Identity configuration supports direct matching or an explicit mapping table,
strict matching, duplicate handling, and unmatched-row handling. Inspect the
processing state when you intentionally allow partial correction.

### Protein-Aware Preparation

`DatasetProteinAwarePreparationConfig` prepares aligned phosphosite and protein
model inputs and diagnostics. It does not change the phosphosite matrix, does
not subtract total protein, does not normalise intensities, and does not run
differential analysis. It also does not claim MSstatsPTM-style inference.

```python
from phospy.advanced import DatasetProteinAwarePreparationConfig

preprocessing = DatasetPreprocessingConfig(
    protein_aware_preparation=DatasetProteinAwarePreparationConfig(
        policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
    )
)
```

After building the dataset:

```python
preparation = dataset.protein_aware_preparation
report = dataset.preprocessing_report.protein_aware_preparation
site_eligibility = report.site_eligibility_dataframe()
```

The current differential workflow does not consume this preparation result.

## Peptide-Evidence Input

Use `site_resolution_mode="peptide_evidence"` when the input contains
peptide-level evidence rather than resolved site rows. Ambiguous multi-site
peptides require an explicit policy.

Supported analysis-ready policies are:

- `reject`: Fail when ambiguous peptide evidence is present.
- `exclude_from_sequence_scoring`: Remove ambiguous peptide evidence before the analysis-ready dataset is built.
- `split`: Allocate ambiguous peptide evidence to strict site-level rows.

The former `keep_joint` value is not supported by
`AnalysisReadyDatasetBuilder`. Joint tokens retain unresolved peptide ambiguity
and cannot provide the strict residue/position identity required by
`AnalysisReadyPhosphoDataset`.

```python
request = DatasetBuildRequest(
    site_resolution_mode="peptide_evidence",
    peptide_evidence=peptide_evidence,
    peptide_evidence_sample_intensity_columns=("sample_a", "sample_b"),
    multi_site_policy="split",  # replaces the removed "keep_joint" builder policy
    organism=Organism.HUMAN,
    input_intensity_scale="linear",
)
```

## Response

`AnalysisReadyDatasetBuilder.run(...)` returns an
`AnalysisReadyPhosphoDataset`. Data-frame properties return defensive snapshots,
so changing a returned table does not mutate the dataset.

| Attribute | Format | Meaning | Presence |
| --- | --- | --- | --- |
| `phospho` | `pandas.DataFrame` | Analysis-ready site-by-sample quantitative matrix indexed by `site_key`. | Always |
| `site_metadata` | `pandas.DataFrame` | Required identity, site, sequence, and scientific metadata aligned to `phospho`. | Always |
| `sample_metadata` | `pandas.DataFrame` or `None` | Aligned sample metadata supplied to the builder. | When supplied |
| `total` | `pandas.DataFrame` or `None` | Aligned total-protein input. | When supplied |
| `comparisons` | `pandas.DataFrame` or `None` | Optional builder-created comparison columns. | When configured |
| `intensity_scale_state` | `IntensityScaleState` | Established intensity scale and supporting evidence. | Always |
| `processing_state` | `DatasetProcessingState` | Typed preprocessing and quantitative-meaning state. | Always |
| `preprocessing_report` | `DatasetPreprocessingReport` or `None` | Row counts, operations, attrition, correction, sequence, and conflict diagnostics. | Built datasets normally include it |
| `imputation_observation_metadata` | Typed metadata or `None` | Originally observed versus imputed-cell evidence. | When imputation metadata exists |
| `imputation_feature_metadata` | `pandas.DataFrame` or `None` | Per-feature imputation summary. | When imputation metadata exists |
| `protein_aware_preparation` | `ProteinAwarePreparationResult` or `None` | Prepared protein-aware model inputs and report. | When configured |
| `organism` | `Organism` or `None` | Dataset species. | When established |
| `opaque_site_values_allowed` | `bool` | Whether the advanced opaque-site policy was enabled. | Always |
| `provenance` | `RunProvenance` or `None` | Inputs, parameters, stage evidence, fingerprints, and caveats. | Built datasets normally include it |
| `reference_context` | Reference context or `None` | Reference compatibility information derived from provenance. | When available |
| `trusted_construction_assertions` | `TrustedDatasetConstructionAssertions` or `None` | Evidence supplied through the advanced/trusted construction route. | Trusted-table construction only |

The most useful preprocessing report tables are:

| Attribute | Meaning |
| --- | --- |
| `row_counts` | Row counts at named preprocessing stages. |
| `operations` | Applied operations, parameters, shapes, and status. |
| `row_audit` | Row-level retention and removal reasons. |
| `duplicate_site_resolution` | How duplicate `site_key` rows were handled. |
| `metadata_conflicts` | Conflicting metadata discovered during row resolution. |
| `batch_correction` | Typed correction status and diagnostics. |
| `protein_aware_preparation` | Typed preparation report. |

## Common Problems

### The intensity scale is not established

Apply an explicit log2 transform or declare the supplied scale with
`input_intensity_scale`. Do not label raw-looking linear values as log2 merely
to pass validation.

### Localisation validation fails

Check the configured confidence column, missing values, and threshold. Use a
waiver only when the missing evidence and its scientific consequences are
understood and recorded.

### Sample metadata does not align

The `sample_metadata` index and `total` columns must match the phospho sample
columns. PhosPy does not reorder ambiguous or duplicate identifiers silently.

### Duplicate sites are rejected

Confirm whether the rows truly represent the same protein-scoped site. Choose a
non-error duplicate policy only when retaining or aggregating those rows is a
deliberate scientific decision.

### Missing values remain

The analysis-ready boundary is complete by design. Filter insufficient rows or
select an explicit imputation policy, then inspect the recorded observation
metadata before downstream inference.

## Related Documentation

- [Run your first analysis](../quickstart.md)
- [Differential analysis](differential-analysis.md)
- [Kinase analysis](kinase.md)
- [Scientific interpretation and limitations](../scientific-interpretation.md)

# Dataset Build Workflow

This page explains the dataset build API in detail. Start here when you have
phosphosite intensity data and want a strict `AnalysisReadyPhosphoDataset` for
differential, kinase, and signalome analysis.

## Purpose

The dataset builder validates table shape, interprets site metadata, applies the
preprocessing policy you request, and returns an `AnalysisReadyPhosphoDataset`.
The dataset that leaves the builder is intentionally strict and must be
missing-value-free. It also uses protein-scoped `site_key` values as
analysis-ready row identity. `display_id` is the human-readable `GENE;SITE;`
label and may repeat when distinct `site_key` values preserve the protein
context.

Direct `AnalysisReadyPhosphoDataset` construction is stricter than builder
ingestion: direct construction requires encoded `site_key` indexes and
auditable identity metadata (`site_key`, `display_id`, `organism`,
`protein_namespace`, `protein_identifier`, `gene_symbol`, `site`, and
`site_sequence`). It does not silently fall back to display-site identity. The
builder may accept legacy display-indexed input only when `site_metadata`
contains enough protein context to derive `site_key` without ambiguity.

`DatasetBuildRequest` is a lightweight command payload. Constructing it stores
the requested inputs and policies, but does not prove the dataset-build request
is valid. `AnalysisReadyDatasetBuilder.run(...)` owns request validation and
rejects invalid source types, unsupported site-resolution modes,
preprocessing/input incompatibilities, and scientific dataset states before
interpretation or dataset construction.

```python
dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LINEAR,
    )
)
```

## Imports

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBatchCorrectionConfig,
    DatasetBuildRequest,
    DatasetComparisonBuildingConfig,
    DatasetGroupCoverageFilterConfig,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetProteinAwarePreparationConfig,
    DatasetRuvReadinessConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
    IntensityScaleKind,
    Organism,
)
```

## Request Parameters

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `phospho` | `pandas.DataFrame`, `str`, or `pathlib.Path` | None | Yes | Site-by-sample intensity matrix. Rows are phosphosites and columns are samples. |
| `site_metadata` | `pandas.DataFrame`, `str`, or `pathlib.Path` | None | Yes | Row metadata aligned to `phospho.index` at ingestion. The analysis-ready boundary requires `site_key`, `display_id`, `organism`, `protein_namespace`, `protein_identifier`, `gene_symbol`, `site`, and `site_sequence`. Builder input may omit `site_key` only when it includes enough protein context, preferably `protein_identifier` plus `protein_namespace`, to derive it. `protein_id` is additionally required for signalome. For site-level scientific workflows, include a localisation-confidence column (default: `localisation_confidence`) and configure explicit localisation policy. |
| `sample_metadata` | `pandas.DataFrame`, `str`, or `pathlib.Path`, or `None` | `None` | No | Descriptive/alignment metadata aligned to phospho columns with unique column names. Required when comparison building uses `sample_metadata_pairs`. It does not automatically define differential-analysis conditions, replicates, batches, or blocks. |
| `total` | `pandas.DataFrame`, `str`, or `pathlib.Path`, or `None` | `None` | No | Total-protein matrix used only when total-protein correction is enabled. Columns must align to phospho sample columns. |
| `organism` | `Organism` or `None` | `None` | No | Species identity for the dataset. Use `Organism.RAT` for the bundled beginner lane. |
| `preprocessing_config` | `DatasetPreprocessingConfig` | `DatasetPreprocessingConfig()` | No | Grouped preprocessing policy for transforms, normalisation, missing data, group-aware coverage filter declaration, optional `linear_residualize_batch` batch residualisation, total-protein correction, protein-aware preparation, site construction, site-sequence resolution, comparisons, and RUV readiness reporting. |
| `input_intensity_scale` | `IntensityScaleKind`, `str`, or `None` | `None` | No | Required when your preprocessing path keeps `intensity_transform.policy="identity"` and you still need a trusted intensity scale (`"linear"` or `"log2"`). |
| `quantitative_meaning` | `QuantitativeMeaning`, `str`, or `None` | `None` | No | Optional explicit scientific meaning for phospho values (for example `phosphosite_abundance` or `phosphosite_log_abundance`). |

Supported file suffixes are `.csv`, `.tsv`, `.txt` as tab-separated text, and
`.parquet`. CSV, TSV, and TXT inputs are read with the first column as the row
index.

Parsing is table-role aware:

- `phospho` and `total` are parsed as numeric matrices with string-preserved row
  and column identifiers. Non-numeric cells fail fast with row/column context.
- `site_metadata` and `sample_metadata` are parsed as string-preserving metadata
  tables (`"NA"` stays the literal string `"NA"`; leading-zero IDs are
  preserved).
- Missing-value interpretation is explicit per table role; dataset loading does
  not rely on pandas default NA inference for metadata.

## Sample Metadata Semantics

`sample_metadata` is passive dataset metadata aligned to `phospho.columns`.

- It is useful for descriptive labels and metadata-driven preprocessing features.
- Its column names must be unique.
- It does not perform scientific design validation.
- It does not automatically define differential-analysis conditions,
  replicates, batches, or blocks.
- Differential analysis consumes explicit workflow design objects
  (`ExperimentalDesign` and `Contrast`) in the differential workflow API.

```python
request = DatasetBuildRequest(
    phospho="./input/phospho.csv",
    site_metadata="./input/site_metadata.tsv",
    sample_metadata="./input/sample_metadata.csv",
    organism=Organism.RAT,
    input_intensity_scale=IntensityScaleKind.LINEAR,
)
```

## Minimum Input Shape

`phospho` should be numeric. Builder input may use human-readable display labels
such as `TSC2;S939;` as the index when enough protein context is available in
`site_metadata` to derive `site_key`.

```python
phospho = pd.DataFrame(
    {
        "sample_a": [1.00, 0.70],
        "sample_b": [1.10, 0.80],
        "sample_c": [0.95, 0.75],
    },
    index=["TSC2;S939;", "GSK3B;S9;"],
)
```

`site_metadata` must align to `phospho.index`.
This builder example intentionally omits `site_key`; direct
`AnalysisReadyPhosphoDataset` construction cannot omit it.

```python
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "ATMSGRPRTTSFAESSSPVQQPSAFGQAAAL",
        ],
        "display_id": ["TSC2;S939;", "GSK3B;S9;"],
        "organism": ["rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id"],
        "protein_identifier": ["TSC2", "GSK3B"],
        "protein_id": ["TSC2", "GSK3B"],
        "localisation_confidence": [0.95, 0.92],
    },
    index=phospho.index.copy(),
)
```

The built dataset will be indexed by `site_key`:

```python
assert dataset.phospho.index.name == "site_key"
assert dataset.site_metadata.index.name == "site_key"
assert dataset.site_metadata["display_id"].tolist() == ["TSC2;S939;", "GSK3B;S9;"]
assert dataset.site_metadata["site_key"].tolist() == dataset.phospho.index.tolist()
assert {"organism", "protein_namespace", "protein_identifier"}.issubset(
    dataset.site_metadata.columns
)
```

### Duplicate Display Labels

`display_id` is not a row key. It is display metadata. Two rows can have the
same human-readable label when protein-scoped identity differs:

```python
phospho = pd.DataFrame(
    {
        "sample_a": [1.0, 2.0],
        "sample_b": [1.1, 2.1],
    },
    index=["MAPK14;Y182;", "MAPK14;Y182;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["MAPK14", "MAPK14"],
        "site": ["Y182", "Y182"],
        "site_sequence": [
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
        ],
        "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
        "organism": ["rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id"],
        "protein_identifier": ["MAPK14_A", "MAPK14_B"],
        "protein_id": ["MAPK14_A", "MAPK14_B"],
        "localisation_confidence": [0.95, 0.95],
    },
    index=phospho.index.copy(),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        input_intensity_scale=IntensityScaleKind.LINEAR,
    )
)

identity = dataset.site_metadata.loc[
    :,
    ["site_key", "display_id", "protein_identifier"],
]
assert identity["display_id"].tolist() == ["MAPK14;Y182;", "MAPK14;Y182;"]
assert identity["site_key"].is_unique
assert identity["protein_identifier"].tolist() == ["MAPK14_A", "MAPK14_B"]
assert dataset.phospho.shape[0] == 2
```

Rows are preserved because their `site_key` values differ. They are not
collapsed, overwritten, aggregated, or deduplicated merely because
`display_id` repeats. A `GENE;SITE;` display label is therefore not sufficient
biological row identity for site-level analysis.

Supported site-metadata aliases are deliberately narrow:

- `gene_name` may stand in for `gene_symbol`.
- `centralized_sequence` may stand in for `site_sequence`.

The builder may derive `gene_symbol` and `site` from index values formatted like
`MAPK14;Y182;`. It does not derive `protein_id` or any protein-scoped identity
from the gene-symbol prefix.

Protein context is used to derive `site_key` when it is available and safe.
`display_id` remains metadata and may repeat after `site_key` becomes the row
identity. If direct construction is used instead of the builder, the caller must
provide matching `site_key` indexes and all required identity metadata up front.

## Preprocessing Configuration

`DatasetPreprocessingConfig` groups the public preprocessing controls.

| Parameter | Type | Default | How to Use It |
| --- | --- | --- | --- |
| `intensity_transform` | `DatasetIntensityTransformConfig` | `policy="identity"`, `pseudocount=1.0` | Controls numeric intensity transformation before downstream preprocessing. |
| `normalisation` | `DatasetNormalisationConfig` | `policy="none"` | Controls sample-wise normalisation. |
| `missing_data` | `DatasetMissingDataConfig` | `policy="forbid"` | Controls missing-value rejection or imputation. |
| `group_coverage_filter` | `DatasetGroupCoverageFilterConfig` | `enabled=False` | Filters phosphosite rows by finite-value coverage within sample groups before missing-data handling. |
| `total_protein_correction` | `DatasetTotalProteinCorrectionConfig` | `policy="none"` | Controls phosphosite-to-total correction. |
| `protein_aware_preparation` | `DatasetProteinAwarePreparationConfig` | `policy="disabled"` | Prepares aligned phosphosite/protein model-input contracts and diagnostics only. |
| `batch_correction` | `DatasetBatchCorrectionConfig` | `method="none"` | Controls optional `linear_residualize_batch` fixed-effect residualisation. |
| `site_matrix` | `DatasetSiteMatrixConfig` | `policy="as_input"` | Controls construction and duplicate-site handling. |
| `site_sequence_resolution` | `DatasetSiteSequenceResolutionConfig` | `mode="validate_existing_and_fill_missing"` | Controls optional local FASTA-backed site-sequence resolution. |
| `comparisons` | `DatasetComparisonBuildingConfig` | `policy="none"` | Controls optional pairwise comparison construction from sample metadata. |
| `localisation` | `DatasetLocalisationConfig` | `mode="require_threshold"`, `min_confidence=0.75`, `confidence_column="localisation_confidence"` | Controls site-level localisation-confidence eligibility at dataset-build time. |
| `ruv_readiness` | `DatasetRuvReadinessConfig` | `enabled=False` | Adds readiness reporting for possible future RUV/SPS/RUV-III preprocessing; it does not select SPS controls or apply correction. |

Use presets for common lanes:

```python
strict = DatasetPreprocessingConfig.strict()
raw_table = DatasetPreprocessingConfig.from_raw_phosphosite_table()
```

If you keep the identity transform (`policy="identity"`), declare
`input_intensity_scale` on `DatasetBuildRequest` so the dataset can establish a
trusted scale state:

```python
request = DatasetBuildRequest(
    phospho=phospho,
    site_metadata=site_metadata,
    input_intensity_scale=IntensityScaleKind.LINEAR,
)
```

Use explicit groups when you need a specific preprocessing policy:

```python
config = DatasetPreprocessingConfig(
    intensity_transform=DatasetIntensityTransformConfig(
        policy="log2",
        pseudocount=1.0,
    ),
    normalisation=DatasetNormalisationConfig(policy="median_center"),
    missing_data=DatasetMissingDataConfig(
        policy="impute_row_median",
        min_observed_values=2,
    ),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=config,
    )
)
```

## Batch-Correction Parameters

`DatasetBatchCorrectionConfig` defaults to `method="none"`. The only supported
executable batch-related preprocessing method is
`method="linear_residualize_batch"`, a
fixed-effect residualisation of batch terms that preserves condition effects by
including condition terms in the design. It is not ComBat, not RUV, not limma
`removeBatchEffect` parity, and not mixed-effects modelling.

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `method` | `str` | `"none"` | `"none"`, `"linear_residualize_batch"` | Selects whether batch residualisation runs. |
| `batch_column` | `str` | `"batch"` | Non-empty string | Column in `sample_metadata` identifying batch labels. |
| `condition_column` | `str` | `"condition"` | Non-empty string | Column in `sample_metadata` identifying condition labels to preserve during residualisation. |
| `preserve_condition_effects` | `True` | `True` | `True` only | Condition preservation is required for `linear_residualize_batch`. |

```python
batch_correction = DatasetBatchCorrectionConfig(
    method="linear_residualize_batch",
    batch_column="batch",
    condition_column="condition",
    preserve_condition_effects=True,
)
```

When this method is requested, `sample_metadata` is required and must align to
the phospho sample columns. Confounded batch/condition designs are rejected
before correction because condition effects cannot be preserved in those
designs. Inspect `dataset.preprocessing_report.batch_correction` after build
for status, observed levels, confounding-check status, matrix shapes, warnings,
and limitations.

## Intensity Transform Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"identity"` | `"identity"`, `"log2"` | `"identity"` keeps values as provided. `"log2"` applies `log2(value + pseudocount)`. |
| `pseudocount` | `float` or `int` | `1.0` | Any finite value `>= 0` | Used only by `policy="log2"`. Choose it to keep all transformed inputs valid. |

Example:

```python
transform = DatasetIntensityTransformConfig(
    policy="log2",
    pseudocount=1.0,
)
```

## Normalisation Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"none"` | `"none"`, `"median_center"`, `"quantile"` | `"none"` keeps sample distributions unchanged. `"median_center"` subtracts sample-wise medians. `"quantile"` forces sample columns to share one empirical distribution. |

Quantile normalisation is dense and sort-heavy. Use it only when matched sample
distributions are scientifically appropriate.

```python
normalisation = DatasetNormalisationConfig(policy="median_center")
```

Normalisation provenance is explicit in both `dataset.provenance.preprocessing_stages`
and `dataset.preprocessing_report.operations`. Stage diagnostics include:

- method and resolved parameters
- input/output matrix shapes
- per-sample summaries before and after normalisation
- row/column drop indicators and counts

## Missing-Data Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"forbid"` | `"forbid"`, `"impute_row_median"`, `"impute_minprob"`, `"impute_knn"` | Selects rejection or imputation behaviour. |
| `min_observed_values` | `int` or `None` | `None` | Integer `>= 1` when `policy="impute_row_median"`; otherwise `None` | Drops rows with too few observed samples before row-median imputation. |
| `q` | `float` or `None` | `None` | `0 < q < 0.5` when `policy="impute_minprob"`; otherwise `None` | Lower-tail quantile for MinProb-style imputation. |
| `width` | `float` or `None` | `None` | `0 < width <= 1.0` when `policy="impute_minprob"`; otherwise `None` | Controls the width of the imputation distribution. |
| `seed` | `int` or `None` | `None` | Integer `>= 0` when `policy="impute_minprob"`; otherwise `None` | Makes random MinProb-style imputation reproducible. |
| `k` | `int` or `None` | `None` | Integer `>= 1` when `policy="impute_knn"`; otherwise `None` | Number of neighbours for KNN imputation. |
| `distance` | `str` or `None` | `None` | `"nan_euclidean"` when `policy="impute_knn"`; otherwise `None` | Distance metric for KNN imputation. |
| `max_missing_fraction_per_row` | `float` or `None` | `None` | `0 < value <= 1` for `"impute_minprob"` and `"impute_knn"`; otherwise `None` | Drops rows with too much missingness before advanced imputation. |

Row-median example:

```python
missing_data = DatasetMissingDataConfig(
    policy="impute_row_median",
    min_observed_values=2,
)
```

MinProb-style example:

```python
missing_data = DatasetMissingDataConfig(
    policy="impute_minprob",
    q=0.01,
    width=0.3,
    seed=12345,
    max_missing_fraction_per_row=0.5,
)
```

KNN example:

```python
missing_data = DatasetMissingDataConfig(
    policy="impute_knn",
    k=5,
    distance="nan_euclidean",
    max_missing_fraction_per_row=0.5,
)
```

## Group Coverage Filter Parameters

`DatasetGroupCoverageFilterConfig` describes a group-aware coverage filter
rule. It can express rules such as "keep sites quantified in at least two
replicates in at least one condition." Use it before analysis when rows with
too little replicate or condition coverage should be removed before imputation,
normalisation, and analysis-ready dataset creation.

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `enabled` | `bool` | `False` | `True`, `False` | Enables filtering. When `False`, existing dataset-building behavior is unchanged. |
| `group_column` | `str` or `None` | `None` | Non-empty string when `enabled=True` | Sample-metadata column that defines groups such as conditions. |
| `min_finite_observations_per_group` | `int` or `None` | `None` | Integer `>= 1`, mutually exclusive with fraction threshold | Minimum finite sample values needed within a group. |
| `min_finite_fraction_per_group` | `float` or `None` | `None` | `0 < value <= 1`, mutually exclusive with count threshold | Minimum finite-value fraction needed within a group. |
| `min_groups_passing_threshold` | `int` | `1` | Integer `>= 1` | Number of groups that must pass the selected threshold. |

The filter uses `sample_metadata[group_column]` to resolve groups. It does not
guess experimental groups from sample names. For each phosphosite row and each
group, it counts finite numeric values in that group's samples. A count
threshold such as `2` means at least two finite values in the group. A fraction
threshold such as `0.67` means at least 67% of samples in that group have finite
values. The row is kept when at least `min_groups_passing_threshold` groups pass.

```python
sample_metadata = pd.DataFrame(
    {"condition": ["control", "control", "control", "treated", "treated", "treated"]},
    index=["c1", "c2", "c3", "t1", "t2", "t3"],
)

coverage_filter = DatasetGroupCoverageFilterConfig(
    enabled=True,
    group_column="condition",
    min_finite_observations_per_group=2,
    min_groups_passing_threshold=1,
)

preprocessing = DatasetPreprocessingConfig(
    group_coverage_filter=coverage_filter,
    missing_data=DatasetMissingDataConfig(
        policy="impute_row_median",
        min_observed_values=1,
    ),
)
```

Filtering runs before analysis-ready dataset creation and before missing-data
imputation. Inspect `dataset.preprocessing_report.row_counts`,
`dataset.preprocessing_report.operations`, and
`dataset.preprocessing_report.row_audit` to see how many rows were retained or
removed and why. If all phosphosite rows are removed, dataset building fails
with a clear input error instead of creating an empty analysis-ready dataset.

## Localisation-Confidence Parameters

Use `DatasetLocalisationConfig` to make localisation policy explicit for
site-level analysis datasets.

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `mode` | `str` | `"require_threshold"` | `"require_threshold"`, `"allow_missing_with_waiver"`, `"ignore"` | `"require_threshold"` is the scientifically strict lane for site-level workflows. |
| `min_confidence` | `float` | `0.75` | `0.0 <= value <= 1.0` | Minimum accepted localisation confidence for `mode="require_threshold"`. |
| `confidence_column` | `str` | `"localisation_confidence"` | Non-empty string | Site-metadata column containing localisation confidence values. |
| `waiver_reason` | `str` or `None` | `None` | Non-empty string when `mode="allow_missing_with_waiver"` | Required reason when localisation strictness is waived. |

Failure behaviour for `mode="require_threshold"`:

- fails dataset build when `site_metadata[confidence_column]` is missing
- fails dataset build when values are missing/invalid
- fails dataset build when any value is below `min_confidence`

Why this matters: low-confidence phosphosite localisation can mis-assign
site-level effects and lead to misleading downstream kinase/signalome
interpretation.

```python
localisation = DatasetLocalisationConfig(
    mode="require_threshold",
    confidence_column="localisation_confidence",
    min_confidence=0.75,
)
```

## Total-Protein Correction Parameters

`DatasetTotalProteinCorrectionConfig` controls whether total-protein correction
runs.

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"none"` | `"none"`, `"subtract_log_total"` | `"subtract_log_total"` computes `log2_phospho - log2_total`. |
| `identity` | `DatasetTotalProteinCorrectionIdentityConfig` | `DatasetTotalProteinCorrectionIdentityConfig()` | Identity config object | Explains how phosphosite rows match total-protein rows. |

`policy="subtract_log_total"` requires `intensity_transform.policy="log2"`, a
`total` table aligned to phospho sample columns, and an explicit identity
mapping.

```python
config = DatasetPreprocessingConfig(
    intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
    total_protein_correction=DatasetTotalProteinCorrectionConfig(
        policy="subtract_log_total",
        identity=DatasetTotalProteinCorrectionIdentityConfig(
            mode="direct",
            phosphosite_key="protein_id",
            total_protein_key="__index__",
        ),
    ),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        total=total,
        preprocessing_config=config,
    )
)
```

### Total-Protein Identity Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `mode` | `str` | `"direct"` | `"direct"`, `"mapping_table"` | `"direct"` matches a phosphosite metadata key directly to a total-protein key. `"mapping_table"` uses an explicit mapping table. |
| `phosphosite_key` | `str` | `"gene_symbol"` | Non-empty string | Column in `site_metadata` used as the phosphosite-side identity. |
| `total_protein_key` | `str` | `"__index__"` | Non-empty string | Total-protein identity key. Use `"__index__"` to match against `total.index`. |
| `mapping_table` | `pandas.DataFrame` or `None` | `None` | Required for `mode="mapping_table"`; otherwise `None` | Two-column table connecting phosphosite identities to total-protein identities. |
| `mapping_phosphosite_key` | `str` or `None` | `None` | Required for `mode="mapping_table"`; otherwise `None` | Column in `mapping_table` that matches `phosphosite_key`. |
| `mapping_total_protein_key` | `str` or `None` | `None` | Required for `mode="mapping_table"`; otherwise `None` | Column in `mapping_table` that matches total-protein identities. |
| `matching_policy` | `str` | `"strict"` | `"strict"`, `"gene_symbol_normalised"` | `"strict"` compares trimmed identity keys exactly. `"gene_symbol_normalised"` uppercases gene-symbol keys and should be chosen only when that is scientifically suitable. |
| `duplicate_policy` | `str` | `"error"` | `"error"` | Fails on duplicate identity matches. |
| `unmatched_policy` | `str` | `"error"` | `"error"`, `"allow_uncorrected"` | `"error"` requires complete matching. `"allow_uncorrected"` keeps unmatched phosphosite rows uncorrected and marks the dataset as mixed. |

Mapping-table example:

```python
mapping_table = pd.DataFrame(
    {
        "phosphosite_protein_id": ["TSC2", "GSK3B"],
        "total_row_id": ["TSC2_total", "GSK3B_total"],
    }
)

identity = DatasetTotalProteinCorrectionIdentityConfig(
    mode="mapping_table",
    phosphosite_key="protein_id",
    total_protein_key="__index__",
    mapping_table=mapping_table,
    mapping_phosphosite_key="phosphosite_protein_id",
    mapping_total_protein_key="total_row_id",
    unmatched_policy="error",
)
```

Quantitative meaning is explicit after preprocessing:

- fully corrected log2 datasets: `phospho_total_log_ratio`
- log2 datasets without total-protein correction: `phosphosite_log_abundance`
- mixed corrected and uncorrected datasets: `mixed_phospho_total_log_ratio_and_phosphosite_log_abundance`

For mixed datasets, row-level correction status is available in
`dataset.processing_state.total_protein_correction.diagnostics`.

## Protein-Aware Preparation Parameters

`DatasetProteinAwarePreparationConfig` is a preparation-only policy for aligned
phosphosite/protein model inputs and diagnostics. It does not change the
phosphosite matrix, does not subtract total protein, does not normalise
intensities, and does not run differential analysis.

Use it only when you provide `total` input data and want an auditable
`ProteinAwarePreparationResult` on the built dataset:

```python
from phospy.api import (
    DatasetPreprocessingConfig,
    DatasetProteinAwarePreparationConfig,
)

preprocessing = DatasetPreprocessingConfig(
    protein_aware_preparation=DatasetProteinAwarePreparationConfig(
        policy="prepare_model_inputs",
        protein_mapping_policy="require_unambiguous",
    )
)
```

The public preparation stage maps explicit protein identifiers from
`site_metadata` (`protein_accession`, `protein_id`, or `protein_group_id`) to
`total.index`. Gene-symbol matching is not the public default. Missing or
ambiguous mappings are reported according to the configured mapping policy.

After build, inspect:

```python
preparation = dataset.protein_aware_preparation
report = dataset.preprocessing_report.protein_aware_preparation
site_eligibility = report.site_eligibility_dataframe()
```

PhosPy does not claim MSstatsPTM-style inference or equivalence for this
preparation stage. Current `DifferentialAnalysisWorkflow` execution does not
consume `ProteinAwarePreparationResult`.

## Site-Matrix Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"as_input"` | `"as_input"`, `"build_from_metadata"` | `"as_input"` preserves interpreted site-matrix-ready rows after `site_key` derivation. `"build_from_metadata"` constructs site identity from metadata before final dataset construction. |
| `duplicate_site_policy` | `str` | `"error"` | `"error"`, `"max_mean_signal"`, `"first"`, `"aggregate_mean"`, `"aggregate_median"` | Controls duplicate rows that resolve to the same `site_key` during site-matrix construction. Duplicate `display_id` values with distinct `site_key` values are valid. |
| `missing_data_policy` | `str` | `"drop_any_missing"` | `"drop_any_missing"` | Keeps only complete rows for strict dataset construction. |
| `minimum_observed_values` | `None` | `None` | `None` | Public strict construction requires this to stay unset. |

When two or more rows resolve to the same `site_key`, the default policy raises
instead of collapsing evidence. Duplicate `site_key` rows are a scientific
ambiguity because choosing one row or aggregating rows changes the
analysis-ready phosphosite evidence model. Non-error policies are deliberate
scientific choices: `max_mean_signal` and `first` retain one source row, while
`aggregate_mean` and `aggregate_median` combine duplicate source rows. When you
use a non-error policy, inspect
`dataset.preprocessing_report.duplicate_site_resolution` and
`dataset.preprocessing_report.metadata_conflicts`.

Rows with repeated `display_id` values can pass when their derived `site_key`
values differ. `display_id` is display metadata, not row identity.

```python
site_matrix = DatasetSiteMatrixConfig(
    policy="build_from_metadata",
    duplicate_site_policy="aggregate_mean",  # intentional row-collapse policy
)
```

## Site-Sequence Resolution Parameters

`DatasetSiteSequenceResolutionConfig` optionally uses a local FASTA file to
validate, fill, or replace site sequences.

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `fasta_path` | `str` or `None` | `None` | Local filesystem path or `None` | Provides the FASTA source. When omitted, no local FASTA file is read. |
| `mode` | `str` | `"validate_existing_and_fill_missing"` | `"validate_existing_and_fill_missing"`, `"fill_missing_only"`, `"validate_existing_only"`, `"replace_existing"` | Controls whether existing sequences are checked, missing sequences are filled, or existing values are replaced. |
| `conflict_policy` | `str` or `None` | `None` | `"error"`, `"preserve_existing"`, `"replace_existing"`, or `None` | Controls conflicts between supplied and FASTA-derived sequences. |
| `flank_size` | `int` | `7` | Integer `>= 1` when `fasta_path` is provided | Number of residues requested on each side of the modified residue. |
| `accession_column` | `str` | `"protein_accession"` | Non-empty string | Metadata column containing protein accessions used for FASTA lookup. |
| `site_column` | `str` | `"site"` | Non-empty string | Metadata column containing site labels such as `S939`. |

Example:

```python
site_sequence_resolution = DatasetSiteSequenceResolutionConfig(
    fasta_path="./references/rat.fasta",
    mode="validate_existing_and_fill_missing",
    conflict_policy="error",
    flank_size=7,
    accession_column="protein_accession",
    site_column="site",
)
```

## Comparison-Building Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"none"` | `"none"`, `"sample_metadata_pairs"` | Selects whether dataset comparisons are built. |
| `sample_group_column` | `str` | `"comparison_group"` | Non-empty string | Column in `sample_metadata` containing group labels. |
| `pairs` | `tuple[tuple[str, str], ...]` or `None` | `None` | Non-empty pair tuples when provided | Explicit comparison pairs. If omitted, pairs are inferred from observed groups. |

When `policy="sample_metadata_pairs"`, `sample_metadata` is required.
This comparison-building feature does not replace differential workflow design
contracts.

```python
sample_metadata = pd.DataFrame(
    {"comparison_group": ["control", "treated", "treated"]},
    index=["sample_a", "sample_b", "sample_c"],
)

comparisons = DatasetComparisonBuildingConfig(
    policy="sample_metadata_pairs",
    sample_group_column="comparison_group",
    pairs=(("control", "treated"),),
)
```

## RUV-Readiness Parameters

`DatasetRuvReadinessConfig` is report-only. It helps audit whether metadata
needed for possible future RUV/SPS/RUV-III preprocessing is present; it does not
select SPS controls or correct the matrix.
It does not make sample metadata scientific design input for differential
analysis.

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `enabled` | `bool` | `False` | `True`, `False` | Enables readiness reporting. |
| `control_feature_column` | `str` | `"is_control_feature"` | Non-empty string | Site-metadata column identifying control features. |
| `replicate_group_column` | `str` | `"replicate_group"` | Non-empty string | Sample-metadata column identifying replicate groups. |
| `batch_column` | `str` or `None` | `"batch"` | Non-empty string or `None` | Sample-metadata column identifying batches. Use `None` when batch information is not available. |

```python
ruv_readiness = DatasetRuvReadinessConfig(
    enabled=True,
    control_feature_column="is_control_feature",
    replicate_group_column="replicate_group",
    batch_column="batch",
)
```

## Full Dataset-Build Example

```python
import pandas as pd

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetLocalisationConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    Organism,
)

phospho = pd.DataFrame(
    {
        "sample_a": [100.0, 70.0],
        "sample_b": [110.0, 80.0],
        "sample_c": [95.0, 75.0],
    },
    index=["TSC2;S939;", "GSK3B;S9;"],
)
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "ATMSGRPRTTSFAESSSPVQQPSAFGQAAAL",
        ],
        "display_id": ["TSC2;S939;", "GSK3B;S9;"],
        "organism": ["rat", "rat"],
        "protein_namespace": ["protein_id", "protein_id"],
        "protein_identifier": ["TSC2", "GSK3B"],
        "protein_id": ["TSC2", "GSK3B"],
        "localisation_confidence": [0.95, 0.92],
    },
    index=phospho.index.copy(),
)

preprocessing = DatasetPreprocessingConfig(
    intensity_transform=DatasetIntensityTransformConfig(
        policy="log2",
        pseudocount=1.0,
    ),
    normalisation=DatasetNormalisationConfig(policy="median_center"),
    missing_data=DatasetMissingDataConfig(policy="forbid"),
    localisation=DatasetLocalisationConfig(
        mode="require_threshold",
        confidence_column="localisation_confidence",
        min_confidence=0.75,
    ),
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=preprocessing,
    )
)

print(dataset.phospho.shape)
print(dataset.intensity_scale_state.label)
```

# Dataset Build Workflow

This page explains the dataset build API in detail. Start here when you have
phosphosite intensity data and want a strict `AnalysisReadyPhosphoDataset` for
kinase and signalome analysis.

## Purpose

The dataset builder validates table shape, interprets site metadata, applies the
preprocessing policy you request, and returns an `AnalysisReadyPhosphoDataset`.
The dataset that leaves the builder is intentionally strict and must be
missing-value-free.

```python
dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
    )
)
```

## Imports

```python
from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetComparisonBuildingConfig,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetRuvReadinessConfig,
    DatasetSiteMatrixConfig,
    DatasetSiteSequenceResolutionConfig,
    DatasetTotalProteinCorrectionConfig,
    DatasetTotalProteinCorrectionIdentityConfig,
    Organism,
)
```

## Request Parameters

| Parameter | Type | Default | Required | How to Use It |
| --- | --- | --- | --- | --- |
| `phospho` | `pandas.DataFrame`, `str`, or `pathlib.Path` | None | Yes | Site-by-sample intensity matrix. Rows are phosphosites and columns are samples. |
| `site_metadata` | `pandas.DataFrame`, `str`, or `pathlib.Path` | None | Yes | Row metadata aligned to `phospho.index`. `gene_symbol`, `site`, and `site_sequence` are required at the analysis-ready boundary; include `protein_id` for signalome. |
| `sample_metadata` | `pandas.DataFrame`, `str`, or `pathlib.Path`, or `None` | `None` | No | Sample metadata aligned to phospho columns. Required when comparison building uses `sample_metadata_pairs`. |
| `total` | `pandas.DataFrame`, `str`, or `pathlib.Path`, or `None` | `None` | No | Total-protein matrix used only when total-protein correction is enabled. Columns must align to phospho sample columns. |
| `organism` | `Organism` or `None` | `None` | No | Species identity for the dataset. Use `Organism.RAT` for the bundled beginner lane. |
| `preprocessing_config` | `DatasetPreprocessingConfig` | `DatasetPreprocessingConfig()` | No | Grouped preprocessing policy for transforms, normalisation, missing data, total-protein correction, site construction, site-sequence resolution, comparisons, and RUV readiness reporting. |

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

```python
request = DatasetBuildRequest(
    phospho="./input/phospho.csv",
    site_metadata="./input/site_metadata.tsv",
    sample_metadata="./input/sample_metadata.csv",
    organism=Organism.RAT,
)
```

## Minimum Input Shape

`phospho` should be numeric. Its index should use stable PhosPy site IDs such as
`TSC2;S939;`.

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

```python
site_metadata = pd.DataFrame(
    {
        "gene_symbol": ["TSC2", "GSK3B"],
        "site": ["S939", "S9"],
        "site_sequence": [
            "FDDTPEKDSFRARSTSLNERPKSLRIARAPK",
            "ATMSGRPRTTSFAESCKPVQQPSAFGQAAAL",
        ],
        "protein_id": ["TSC2", "GSK3B"],
    },
    index=phospho.index.copy(),
)
```

Supported site-metadata aliases are deliberately narrow:

- `gene_name` may stand in for `gene_symbol`.
- `centralized_sequence` may stand in for `site_sequence`.

The builder may derive `gene_symbol` and `site` from index values formatted like
`MAPK14;Y182;`. It does not derive `protein_id`.

## Preprocessing Configuration

`DatasetPreprocessingConfig` groups the public preprocessing controls.

| Parameter | Type | Default | How to Use It |
| --- | --- | --- | --- |
| `intensity_transform` | `DatasetIntensityTransformConfig` | `policy="identity"`, `pseudocount=1.0` | Controls numeric intensity transformation before downstream preprocessing. |
| `normalisation` | `DatasetNormalisationConfig` | `policy="none"` | Controls sample-wise normalisation. |
| `missing_data` | `DatasetMissingDataConfig` | `policy="forbid"` | Controls missing-value rejection or imputation. |
| `total_protein_correction` | `DatasetTotalProteinCorrectionConfig` | `policy="none"` | Controls phosphosite-to-total correction. |
| `site_matrix` | `DatasetSiteMatrixConfig` | `policy="as_input"` | Controls construction and duplicate-site handling. |
| `site_sequence_resolution` | `DatasetSiteSequenceResolutionConfig` | `mode="validate_existing_and_fill_missing"` | Controls optional local FASTA-backed site-sequence resolution. |
| `comparisons` | `DatasetComparisonBuildingConfig` | `policy="none"` | Controls optional pairwise comparison construction from sample metadata. |
| `ruv_readiness` | `DatasetRuvReadinessConfig` | `enabled=False` | Adds readiness reporting for future RUV-compatible workflows; it does not apply correction. |

Use presets for common lanes:

```python
strict = DatasetPreprocessingConfig.strict()
raw_table = DatasetPreprocessingConfig.from_raw_phosphosite_table()
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

## Site-Matrix Parameters

| Parameter | Type | Default | Allowed Values | How to Use It |
| --- | --- | --- | --- | --- |
| `policy` | `str` | `"as_input"` | `"as_input"`, `"build_from_metadata"` | `"as_input"` preserves interpreted site-matrix-ready rows. `"build_from_metadata"` constructs site IDs from metadata. |
| `duplicate_site_policy` | `str` | `"max_mean_signal"` | `"max_mean_signal"`, `"first"`, `"aggregate_mean"`, `"aggregate_median"`, `"error"` | Controls duplicate-site collapse when `policy="build_from_metadata"`. |
| `missing_data_policy` | `str` | `"drop_any_missing"` | `"drop_any_missing"` | Keeps only complete rows for strict dataset construction. |
| `minimum_observed_values` | `None` | `None` | `None` | Public strict construction requires this to stay unset. |

Duplicate-site handling is a scientific choice. `"error"` is cautious,
`"first"` is input-order dependent, `"max_mean_signal"` favours stronger rows,
and aggregate policies can blur peptide context.

```python
site_matrix = DatasetSiteMatrixConfig(
    policy="build_from_metadata",
    duplicate_site_policy="error",
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
needed for future RUV-compatible correction is present; it does not correct the
matrix.

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
            "ATMSGRPRTTSFAESCKPVQQPSAFGQAAAL",
        ],
        "protein_id": ["TSC2", "GSK3B"],
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

# Validation Guide

PhosPy validates early and loudly so scientific assumptions are not hidden. Most
errors are fixable once you know which boundary rejected the input.

## Dataset Input Rules

`phospho` must be a non-empty numeric pandas `DataFrame` or supported file path.
Rows are phosphosites and columns are samples. Builder input may use display
labels such as `MAPK14;Y182;` as the index when `site_metadata` provides enough
protein context to derive `site_key`. Direct `AnalysisReadyPhosphoDataset`
construction must already use encoded `site_key` row indexes; display-indexed
direct construction is invalid. Missing values are rejected by default.

`site_metadata` must be a non-empty table aligned to `phospho.index`. It must
include non-empty `gene_symbol`, `site`, and `site_sequence` columns at the
analysis-ready boundary, plus auditable protein context (`organism`,
`protein_namespace`, and `protein_identifier`). `site_sequence` may be omitted
at ingestion only when preprocessing can derive it before final dataset
construction. `protein_id` is not part of dataset row identity; it is optional
at the analysis-ready dataset boundary, may be absent, and may contain missing
or blank values until a workflow explicitly requires it. Signalome is the
workflow that requires complete `protein_id` values as algorithm-specific
protein grouping metadata.

`sample_metadata`, when provided, must align to the phospho sample columns and
must have unique column names.

`total`, when provided, must be numeric, missing-value-free, and aligned to the
phospho sample columns.

## Site Metadata Conventions

Accepted column aliases are narrow:

| Accepted alias | Normalised column |
| --- | --- |
| `gene_name` | `gene_symbol` |
| `centralized_sequence` | `site_sequence` |

If `gene_symbol` or `site` is missing, the builder can derive them from an
input index like `TSC2;S939;`. It does not derive `protein_id` from the
gene-symbol prefix, and it does not treat the display prefix as protein context
or signalome grouping metadata. Builder ingestion may accept legacy
display-indexed input only when enough protein context exists to derive
`site_key`. Direct analysis-ready construction must provide `site_key`; it does
not silently fall back to `GENE;SITE;` display labels.

## Analysis-Ready Dataset Boundary

A built `AnalysisReadyPhosphoDataset` must have:

- numeric, non-empty, missing-value-free `phospho`
- unique sample columns
- unique `site_key` values
- `phospho.index.name == "site_key"`
- `site_metadata.index` exactly matching `phospho.index`
- `site_metadata["site_key"]` exactly matching `site_metadata.index`
- required `display_id`; repeated `display_id` values are valid when
  `site_key` values differ
- required non-empty `organism`, `protein_namespace`, `protein_identifier`,
  `gene_symbol`, `site`, and `site_sequence`
- `site_metadata["site_key"]` matching the metadata-derived
  (`organism`, `protein_namespace`, `protein_identifier`, `site`) key
- optional `protein_id`, which is signalome grouping metadata when signalome is
  run, may be incomplete at the dataset boundary, and is not a replacement for
  `protein_identifier` or `site_key`
- `sample_metadata.index` exactly matching `phospho.columns` when provided
- `sample_metadata.columns` unique when provided
- `total.columns` exactly matching `phospho.columns` when provided
- an `Organism` enum value or `None`
- explicit intensity-scale and processing-state metadata

## Preprocessing Rules

Defaults are intentionally strict: no transform, no normalisation, no imputation,
no total-protein correction, and no comparison construction.

Common cross-field checks:

- missing-data handling runs before normalisation in preprocessing stage order.
- when `intensity_transform.policy="identity"`, declare
  `input_intensity_scale` on `DatasetBuildRequest` (`"linear"` or `"log2"`), or
  use an explicit scale-changing transform (for example `policy="log2"`).
- `missing_data.policy="impute_row_median"` is deterministic.
- row-median imputation is not left-censored imputation.
- imputed row-median values are replacements and must not be treated as evidence that the original values were observed.
- `missing_data.policy="impute_minprob"` requires `intensity_transform.policy="log2"`.
- `impute_minprob` requires explicit `q`, `width`, `seed`, and `max_missing_fraction_per_row`.
- `impute_minprob` is left-censored random imputation with deterministic seeded draws and row-drop reporting above the configured missing-fraction threshold.
- `missing_data.policy="impute_knn"` requires explicit `k`, `distance="nan_euclidean"`, and `max_missing_fraction_per_row`.
- `impute_knn` requires `min_observed_values=None` and does not support alternative distance metrics in the public contract.
- `impute_knn` drops rows above `max_missing_fraction_per_row`, reports dropped rows as not imputable, and must produce a complete matrix.
- `subtract_log_total` requires `total` input data.
- `subtract_log_total` requires `intensity_transform.policy="log2"`.
- When `subtract_log_total` runs with `unmatched_policy="allow_uncorrected"` and
  unmatched phosphosite rows are retained, dataset quantitative meaning is set to
  `mixed_phospho_total_log_ratio_and_phosphosite_log_abundance`.
- `sample_metadata_pairs` requires `sample_metadata`.
- site-matrix construction may drop incomplete rows because the public output
  dataset must be complete.
- duplicate rows resolving to the same `site_key` are a scientific ambiguity and
  fail by default.
- non-error duplicate-site policies are deliberate row-retention or
  row-collapse choices; when they run, inspect
  `dataset.preprocessing_report.duplicate_site_resolution` and
  `dataset.preprocessing_report.metadata_conflicts`.
- `ruv_readiness.enabled=True` records readiness signals for possible future
  RUV/SPS/RUV-III preprocessing, including complete-matrix status, missingness
  mask provenance, control-feature availability, replicate groups, and optional
  batch metadata. It does not select SPS controls or apply correction.
- RUV readiness is informational: dataset construction is not rejected when
  readiness is false. Native SPS/RUV-style correction is separate and requires
  explicit `SpsRuvBatchCorrectionConfig`.

## Reference Validation

`ReferenceBundle` requires:

- `organism` as an `Organism` enum value
- `kinase_substrate_map` with non-empty `kinase` and `substrate_site`
- `site_sequences` indexed by display site ID with non-empty `site_sequence`
- no duplicate `(kinase, substrate_site)` pairs

Reference `substrate_site` and `site_sequences.index` values are display IDs at
the reference boundary. Kinase workflow interpretation maps those display IDs
through dataset `display_id` metadata onto internal `site_key` rows. Reference
validation does not convert display IDs into analysis-ready row identity.

`ReferencePreset.AUTO` uses `dataset.organism`. In the current release, bundled runtime
references are rat-only.

## Workflow Validation

### Differential Workflow

`DifferentialAnalysisRequest.dataset` must be an
`AnalysisReadyPhosphoDataset` whose phospho matrix is already indexed by
`site_key`. Differential result tables are strict public contracts: each
per-contrast table is indexed by encoded protein-scoped `site_key` values and
includes non-empty `site_key`, `display_id`, `organism`, `protein_namespace`,
`protein_identifier`, `gene_symbol`, and `site` columns, with `site_key`
exactly matching the index. Workflow-created results preserve that required
protein context from dataset `site_metadata` and optional protein metadata such
as `protein_id` when present.

Display-indexed, `GENE;SITE;`-keyed, arbitrary non-encoded, or stat-only
differential result tables are invalid in public `DifferentialAnalysisResult`
construction. Result validation does not infer identity from `gene_symbol`,
`site`, or display labels, and it does not repair missing identity columns.

### Kinase Workflow

`KinaseWorkflowRequest.dataset` must be an `AnalysisReadyPhosphoDataset`.
References must be compatible with the dataset organism when organism information
is present. Kinase scoring and prediction operate on `site_key`; display IDs are
used only through the explicit reference-mapping layer described above.

`KinaseScoringConfig.min_substrates` must be at least `2`. The activity stage can
be disabled with `activity_config=None`, which is useful for tiny examples.
Mixed corrected/uncorrected quantitative meaning is rejected by default; set
`scoring_config.allow_mixed_total_protein_quantitative_meaning=True` to opt in.

### Signalome Workflow

`SignalomeWorkflowRequest.kinase_result` must be a `KinaseWorkflowResult`.
Signalome also requires explicit, non-empty `protein_id` values for every
interpreted site as signalome-specific protein grouping metadata. Gene-symbol
prefixes in display labels are not treated as protein grouping metadata or
protein identity.
Signalome aligns dataset, prediction, and score tables by `site_key` and does
not reinterpret display IDs as row identity.
Mixed corrected/uncorrected quantitative meaning is rejected by default; set
`config.validation.allow_mixed_total_protein_quantitative_meaning=True` to opt in.

Signalome scale guards protect expensive clustering work:

- `performance.max_exact_tree_sites` limits exact tree construction.
- `performance.max_full_candidate_scoring_sites` limits full candidate
  correlation scoring.
- `clustering.candidate_scoring_policy="sampled"` can reduce candidate-scoring
  cost but still needs exact tree construction.

After a successful run, `result.provenance.workflow_parameters["scale_guard"]`
shows exact tree-generation details and candidate-scoring details separately.

## Quick Fix Table

| Error shape | What to check first |
| --- | --- |
| unsupported file format | Use `.csv`, `.tsv`, `.txt`, or `.parquet`; install parquet support for `.parquet`. |
| missing `gene_symbol` or `site` | Add those columns or, for builder input only, use index labels formatted as `GENE;SITE;` with sufficient protein context. |
| missing protein-scoped identity metadata | Add non-empty `organism`, `protein_namespace`, `protein_identifier`, and `site`, or use a builder-compatible protein-context source that derives them before final construction. |
| display-indexed direct construction | Construct through the builder with enough protein context, or provide encoded `site_key` indexes and matching `site_metadata.site_key` directly. |
| signalome protein grouping metadata error | Add non-empty `protein_id` grouping metadata for every interpreted site; do not use `gene_symbol` or `display_id` as a fallback. |
| reference resolution error | Use rat with `AUTO`, or pass an explicit `ReferenceBundle`. |
| total-protein correction error | Provide `total`, set `intensity_transform.policy="log2"`, and configure identity mapping. |
| mixed quantitative meaning rejected | Use `unmatched_policy="error"` or complete total-protein mapping; if mixed inputs are intentional, set the workflow mixed-state opt-in flag. |
| activity error on a tiny example | Disable activity or provide enough supported substrates. |
| signalome scale error | Reduce sites, use `clustering.candidate_scoring_policy="sampled"` where appropriate, or raise `performance` guards deliberately. |

## Error Families

Common public exception families are:

- `PhosPyInputError`: file, table, or request input problem
- `PhosPyValidationError`: validated object does not satisfy its contract
- `PhosPyReferenceError`: reference resolution or compatibility problem
- `PhosPyWorkflowError`: workflow boundary or execution problem

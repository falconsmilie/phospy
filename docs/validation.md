# Validation Guide

PhosPy validates early and loudly so scientific assumptions are not hidden. Most
errors are fixable once you know which boundary rejected the input.

## Dataset Input Rules

`phospho` must be a non-empty numeric pandas `DataFrame` or supported file path.
Rows are phosphosites and columns are samples. Builder input may use display
labels such as `MAPK14;Y182;` as the index when `site_metadata` provides enough
protein context to derive `site_key`. The direct
`AnalysisReadyPhosphoDataset` constructor raises immediately. Advanced
trusted reconstruction of already prepared tables must use
`AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with typed evidence or
explicit waivers for identity, intensity scale, quantitative meaning,
localisation, sequence, and reference context, plus non-waivable aligned-table
structure evidence. Supplied trusted provenance must match the actual table
fingerprints.
Display-indexed direct construction is invalid. Missing values are rejected by
default.

`site_metadata` must be a non-empty table aligned to `phospho.index`. It must
include non-empty `gene_symbol`, `site`, and `site_sequence` columns at the
analysis-ready boundary, plus auditable protein context (`organism`,
`protein_namespace`, and `protein_identifier`). `site_sequence` may be omitted
at ingestion only when preprocessing can derive it before final dataset
construction. `protein_group_id` is not part of dataset row identity; it is
optional at the analysis-ready dataset boundary, may be absent, and may contain
missing or blank values until a workflow explicitly requires it. Signalome is
the workflow that requires complete `protein_group_id` values as
algorithm-specific protein grouping metadata. Legacy `protein_id` is accepted
only as a Signalome migration alias.

At the analysis-ready dataset boundary, `site_sequence` means required sequence
evidence for the row: a non-empty, plausible amino-acid context string aligned
to the phosphosite identity. It does not mean the value is already suitable for
motif scoring or any specific sequence-aware workflow.

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
input index like `TSC2;S939;`. It does not derive `protein_group_id` from the
gene-symbol prefix, and it does not treat the display prefix as protein context
or Signalome grouping metadata. Builder ingestion may accept legacy
display-indexed input only when enough protein context exists to derive
`site_key`. The trusted analysis-ready factory must receive `site_key`; it does
not silently fall back to `GENE;SITE;` display labels and cannot prove the
biological correctness of user-asserted provenance.

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
- optional `protein_group_id`, which is Signalome grouping metadata when
  Signalome is run, may be incomplete at the dataset boundary, and is not a
  replacement for `protein_identifier` or `site_key`. Legacy `protein_id` is a
  migration alias for this Signalome field.
- `sample_metadata.index` exactly matching `phospho.columns` when provided
- `sample_metadata.columns` unique when provided
- `total.columns` exactly matching `phospho.columns` when provided
- an `Organism` enum value or `None`
- explicit intensity-scale and processing-state metadata

Sequence readiness has two levels:

- Analysis-ready sequence evidence: `AnalysisReadyPhosphoDataset` requires
  `site_sequence` for every row and validates it as a plausible amino-acid
  context string.
- Workflow-specific sequence-context readiness: a sequence-aware workflow
  validates the selected dataset/reference sequence after request-specific
  resolution and conflict policy. That stricter contract may require centered
  phosphosite context, exact window length, center index, center residue,
  alphabet, terminal-padding policy, lowercase policy,
  modified-residue-symbol policy, and known sequence source.

The base dataset check is deliberately plausibility-level. Dataset construction
does not know which workflow, scoring mode, reference bundle, conflict policy, or
motif resource will later be selected, so it cannot prove biological correctness
or motif readiness for every `site_sequence` value.

## Preprocessing Rules

Defaults are intentionally strict: no transform, no normalisation, no imputation,
no total-protein correction, and no comparison construction.

Common cross-field checks:

- missing-data handling runs before normalisation in preprocessing stage order.
  Its placement relative to intensity transformation is method/scale policy:
  linear imputations run before a configured log2 transform, and log2
  imputations run after the transform or on input declared as already log2.
- when `intensity_transform.policy="identity"`, declare
  `input_intensity_scale` on `DatasetBuildRequest` (`"linear"` or `"log2"`), or
  use an explicit scale-changing transform (for example `policy="log2"`).
- high-confidence suspicious declared `log2` input scales fail by default when
  matrix values strongly resemble raw linear intensities. If the declaration is
  scientifically justified despite diagnostics, set
  `allow_suspicious_declared_input_intensity_scale=True`; the override flag,
  effective policy, input declaration source, and diagnostic warnings are
  recorded in provenance. Other declared-scale diagnostics are recorded as
  provenance warnings.
- `missing_data.policy="impute_row_median"` is deterministic and requires
  explicit `missing_data.input_scale` (`"linear"` or `"log2"`) during
  preprocessing plan interpretation.
- row-median imputation is not left-censored imputation.
- imputed row-median values are replacements and must not be treated as evidence that the original values were observed.
- `missing_data.policy="impute_minprob"` has method-required
  `missing_data.input_scale="log2"`. It requires either a configured log2
  transform before imputation or a dataset input scale declared as already
  log2.
- `impute_minprob` requires explicit `q`, `width`, `seed`, and `max_missing_fraction_per_row`.
- `impute_minprob` is left-censored random imputation with deterministic seeded draws and row-drop reporting above the configured missing-fraction threshold.
- `missing_data.policy="impute_knn"` requires explicit
  `missing_data.input_scale` (`"linear"` or `"log2"`), `k`,
  `distance="nan_euclidean"`, and `max_missing_fraction_per_row`.
- `impute_knn` requires `min_observed_values=None` and does not support alternative distance metrics in the public contract.
- `impute_knn` is a deterministic PhosPy-owned implementation, not a
  scikit-learn delegation. It drops rows above
  `max_missing_fraction_per_row`, reports dropped rows as not imputable, and
  must produce a complete matrix.
- for each retained missing cell, `impute_knn` considers only retained donor
  rows that have an observed value in that cell's column. Distances are
  nan-euclidean over target/donor shared observed columns and scaled by
  `n_columns / shared_observed_column_count`.
- `impute_knn` resolves equal-distance donors by `(str(row_id),
  original_position)`, averages selected donor values without distance
  weighting, and falls back to the retained-column mean when no donor has any
  shared observed value with the target row.
- `impute_knn` is chunked but guarded: retained matrices above the documented
  row, sample, or distance-work budgets fail with a `PhosPyInputError` that
  reports the shape and suggests reducing retained missing rows or choosing a
  simpler missing-data policy.
- Missing-data diagnostics and processing state preserve the observation mask,
  imputation input scale, and imputation operation order. Missing-data and
  total-protein correction diagnostics use the shared immutable JSON policy
  internally. Bundle payloads still serialize as schema v1 `dict`/`list` JSON,
  and each serialization returns fresh detached containers.
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
- `ruv_readiness.enabled=True` records report-only readiness signals for
  possible future RUV-family work, including complete-matrix status,
  missingness mask provenance, control-feature availability, replicate groups,
  and optional batch metadata. It does not select SPS controls or apply
  correction, and it does not make RUV-III correction executable.
- RUV readiness is informational: dataset construction is not rejected when
  readiness is false. Native SPS/RUV-style correction is separate and requires
  explicit `SpsRuvBatchCorrectionConfig`.
- Native SPS/RUV-style correction runs only through `batch_correction` with an
  explicit `SpsRuvBatchCorrectionConfig`. The native PhosPy SPS/RUV-style
  preprocessing correction estimates unwanted factors from eligible
  control-site residuals after protected-design handling. Batch terms are
  resolved for validation and diagnostics, including batch-associated-variance
  summaries; they are not directly residualized as fixed effects by the native
  correction. Validation rejects missing or misaligned `sample_metadata`,
  missing `batch_column`, missing protected `condition_columns`, confounded or
  rank-inadequate batch/condition designs, too few eligible caller-supplied
  control sites, incompatible control-site metadata, caller controls missing
  audit metadata without explicit `metadata_missing_reason`, conflicting
  control-source declarations such as disagreeing `source_type`/`source`,
  incomplete
  packaged-control metadata, duplicate or ambiguous accepted-control metadata,
  missing values without a supported `CorrectionMissingnessPolicy`, unsupported
  temporary-imputation policies, missing observation masks when temporary
  imputation is allowed, and attempts to run correction without provenance.
  Optional `replicate_column` metadata is validated and recorded for native-lane
  provenance and diagnostics only; it is not used for numerical unwanted-factor
  estimation and does not make RUV-III or replicate-aware RUV-III correction
  executable.
  Validation does not infer organism or identifier namespace from `site_key`
  strings and does not fetch metadata online.

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

`KinaseWorkflowRequest.scoring_config` must be explicit: use
`KinaseScoringConfig.exploratory()`, `KinaseScoringConfig.production(...)`, or
direct custom construction with `reliability_profile="custom"`.
`KinaseScoringConfig.min_substrates` must be at least `2`. The activity stage is
disabled by default with `activity_config=None`, which is useful for tiny
examples.
Mixed corrected/uncorrected quantitative meaning is rejected by default; set
`scoring_config.allow_mixed_total_protein_quantitative_meaning=True` to opt in.

Every current kinase scoring mode requires workflow-specific centered
phosphosite sequence context before execution. Kinase Library-style motif modes
add a fixed centered-window contract from the local
`KinaseLibraryResource.sequence_window`: the selected sequence must be present,
non-empty, exactly the configured upstream-plus-site-plus-downstream length,
centered at the configured index, and centered on `S`, `T`, or `Y`. Unsupported
characters, lowercase residues, modified-residue symbols, and terminal padding
are rejected unless the selected workflow contract explicitly allows them.
Opaque display labels are not sequence evidence, unknown sequence source is
invalid for Kinase Library motif scoring, and incompatible dataset/reference
sequences require an explicit conflict policy such as `prefer_dataset` or
`prefer_reference`.

### Signalome Workflow

`SignalomeWorkflowRequest.kinase_result` must be a `KinaseWorkflowResult`.
Signalome also requires explicit, non-empty `protein_group_id` values for every
interpreted site as Signalome-specific protein grouping metadata. Legacy
`protein_id` is accepted only as a migration alias, and conflicting
`protein_group_id`/`protein_id` values fail validation. Signalome uses this
field to group retained phosphosites into protein-level module and protein-site
context summaries. `protein_group_id` is not core protein identity; core protein
identity remains the dataset-level `organism`,
`protein_namespace`, and `protein_identifier` metadata. Gene-symbol prefixes in
display labels are not treated as protein grouping metadata or protein
identity.
Signalome aligns dataset, prediction, and score tables by `site_key` and does
not reinterpret display IDs as row identity.
Signalome validates that sequence-aware upstream site identity still provides a
centered phosphosite sequence context, including an `S`, `T`, or `Y` center that
matches the site residue, but it does not apply the fixed Kinase Library
motif-window length contract.
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
| display-indexed direct construction | Use the builder with enough protein context, or, for advanced/trusted construction, use `AnalysisReadyPhosphoDataset.from_trusted_tables(...)` with encoded `site_key` indexes and matching `site_metadata.site_key`. |
| signalome protein grouping metadata error | Add non-empty `protein_group_id` grouping metadata for every interpreted site. Legacy `protein_id` is accepted only as a migration alias. Keep core protein identity in `protein_namespace` and `protein_identifier`; do not use `gene_symbol` or `display_id` as a fallback. |
| workflow-specific sequence context error | Check selected `site_sequence`, sequence source, required window length, center index, center residue, alphabet, padding policy, and dataset/reference conflict policy for the workflow/scoring mode. |
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

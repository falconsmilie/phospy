# Validation Guide

This guide documents validation in the supported rewrite contract.

For API shapes, see [`api.md`](api.md).

## Boundary Model

Validation is split across three boundaries:

- Builder request boundary:
  validates supported source types (`DataFrame` or file path) and builder conventions.
- Final dataset boundary:
  validates strict `AnalysisReadyPhosphoDataset` invariants.
- Workflow boundaries:
  validate request DTOs, config ranges, reference compatibility, and runtime/science seams.

## Builder Flexibility vs Dataset Strictness

`DatasetBuildRequest` accepts:

- `pandas.DataFrame`
- file path (`str`, `Path`, `PathLike`) to `.csv`, `.tsv`/`.txt`, `.parquet`

Builder preprocessing config is nested under
`DatasetBuildRequest.preprocessing_config` and validated as one policy object
rather than root-level scalar fields.

This flexibility exists only at the builder input boundary.

After source loading, both routes pass through the same normalization path:

- site ID canonicalization
- site-metadata alias resolution
- optional derivation of `gene_symbol` and `site` from strict index format
- fail-fast rejection of ambiguous/unsupported legacy aliases

Before final dataset construction, supported preprocessing policy is applied to
`phospho` according to `preprocessing_config`:

- `"forbid"` (default): no missing-value preprocessing
- `"impute_row_median"`: drop rows below
  `missing_data.min_observed_values`, then row-median imputation for remaining
  missing phospho values
- `"ratio_to_total"` under `preprocessing_config.total_protein_correction.policy`:
  subtract matched total-protein abundance from phosphosite abundance before
  transformation-state establishment
- `"build_from_metadata"` under `preprocessing_config.site_matrix.policy`:
  construct site-matrix-ready rows from metadata after upstream missing-data and
  total/protein-correction stages

`AnalysisReadyPhosphoDataset` itself is strict and DataFrame-only.
Workflows consume only this dataset type.

## AnalysisReady Dataset Validation

`AnalysisReadyPhosphoDataset` constructor composes:

- `AnalysisReadyDatasetValidator`
- `TransformationStateValidator`

Enforced dataset invariants:

- `phospho`: non-empty numeric DataFrame, unique index/columns,
  canonical site IDs.
- `site_metadata`: non-empty DataFrame, exact index alignment with `phospho.index`,
  required columns `gene_symbol`, `site` with non-empty strings.
- `site_metadata.site_sequence` is optional; if present it must contain non-empty strings.
- `sample_metadata` (if present): index aligns to `phospho.columns`.
- `total` (if present): non-empty numeric DataFrame, unique index,
  columns align to `phospho.columns`.
- `organism` (if present): `Organism` enum.
- `transformation_state`: typed `TransformationState`, established through a
  supported PhosPy path, coherent with presence/absence of `total`, and coherent
  transformation kind between matrices when both exist.

Boundary constructors validate; they do not silently repair invalid data.

Transformation-state establishment is enforced at the final dataset boundary.
A coherent but externally declared state object is rejected unless it was
established through a supported PhosPy path (builder/resolver transformer lane
or supported bundle reconstruction lane).

In the current public builder lane, this establishment path is intentionally
narrow: builder execution establishes pass-through `linear` state only.
No additional transformation mode is publicly selectable.

## Builder Preprocessing Policy Rules

`DatasetPreprocessingConfig` validation enforces:

- grouped config type checks:
  `missing_data`, `total_protein_correction`, `site_matrix`, `comparisons`
- `missing_data.policy` must be one of `forbid`, `impute_row_median`
- `missing_data.min_observed_values` must be `None` when policy is `forbid`
- `missing_data.min_observed_values` must be an integer `>= 1` when policy is
  `impute_row_median`
- supported lane policy restrictions:
  - `total_protein_correction.policy` must be one of `none`, `ratio_to_total`
  - `site_matrix.policy` must be one of `as_input`, `build_from_metadata`
  - `comparisons.policy` must stay `none`
- at execution, `missing_data.min_observed_values` must not exceed phospho
  sample count
- at request validation, `total_protein_correction.policy='ratio_to_total'`
  requires `request.total` to be present
- at execution, `total_protein_correction.policy='ratio_to_total'` requires:
  - `total.columns` exactly matching `phospho.columns`
  - numeric phospho/total columns
  - unique normalized `total.index` identifiers
  - complete matching between `site_metadata.gene_symbol` and `total.index`
- at execution, `site_matrix.policy='build_from_metadata'` requires:
  - `site_metadata` columns `gene_symbol`, `site`, `site_sequence`
  - non-empty `gene_symbol`/`site` values for all rows
  - at least one retained row after dropping rows with missing `site_sequence`
    and rows with incomplete phospho values
  - duplicate constructed site IDs are collapsed by keeping the row with the
    strongest mean phospho signal

## Builder Convention Rules

Supported site-metadata aliases:

- `gene_symbol`: `gene_symbol`, `gene_name`
- `site`: `site`
- `site_sequence`: `site_sequence`, `centralized_sequence`
- `protein_id`: `protein_id`

Unsupported legacy aliases are rejected:

- `gene`
- `residue`
- `phosphosite`
- `site_position`
- `sequence`
- `protein`

If `gene_symbol` and/or `site` are missing, builder derivation is allowed only from
index values exactly matching `"<gene_symbol>;<site>;"`.
This derivation does not produce `protein_id`.

If both `site_metadata.index` and `site_metadata.site_id` are present, they must match
after canonicalization.

## Reference Validation

Reference resolution/compatibility enforces:

- `ReferencePreset.AUTO` requires `dataset.organism`.
- Explicit preset must match dataset organism when both are set.
- Explicit `ReferenceBundle.organism` must match dataset organism when both are set.
- Bundled runtime references are currently rat-only.
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` are valid enum values but are
  not bundled runtime lanes in this release.

`ReferenceBundle` structure/content validation enforces non-empty, canonical,
internally consistent reference tables.

## Workflow Validation

`KinaseWorkflowValidator` enforces:

- request type is `KinaseWorkflowRequest`
- `dataset` is `AnalysisReadyPhosphoDataset`
- `references` is `ReferencePreset | ReferenceBundle`
- config types and numeric floors/ranges:
  `scoring_config.min_substrates >= 2`,
  `scoring_config.include_diagnostic_scoring_tables` is bool,
  `prediction_config.top_k >= 1`,
  `prediction_config.ensemble_size >= 1`,
  activity config bounds when activity is enabled/configured

`SignalomeWorkflowValidator` enforces:

- request type is `SignalomeWorkflowRequest`
- `kinase_result` is `KinaseWorkflowResult`
- signalome config bounds in `[0.0, 1.0]`
- upstream downstream-score/prediction matrices are usable numeric matrices for
  signalome execution (`combined_scores` preferred when available, otherwise
  `profile_scores`)
- `kinase_result.dataset.site_metadata` is present and index-aligned to
  `kinase_result.dataset.phospho.index`
- `kinase_result.dataset.site_metadata.protein_id` is required and must contain
  non-empty string values
- missing values in the upstream downstream score matrix are allowed; this is a
  normal outcome for correlation-based kinase scoring in low-information rows
- infinite values in upstream score/prediction matrices remain hard failures

`SignalomeWorkflowInterpreter` preconditions the downstream score lane before
execution:

- rows with no finite kinase score support (all-missing rows) are excluded from
  score-driven network correlation inputs
- partially missing rows are retained and consumed with pairwise-complete
  correlation handling
- prediction rows remain available for module assignment logic
- drop policy is `allow_and_report`: non-zero dropped-row counts are surfaced in
  `signalome_result.score_preconditioning_diagnostics` instead of causing a
  boundary failure

Interpreters/executors enforce seam-level scientific/runtime boundary checks and raise
`WorkflowBoundaryError` with seam names, concrete counts, and `next_action` hints.

For the supported kinase scoring lane, motif sequence inputs come from
`references.site_sequences` (resolved reference bundle), not from
`dataset.site_metadata.site_sequence`.

## Nested Result Access (Validation-Relevant Contract)

Stable access paths are nested by stage:

- `result.scoring_result.profile_scores`
- `result.scoring_result.combined_scores`
- `result.scoring_result.motif_scores` (optional diagnostics)
- `result.scoring_result.weights` (optional diagnostics)
- `result.prediction_result.pred_mat`
- `result.activity_result.weighted_activity` (when activity is enabled)
- `signalome_result.kinase_result.prediction_result.pred_mat`

Optional outputs must be checked before dereference:

- `result.activity_result`
- `result.prediction_result.substrate_list`
- `signalome_result.kinase_network.nodes`
- `signalome_result.expanded_signalome`

## Validation Ownership Summary

| Invariant | Owner |
| --- | --- |
| Builder input source type checks | `DatasetBuildRequestValidator` + `DatasetInputSourceValidator` |
| Builder preprocessing config policy | `DatasetPreprocessingConfigValidator` |
| Builder convention normalization/derivation | `DatasetBuildRequestInterpreter` collaborators |
| Analysis-ready dataset structure/content | `AnalysisReadyDatasetValidator` |
| Transformation-state coherence and establishment | `TransformationStateValidator` |
| Reference compatibility (dataset vs preset/bundle) | `ReferenceCompatibilityValidator` / `ReferenceResolver` |
| Reference bundle structure/content | `ReferenceBundleValidator` |
| Kinase workflow request/config validity | `KinaseWorkflowValidator` + `WorkflowConfigValidator` |
| Signalome workflow request/config validity | `SignalomeWorkflowValidator` + `WorkflowConfigValidator` |
| Kinase/signalome runtime seam diagnostics | workflow interpreters/executors (`WorkflowBoundaryError`) |

## Quick Troubleshooting

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Builder rejects input format | Field is neither DataFrame nor supported file path | Pass DataFrame or path to `.csv`/`.tsv`/`.txt`/`.parquet` |
| Dataset constructor fails on site metadata | Required strict boundary columns/values are missing | Provide `gene_symbol` and `site` with non-blank strings |
| `ReferencePreset.AUTO` fails | Dataset organism is missing | Set `organism` in `DatasetBuildRequest` |
| Bundled human/mouse preset fails | Bundled references are rat-only in this release | Provide explicit non-rat `ReferenceBundle` |
| Kinase boundary seam fails | Overlap/support constraints were not met | Read seam details and adjust dataset/references/config |
| Signalome validation fails on `site_metadata.protein_id` | Explicit protein identity required for supported signalome execution | Provide non-empty `site_metadata.protein_id` for all interpreted sites; gene-symbol site-ID prefixes are not a substitute |

# Validation Guide

This guide documents validation in the supported rewrite contract.

> Audience: advanced users who need strict contract and invariant behavior.
> New users should begin with
> [Quickstart](getting-started/quickstart-first-workflow.md),
> [Troubleshooting: first-run and supported-lane failures](getting-started/troubleshooting-first-run.md), and
> [Core concepts](concepts/core-concepts.md).

For API shapes, see [`api.md`](api.md).

## Boundary Model

Validation is split across three boundaries:

- Builder request boundary:
  validates supported source types (`DataFrame` or file path) and builder conventions.
- Final dataset boundary:
  validates strict `AnalysisReadyPhosphoDataset` invariants.
- Workflow boundaries:
  validate request DTOs, config ranges, reference compatibility, and runtime/science seams.

## Signalome Protein-Identity Prerequisite

Supported signalome execution has an explicit protein-identity contract:

- signalome requires explicit, non-empty `dataset.site_metadata.protein_id`
  for every interpreted site
- gene-symbol site IDs (for example `"<gene_symbol>;<site>;"`) encode site
  identity, not protein identity, and are not a fallback substitute
- this is an intentional scientific boundary for protein-aware signalome
  grouping and module assignment
- builder flexibility at ingestion does not weaken this downstream workflow
  contract

## Builder Flexibility vs Dataset Strictness

`DatasetBuildRequest` accepts:

- `pandas.DataFrame`
- file path (`str`, `Path`, `PathLike`) to `.csv`, `.tsv`/`.txt`, `.parquet`
  (`.parquet` requires optional parquet dependencies; install with
  `pip install "phospy[parquet]"`)

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

Important contract distinction:

- Final dataset boundary: `site_metadata.site_sequence` is optional in
  `AnalysisReadyPhosphoDataset` (validated when present).
- Policy-specific preprocessing boundary:
  `site_matrix.policy="build_from_metadata"` still requires usable
  sequence-bearing rows for that path and excludes rows without usable sequence
  from that construction path.

`AnalysisReadyPhosphoDataset` itself is strict, missing-value-free, and
DataFrame-only. Workflows consume only this dataset type.

## AnalysisReady Dataset Validation

`AnalysisReadyPhosphoDataset` constructor composes:

- `AnalysisReadyDatasetValidator`
- `TransformationStateValidator`

Enforced dataset invariants:

- `phospho`: non-empty numeric DataFrame, missing-value-free, unique
  index/columns, canonical site IDs.
- `site_metadata`: non-empty DataFrame, exact index alignment with `phospho.index`,
  required columns `gene_symbol`, `site` with non-empty strings.
- site-identity coherence: each `phospho.index` row ID must be parseable as
  `"<gene_symbol>;<site>;"` and the parsed values must exactly match
  `site_metadata.gene_symbol` and `site_metadata.site` for the same row.
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
Direct caller-side establishment attempts are rejected at construction time:
`TransformationState.established_raw(...)` and
`establish_transformation_state(...)` only accept approved internal
establishment authority from supported lanes.

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
  - `site_matrix.duplicate_site_strategy` must be one of:
    `max_mean_signal`, `first`, `aggregate_mean`, `aggregate_median`, `error`
  - `site_matrix.missing_data_policy` must be `drop_any_missing` for strict
    `AnalysisReadyPhosphoDataset` construction in the public complete-case lane
  - `site_matrix.minimum_observed_values` is not supported in the public
    strict analysis-ready lane and must be `None`
  - when `site_matrix.policy='as_input'`, site-matrix execution-only fields
    (`duplicate_site_strategy`, `missing_data_policy`,
    `minimum_observed_values`) must remain unset/default
  - `comparisons.policy` must be one of `none`, `sample_metadata_pairs`
  - `comparisons.sample_group_column` must be a non-empty string
  - when `comparisons.policy='none'`, `comparisons.pairs` must be unset
  - when `comparisons.pairs` is provided:
    - pairs must be `(left_group, right_group)` tuples
    - pair group names must be non-empty
    - self pairs and reverse-direction duplicates are rejected
- at execution, `missing_data.min_observed_values` must not exceed phospho
  sample count
- at request validation, `total_protein_correction.policy='ratio_to_total'`
  requires `request.total` to be present
- at request validation, `comparisons.policy='sample_metadata_pairs'` requires
  `request.sample_metadata` to be present
- at execution, `total_protein_correction.policy='ratio_to_total'` requires:
  - `total.columns` exactly matching `phospho.columns`
  - numeric phospho/total columns
  - unique normalized `total.index` identifiers
  - complete matching between `site_metadata.gene_symbol` and `total.index`
- at execution, `site_matrix.policy='build_from_metadata'` requires:
  - `site_metadata` columns `gene_symbol`, `site`
  - row-level sequence support established from supplied
    `site_metadata.site_sequence` values and/or bundled derivation when available
    (bundled derivation keys each row from `gene_symbol` + `site`, with
    row-index fallback)
  - usable row-level `site_sequence` values for rows that should participate in
    site-matrix construction; rows lacking usable sequence are excluded from
    this path rather than auto-filled or inferred
  - non-empty values for all rows in `gene_symbol`/`site` and site tokens that
    normalize to canonical `SITE_TOKEN` form (for example `S123`)
  - deterministic site identity construction as canonical
    `GENE_SYMBOL;SITE_TOKEN;`
  - the supported public row-retention rule is
    `site_matrix.missing_data_policy='drop_any_missing'`, so rows with
    incomplete phospho values are excluded before final dataset construction
  - retained-missingness site-matrix behaviour remains internal-only and is
    rejected at public config validation as incompatible with strict
    `AnalysisReadyPhosphoDataset` construction
  - duplicate-site handling by explicit `site_matrix.duplicate_site_strategy`:
    - `max_mean_signal`, `first`, `aggregate_mean`, `aggregate_median`, or `error`
  - at least one retained row after sequence filtering, missing-data policy, and
    duplicate-site handling (otherwise a diagnostic-rich input error is raised)
  - effective row retention can be narrower than the original metadata table as
    a supported consequence of the selected preprocessing policy
- at execution, `comparisons.policy='sample_metadata_pairs'` requires:
  - `sample_metadata.index` exactly matching `phospho.columns`
  - `sample_metadata[comparisons.sample_group_column]` with non-empty values
  - explicit pairs (when provided) to reference only observed sample groups
  - inferred mode (no explicit pairs) to have at least two observed groups

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
  non-empty string values for every interpreted site
- gene-symbol-prefixed site IDs are site identity and are not accepted as
  protein-identity fallback
- missing values in the upstream downstream score matrix are allowed; this is a
  normal outcome for correlation-based kinase scoring in low-information rows
- infinite values in upstream score/prediction matrices remain hard failures

Why this contract is strict:

- signalome grouping/module assignment is protein-identity-aware
- deriving only `gene_symbol`/`site` at builder input does not establish protein
  identity
- builder flexibility at ingestion does not weaken the downstream supported
  signalome contract

`SignalomeWorkflowInterpreter` preconditions the downstream score lane before
execution:

- rows with no finite kinase score support (all-missing rows) are identified as
  unsupported evidence for score-driven network correlation inputs
- partially missing rows are retained and consumed with pairwise-complete
  correlation handling
- prediction rows remain available for module assignment logic
- drop policy is caller-owned via
  `SignalomeConfig.score_preconditioning_policy`:
  - `allow_and_report` (default): continue and surface dropped-row counts in
    `signalome_result.score_preconditioning_diagnostics`
  - `error_on_drop`: fail at the interpreter boundary when any all-missing row
    would be dropped, with explicit policy/count details in the boundary error

This policy is a scientific contract decision: it controls whether the effective
score matrix used for signalome clustering/network stages may be narrower than
the aligned upstream matrix.

Interpreters/executors enforce seam-level scientific/runtime boundary checks and raise
`WorkflowBoundaryError` with seam names, concrete counts, and `next_action` hints.

For the supported kinase scoring lane, motif sequence inputs come from
`references.site_sequences` (resolved reference bundle), not from
`dataset.site_metadata.site_sequence`.

## Kinase Scoring/Prediction Contract Classification (2026-04-22)

Validation and workflow execution intentionally enforce rewrite contract
behavior that differs from legacy defaults in the kinase scoring/prediction
path:

- `KinaseWorkflowExecutor` uses profile+motif combine with profile fallback
  enabled and preserves profile evidence when motif values are missing.
- Workflow prediction candidate filtering uses fixed execution semantics:
  `score_threshold=0.0`, `inclusion=1`, and caller-owned `top_k`.
- Public request/config validation accepts rewrite knobs only:
  `mode`, `adaptive_policy`, `top_k`, `ensemble_size`, `n_iterations`,
  `random_state`, and scoring `profile_missing_value_strategy`.
- Legacy prediction/scoring knobs are intentionally out of contract and should
  fail fast at request-object construction (for example `svm_mode`,
  `allow_profile_only_fallback`, `score_threshold`, `inclusion`,
  `profile_policy`).

Classification for these differences:

- Intentional and supported: all runtime differences listed above.
- Temporary and should be removed for parity: none in runtime behavior.
- Unresolved and requires design decision: none currently tracked in this lane.

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
| Kinase workflow request/config validity | `KinaseWorkflowValidator` + `KinaseWorkflowConfigValidator` |
| Signalome workflow request/config validity | `SignalomeWorkflowValidator` + `SignalomeConfigValidator` |
| Kinase/signalome runtime seam diagnostics | workflow interpreters/executors (`WorkflowBoundaryError`) |

## Quick Troubleshooting

Start with [Troubleshooting: first-run and supported-lane failures](getting-started/troubleshooting-first-run.md) if you want symptom-first guidance. Use the table below when you already know you need contract-level detail.

| Problem | Usually means | Good next step |
| --- | --- | --- |
| Builder rejects input format | Field is neither DataFrame nor supported file path | Pass DataFrame or path to `.csv`/`.tsv`/`.txt`, or install parquet extras before using `.parquet` (`pip install "phospy[parquet]"`) |
| Dataset constructor fails on site metadata | Required strict boundary columns/values are missing | Provide `gene_symbol` and `site` with non-blank strings |
| Dataset constructor fails on site-identity coherence | `phospho.index` site IDs disagree with row-level `site_metadata.gene_symbol` / `site_metadata.site` (or site IDs are not parseable as `"<gene_symbol>;<site>;"`) | Ensure each row ID and metadata row describe the same site; fix source data rather than mutating at runtime |
| Row count unexpectedly drops with `site_matrix.policy='build_from_metadata'` | Rows without usable `site_sequence` cannot participate in sequence-derived site-matrix construction and are excluded | Compare input row count vs `dataset.phospho.shape[0]`, review `site_metadata.site_sequence` completeness, and choose policy intentionally |
| `ReferencePreset.AUTO` fails | Dataset organism is missing | Set `organism` in `DatasetBuildRequest` |
| Bundled human/mouse preset fails | Bundled references are rat-only in this release | Provide explicit non-rat `ReferenceBundle` |
| Kinase boundary seam fails | Overlap/support constraints were not met | Read seam details and adjust dataset/references/config |
| Signalome validation fails on `site_metadata.protein_id` | Explicit protein identity required for supported signalome execution | Provide non-empty `site_metadata.protein_id` for all interpreted sites; gene-symbol site-ID prefixes are not a substitute |

## Where Next

- Public API contract details: [API Guide](api.md)
- User-level workflow navigation: [Workflow guides](workflow-guides/index.md)
- Governance-level science interpretation: [Parity to PhosR](parity.md)

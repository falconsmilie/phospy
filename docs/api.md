# API Guide

This guide describes the current supported public contract only.

## Supported Lanes

PhosPy has one public dataset boundary and two public workflow stories:

- Dataset construction:
  `DatasetBuildRequest -> AnalysisReadyDatasetBuilder.run(request) -> AnalysisReadyPhosphoDataset`
- Kinase workflow:
  `KinaseWorkflow.run(KinaseWorkflowRequest(...)) -> KinaseWorkflowResult`
- Signalome workflow:
  `SignalomeWorkflow.run(SignalomeWorkflowRequest(...)) -> SignalomeWorkflowResult`

All public executors use `run(request)`.

## Package Boundary

- `src/phospy/`: supported package
- `legacy_archive/phospy_legacy/`: historical migration reference only

## Installation Contract

Standard installation includes all dependencies required by supported prediction
lanes, including `mode="adaptive_ensemble"` (scikit-learn is part of base
dependencies). No extra install step is required for adaptive mode.

## Import Contract

`phospy.api` is the canonical namespace where public API types are defined and
organised in source.

Top-level `phospy` is the primary supported import route for user-facing code,
examples, and quickstarts.

This split is intentional:

- `phospy.api` owns API definition and package structure.
- top-level `phospy` is the stable curated facade users are expected to import
  from.
- detailed authored types that are not part of the curated facade remain under
  `phospy.api` (for example `phospy.api.results` stage result containers).

## Public Types

Import from top-level `phospy`.

- Dataset and references:
  `AnalysisReadyPhosphoDataset`, `Organism`, `ReferencePreset`, `ReferenceBundle`
- Builder:
  `DatasetBuildRequest`, `AnalysisReadyDatasetBuilder`
- Workflows and requests:
  `KinaseWorkflow`, `KinaseWorkflowRequest`,
  `SignalomeWorkflow`, `SignalomeWorkflowRequest`
- Config models:
  `DatasetPreprocessingConfig`, `DatasetMissingDataConfig`,
  `DatasetTotalProteinCorrectionConfig`, `DatasetSiteMatrixConfig`,
  `DatasetComparisonBuildingConfig`, `KinaseScoringConfig`,
  `KinasePredictionConfig`, `KinaseActivityConfig`, `SignalomeConfig`
- Result models:
  `KinaseWorkflowResult`, `SignalomeWorkflowResult`

Stage-level result containers (`KinaseScoringResult`, `KinasePredictionResult`,
`KinaseActivityResult`) remain public in the canonical authored namespace:
`phospy.api.results`.

## Builder Contract

There is one public builder story: `AnalysisReadyDatasetBuilder.run(DatasetBuildRequest(...))`.

`DatasetBuildRequest` supports both public input routes:

- pandas `DataFrame` values
- file paths (`str`, `pathlib.Path`, `os.PathLike`) for supported table formats

Required request fields:

- `phospho`
- `site_metadata`

Optional request fields:

- `sample_metadata`
- `total`
- `organism`
- `preprocessing_config` (`DatasetPreprocessingConfig`)

After loading, both routes share the same normalization and validation path.

Supported site-metadata alias mapping is explicit and narrow:

- `gene_symbol`: `gene_symbol`, `gene_name`
- `site`: `site`
- `site_sequence`: `site_sequence`, `centralized_sequence`
- `protein_id`: `protein_id`

Unsupported legacy aliases (`gene`, `residue`, `phosphosite`, `site_position`,
`sequence`, `protein`) are rejected instead of guessed.

If `gene_symbol` and/or `site` are absent, one derivation convention is supported:
`site_metadata.index` values formatted exactly as `"<gene_symbol>;<site>;"`.
This derivation populates only `gene_symbol` and `site`; it does not infer
`protein_id`.

Current transformation establishment policy in this public builder lane is
intentionally narrow:

- builder execution establishes only the supported pass-through `linear` state
- quantitative matrix values are preserved as provided after builder normalization,
  except when an explicit preprocessing policy modifies them
- no additional log or heuristic transformation path is exposed in the public builder

`DatasetPreprocessingConfig` is grouped and builder-owned:

- `missing_data` (`DatasetMissingDataConfig`)
- `total_protein_correction` (`DatasetTotalProteinCorrectionConfig`)
- `site_matrix` (`DatasetSiteMatrixConfig`)
- `comparisons` (`DatasetComparisonBuildingConfig`)

Current supported policies:

- `missing_data.policy="forbid"` (default): no missing-value preprocessing.
- `missing_data.policy="impute_row_median"`: drop rows below
  `missing_data.min_observed_values`, then row-median imputation.
- `total_protein_correction.policy="none"` (default): no total/protein correction.
- `total_protein_correction.policy="ratio_to_total"`: subtract matched
  `total` abundance from `phospho` abundance per sample column in builder
  preprocessing. This policy requires `total` input, exact
  `total.columns == phospho.columns`, and complete
  `site_metadata.gene_symbol` vs `total.index` matching.
- `site_matrix.policy="as_input"` (default): preserve interpreted phospho/site rows.
- `site_matrix.policy="build_from_metadata"`: construct site-matrix-ready rows
  from `site_metadata.gene_symbol`, `site_metadata.site`, and
  `site_metadata.site_sequence` after upstream missing-data and
  total/protein-correction stages.
  This path is sequence-dependent at preprocessing time: rows must have usable
  `site_sequence` values to participate in construction. Rows lacking usable
  sequence are excluded from this path, so retained rows can be narrower than
  the original metadata table.
  This requirement is specific to the selected preprocessing policy and does
  not change the final dataset boundary where `site_sequence` remains optional.
  Additional supported policy controls under `DatasetSiteMatrixConfig`:
  - `missing_data_policy="drop_any_missing"` (default): keep only complete rows.
  - strict analysis-ready boundary alignment: retained-missingness site-matrix
    modes are not supported in the public builder lane.
  - `minimum_observed_values` is internal-only compatibility state and must
    remain unset in the supported public builder lane.
  - `duplicate_site_strategy="max_mean_signal"` (default): keep strongest row.
  - `duplicate_site_strategy="first"`: keep first duplicate row.
  - `duplicate_site_strategy="aggregate_mean"`: aggregate duplicate rows by mean.
  - `duplicate_site_strategy="aggregate_median"`: aggregate duplicate rows by median.
  - `duplicate_site_strategy="error"`: fail on duplicate constructed site IDs.
  Site IDs are constructed deterministically as canonical
  `GENE_SYMBOL;SITE_TOKEN;` identifiers.
  When `site_matrix.policy="as_input"`, these execution-only fields must remain
  at defaults and are rejected if overridden.
- `comparisons.policy="none"` (default): do not construct comparison columns.
- `comparisons.policy="sample_metadata_pairs"`: build dataset-level pairwise
  comparison columns in `dataset.comparisons` from grouped sample metadata.
  Required and supported inputs:
  - `sample_metadata` must be provided
  - `sample_metadata.index` must align to `phospho.columns`
  - `sample_metadata[comparisons.sample_group_column]` must contain one
    non-empty group label per sample
  - `comparisons.pairs` can pass explicit `(left_group, right_group)` pairs;
    when omitted, all unique pairwise combinations are inferred from observed
    groups.

## Final Dataset Boundary

`AnalysisReadyPhosphoDataset` is strict and workflow-facing.

- It owns validated tables, not input files.
- It requires DataFrame values for `phospho` and `site_metadata` at construction time.
- `site_metadata` must contain `gene_symbol`, `site` with non-empty strings.
- Site identity coherence is strict: each `phospho.index` row ID must be
  interpretable as `"<gene_symbol>;<site>;"`, and parsed values must exactly
  match the corresponding `site_metadata.gene_symbol` and `site_metadata.site`
  row values.
- `site_sequence` is optional at this boundary; when present it must be non-empty.
- Final-boundary optionality does not remove policy-specific preprocessing
  requirements: `site_matrix.policy="build_from_metadata"` still requires
  usable sequence-bearing rows for that construction path.
- Site identifiers must already be canonical and non-colliding.
- `sample_metadata` (if present) must align to `phospho.columns`.
- `total` (if present) must be numeric and column-aligned to `phospho`.
- `comparisons` (if present) must be numeric, non-missing, and aligned to
  `phospho.index`.
- `transformation_state` is mandatory and must be both coherent and
  established through a supported PhosPy path.

Supported establishment paths are:

- `AnalysisReadyDatasetBuilder.run(...)` (default supported lane)
- supported transformer execution through the dataset builder executor/resolver
- supported bundle reconstruction paths

For `AnalysisReadyDatasetBuilder.run(...)`, the established state is currently
the pass-through `linear` path.
When `total_protein_correction.policy="ratio_to_total"` is used, corrected
quantitative values are represented directly in `dataset.phospho`.

Directly declared transformation-state objects are not accepted as authoritative
at the dataset boundary unless they were established through one of the supported
PhosPy paths above.
Direct caller-side minting is explicitly rejected: `TransformationState.established_raw(...)`
and `establish_transformation_state(...)` are accepted only when an approved
internal establishment authority is provided by the supported
builder/transformer or bundle-reconstruction lanes.

Builder flexibility does not weaken this final dataset strictness.
Workflows consume only `AnalysisReadyPhosphoDataset`.

## Workflow Contract

`KinaseWorkflowRequest` fields:

- `dataset: AnalysisReadyPhosphoDataset`
- `references: ReferencePreset | ReferenceBundle`
- `scoring_config: KinaseScoringConfig`
- `prediction_config: KinasePredictionConfig`
- `activity_config: KinaseActivityConfig | None`

`SignalomeWorkflowRequest` fields:

- `kinase_result: KinaseWorkflowResult`
- `config: SignalomeConfig`
- supported prerequisite: `kinase_result.dataset.site_metadata.protein_id` must be
  present and non-empty for all interpreted sites

`KinaseScoringConfig` fields:

- `min_substrates` (validated floor: `>= 2`)
- `include_diagnostic_scoring_tables`
- `profile_missing_value_strategy` (`"strict"` or `"median_skipna"`)

`KinasePredictionConfig` fields:

- `top_k` (validated floor: `>= 1`)
- `ensemble_size` (validated floor: `>= 1`)
- `mode` (`"deterministic_ranking"` or `"adaptive_ensemble"`)
- `adaptive_policy` (`"stable"` or `"r_parity"`)
- `n_iterations` (validated floor: `>= 1`; used by adaptive lane)
- `random_state` (`None` or integer `>= 0`)

`"adaptive_ensemble"` is part of the normal supported contract and is expected
to work after a standard package install.

Kinase scoring authority in the supported lane:

- Scoring is a function of:
  - `request.dataset` (analysis-ready phospho matrix)
  - resolved `ReferenceBundle` content (`kinase_substrate_map`, `site_sequences`)
  - `request.scoring_config`
- Scoring is not changed by:
  - `prediction_config.mode` (`deterministic_ranking` vs `adaptive_ensemble`)
  - whether references were supplied via `ReferencePreset` or explicit
    `ReferenceBundle` when both resolve to equivalent reference content
- Adaptive prediction is downstream of scoring and consumes the same authoritative
  scoring outputs as deterministic prediction.
- Adaptive-policy settings (`adaptive_policy`) are prediction-stage controls and
  do not opt into a separate scoring contract.

`SignalomeConfig` fields:

- `substrate_support_cutoff`
- `network_correlation_threshold`
- `network_policy` (`"positive_only"`, `"absolute_threshold"`, or `"signed"`)
- `assignment_policy` (`"cutoff_binary"` or `"weighted_top"`)
- `score_preconditioning_policy` (`"allow_and_report"` default, or `"error_on_drop"`)
- `module_count` (`None` for automatic module-count selection)
- `module_selection_primary_correlation_threshold`
- `module_selection_fallback_correlation_threshold`
- `module_selection_max_clusters`

## Result Contract (Nested Stage Outputs)

`KinaseWorkflowResult`:

- `result.dataset`
- `result.references`
- `result.scoring_result.profile_scores`
- `result.scoring_result.combined_scores`
- `result.scoring_result.motif_scores` (optional diagnostic field)
- `result.scoring_result.weights` (optional diagnostic field)
- `result.prediction_result.pred_mat`
- `result.prediction_result.substrate_list` (optional)
- `result.activity_result` (`None` when activity is disabled)

`SignalomeWorkflowResult`:

- `result.dataset`
- `result.kinase_result` (full upstream nested lineage)
- `result.module_assignments.table`
- `result.signalome_modules.table`
- `result.kinase_network.edges`
- `result.kinase_network.nodes` (optional)
- `result.module_selection_diagnostics`
- `result.score_preconditioning_diagnostics`
- `result.expanded_signalome` (official supported output; optional by type for
  compatibility, populated in the supported executor lane)

`module_selection_diagnostics` fields:

- `strategy` (`"correlation_thresholds"` or `"explicit_module_count"`)
- `selected_module_count`
- `requested_module_count`
- `threshold_used`
- `max_clusters_evaluated`
- `candidate_scores` (per-candidate `min_median_correlation` and `mean_median_correlation`)
- `reason`
- `zero_variance_profile_count`
- `near_constant_profile_count`
- `excluded_from_correlation_count`

`score_preconditioning_diagnostics` fields:

- `policy` (`"allow_and_report"` or `"error_on_drop"` from
  `SignalomeConfig.score_preconditioning_policy`)
- `input_row_count` (aligned downstream-score rows before preconditioning)
- `dropped_all_missing_row_count` (rows removed because all kinase scores are missing)
- `retained_row_count` (rows retained for signalome score-driven stages)

`module_assignments.table` includes site-level assignment metadata:

- structural fields: `protein_id`, `module_id`, `top_kinase`, `top_score`
- tie/ambiguity fields: `top_kinase_candidates`, `top_kinase_weights`,
  `top_kinase_tie_count`, `top_kinase_is_ambiguous`,
  `top_kinase_selection_policy`
- module-level attribution fields: `module_top_kinase`,
  `module_top_kinase_candidates`, `module_top_kinase_tie_count`,
  `module_top_kinase_is_ambiguous`, `module_top_kinase_selection_policy`

In the supported signalome lane, `protein_id` is always taken from validated
`dataset.site_metadata.protein_id` (required and non-empty for interpreted
sites). Gene-symbol-prefixed site IDs are not used as protein-identity
fallbacks.

`assignment_policy` behavior:

- `"cutoff_binary"`: support is binary from `prediction_result.pred_mat >
  substrate_support_cutoff`
- `"weighted_top"`: support is fractional from per-site `top_kinase_weights`
  tie metadata propagated through module/expanded outputs

`expanded_signalome` flattened schema:

- `kinase`: focal kinase
- `row_kind`: `"site"` or `"summary"`
- `assignment_policy`: emitted assignment policy
- `linked_kinases`: JSON array string of focal + linked kinases
- `regulated_module_ids`: JSON array string of module IDs where focal-kinase
  share is strictly greater than `1.0` percent
- `site_id`: selected phosphosite ID (`""` on summary rows)
- `site_order`: zero-based position in the original module-assignment order
- `protein_id`, `module_id`, `top_kinase`, `top_score`: selected site metadata
- `support_kinases`: JSON array string of linked kinases supporting the site
- `support_weight`: site support weight under the selected assignment policy

`expanded_signalome` row-selection conditions:

- The executor emits one focal-kinase block per kinase in `signalome_modules.columns`.
- A `row_kind="site"` row is emitted when:
  - the site belongs to a regulated module for that focal kinase
    (`signalome_modules[module_id, kinase] > 1.0`), and
  - at least one linked kinase provides positive site support under the active
    `assignment_policy`.
- If no site rows qualify for a focal kinase, one `row_kind="summary"` row is
  emitted for that kinase to preserve linked-kinase and regulated-module metadata.

No top-level convenience mirrors flatten nested stage outputs.

## Supported Science vs Deferred Science

Supported public lane (stable and recommended):

- Kinase scoring stage always outputs `profile_scores` and `combined_scores`.
- Diagnostic scoring tables (`motif_scores`, `weights`) are opt-in via
  `KinaseScoringConfig(include_diagnostic_scoring_tables=True)`.
- Motif scoring in the supported kinase lane uses `references.site_sequences`
  from the resolved `ReferenceBundle`.
- Prediction stage uses a downstream score matrix that resolves to
  `combined_scores` when present and falls back to `profile_scores` only when
  combined scores are unavailable.
- Activity stage is supported and optional inside `KinaseWorkflow`.
- Signalome stage consumes the same downstream score-matrix lane as prediction
  and outputs module assignments, module matrix, kinase network, and
  `expanded_signalome`.
- Signalome assignment policy is explicit:
  `assignment_policy="cutoff_binary"` keeps cutoff/binary support semantics;
  `assignment_policy="weighted_top"` propagates fractional
  `top_kinase_weights` support into module shares.
- Signalome network policy is explicit:
  `network_policy="positive_only"` keeps only positive correlations above
  threshold; `network_policy="absolute_threshold"` keeps absolute correlations
  above threshold and emits unsigned edge correlations; `network_policy="signed"`
  keeps absolute correlations above threshold and emits signed edge
  correlations.
- Downstream score missingness is part of the supported scientific contract:
  all-missing score rows are explicitly policy-governed preconditioning
  candidates, partially missing rows are retained, and infinite values remain
  invalid.
- Preconditioning drop policy is caller-owned in the supported lane:
  `score_preconditioning_policy="allow_and_report"` drops all-missing rows and
  reports counts via `result.score_preconditioning_diagnostics`;
  `score_preconditioning_policy="error_on_drop"` fails interpretation when any
  all-missing row would be dropped.
- Active parity regression for this lane is rewrite execution against committed
  rewrite-owned fixtures; no live `legacy_archive` module execution is part of
  the active parity gate.

Deferred/experimental/not yet ported into the public lane:
- Additional legacy science lanes listed as roadmap follow-ons.

## Reference Resolution

- `ReferencePreset.AUTO` requires `dataset.organism`.
- Preset/dataset organism compatibility is enforced.
- Bundled runtime references are currently rat-only.
- `ReferencePreset.HUMAN` and `ReferencePreset.MOUSE` remain public enum lanes,
  but bundled resolution for those presets is intentionally unsupported in this release.
- Non-rat execution uses explicit caller-provided `ReferenceBundle`.

## User-Handleable Exceptions

Top-level `phospy` exports a focused, user-handleable error facade:

- `PhosPyError`
- `PhosPyInputError`, `UnsupportedInputFormatError`
- `PhosPyBuildError`
- `PhosPyValidationError`
- `PhosPyReferenceError`, `UnsupportedOrganismError`
- `PhosPyTransformationError`
- `PhosPyWorkflowError`, `WorkflowBoundaryError`

The complete public exception taxonomy remains available under
`phospy.errors`.

## Quick Usage Pattern

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    Organism,
    ReferencePreset,
)

dataset = AnalysisReadyDatasetBuilder().run(
    DatasetBuildRequest(
        phospho="./input/phospho.csv",
        site_metadata="./input/site_metadata.csv",
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            missing_data=DatasetMissingDataConfig(policy="forbid"),
        ),
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)

pred_mat = kinase_result.prediction_result.pred_mat
if kinase_result.activity_result is not None:
    weighted_activity = kinase_result.activity_result.weighted_activity
```

If you choose `site_matrix.policy="build_from_metadata"`, inspect row-retention
counts after builder execution (`dataset.phospho.shape[0]` versus input row
count). Rows without usable `site_sequence` do not participate in that
construction path and are excluded there. See
[`examples/dataset_builder_demo.py`](../examples/dataset_builder_demo.py) for a
concrete retained-vs-excluded example.

For CLI and bundle persistence details, see:

- [`cli.md`](cli.md)
- [`output_bundles.md`](output_bundles.md)

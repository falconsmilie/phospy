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
  `KinaseScoringConfig`, `KinasePredictionConfig`, `KinaseActivityConfig`,
  `SignalomeConfig`
- Result models:
  `KinaseWorkflowResult`, `SignalomeWorkflowResult`,
  `KinaseScoringResult`, `KinasePredictionResult`, `KinaseActivityResult`

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

## Final Dataset Boundary

`AnalysisReadyPhosphoDataset` is strict and workflow-facing.

- It owns validated tables, not input files.
- It requires DataFrame values for `phospho` and `site_metadata` at construction time.
- `site_metadata` must contain `gene_symbol`, `site` with non-empty strings.
- `site_sequence` is optional at this boundary; when present it must be non-empty.
- Site identifiers must already be canonical and non-colliding.
- `sample_metadata` (if present) must align to `phospho.columns`.
- `total` (if present) must be numeric and column-aligned to `phospho`.
- `transformation_state` is mandatory on direct dataset construction and must be coherent.

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

`SignalomeConfig` policy fields:

- `substrate_support_cutoff`
- `network_correlation_threshold`
- `assignment_policy` (`"cutoff_binary"` or `"weighted_top"`)

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
- `result.expanded_signalome` (optional DataFrame; populated in the supported lane)

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

No top-level convenience mirrors flatten nested stage outputs.

## Supported Science vs Deferred Science

Supported public lane today:

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
- Downstream score missingness is part of the supported scientific contract:
  all-missing score rows are preconditioned out of score-driven network inputs,
  partially missing rows are retained, and infinite values remain invalid.

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

Top-level `phospy` exports the public exception taxonomy:

- Base:
  `PhosPyError`
- Input/build:
  `PhosPyInputError`, `UnsupportedInputFormatError`, `PhosPyBuildError`,
  `DatasetBuildError`
- Validation:
  `PhosPyValidationError`, `DatasetValidationError`, `ReferenceValidationError`,
  `TransformationValidationError`, `WorkflowValidationError`
- Reference:
  `PhosPyReferenceError`, `ReferenceResolutionError`,
  `ReferenceCompatibilityError`, `UnsupportedOrganismError`
- Transformation:
  `PhosPyTransformationError`, `InvalidTransformationStateError`,
  `TransformationStateEstablishmentError`, `TransformerExecutionError`
- Workflow:
  `PhosPyWorkflowError`, `WorkflowBoundaryError`, `WorkflowStageError`

## Quick Usage Pattern

```python
from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
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
    )
)

kinase_result = KinaseWorkflow().run(
    KinaseWorkflowRequest(dataset=dataset, references=ReferencePreset.AUTO)
)

pred_mat = kinase_result.prediction_result.pred_mat
if kinase_result.activity_result is not None:
    weighted_activity = kinase_result.activity_result.weighted_activity
```

For CLI and bundle persistence details, see:

- [`cli.md`](cli.md)
- [`output_bundles.md`](output_bundles.md)

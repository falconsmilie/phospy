# API Guide

This guide describes the current supported rewrite contract only.

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

- `src/phospy/`: supported rewrite package
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
- `site_metadata` must contain `gene_symbol`, `site`, `site_sequence` with non-empty strings.
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

## Result Contract (Nested Stage Outputs)

`KinaseWorkflowResult`:

- `result.dataset`
- `result.references`
- `result.scoring_result.profile_scores`
- `result.scoring_result.motif_scores` (optional field in model)
- `result.scoring_result.combined_scores` (optional field in model)
- `result.scoring_result.weights` (optional field in model)
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
- `result.expanded_signalome` (optional; currently `None` in the default lane)

No top-level convenience mirrors flatten nested stage outputs.

## Supported Science vs Deferred Science

Supported public lane today:

- Kinase scoring stage outputs profile, motif, and combined scoring tables.
- Prediction stage remains profile-driven for ranking and prediction matrix assembly.
- Activity stage is supported and optional inside `KinaseWorkflow`.
- Signalome stage outputs module assignments, module matrix, and kinase network.

Deferred/experimental/not yet ported into the public lane:

- Population of `SignalomeWorkflowResult.expanded_signalome`.
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
